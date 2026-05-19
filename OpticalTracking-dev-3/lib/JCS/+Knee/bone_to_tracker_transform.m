function transforms = bone_to_tracker_transform(trackers, config)

tibia = trackers.tibia;
femur = trackers.femur;
patella = trackers.patella;
right = config.is_right_knee;

%% Define coordinate systems for each bone using digitised points
% For tibia
gTt0 = defineBodyFixedFrameTibia(tibia, right);
grt0 = origins(gTt0);
gTf0 = defineBodyFixedFrameFemur(femur, right);
grf0 = origins(gTf0);
gTp0 = defineBodyFixedFramePatella(patella, right);
grp0 = origins(gTp0);

%% Define frame of reference for each of the trackers in global coordinates
gTft0 = femur.tracker.map(@(x) defineTrackerFixedFrame(x.rotations_mean, x.translations_mean));
gTtt0 = tibia.tracker.map(@(x) defineTrackerFixedFrame(x.rotations_mean, x.translations_mean));
gTpt0 = patella.tracker.map(@(x) defineTrackerFixedFrame(x.rotations_mean, x.translations_mean));

gTft0 = gTft0.unwrap_or([]);
gTtt0 = gTtt0.unwrap_or([]);
gTpt0 = gTpt0.unwrap_or([]);
%% Relate body fixed frames and origin to the tracker rigid body
% A constant transform of the body fixed frame in the tracker frame of reference (assumes rigid body)
% As in, tibia in the tibial tracker's frame of reference.
transforms.tibia.transform = gTtt0\gTt0; %ttTtc
transforms.tibia.origin = gTtt0\grt0; %ttrtc

transforms.femur.transform = gTft0\gTf0; %ftTfc
transforms.femur.origin = gTft0\grf0; %ftrfc

transforms.patella.transform = gTpt0\gTp0; %ttTpc
transforms.patella.origin = gTpt0\grp0; %ptrpc

% trackers.tibia.transform = gTtt0.map(@(x) x\gTt0); %ttTtc
% trackers.tibia.origin = gTtt0.map(@(x) x\grt0); %ttrtc
%
% trackers.femur.transform = gTft0.map(@(x) x\gTf0); %ftTfc
% trackers.femur.origin = gTft0.map(@(x) x\grf0); %ftrfc
%
% trackers.patella.transform = gTpt0.map(@(x) x\gTp0); %ttTpc
% trackers.patella.origin = gTpt0.map(@(x) x\grp0); %ptrpc

end

function o = origins(tracker)
    if isempty(tracker)
        o = [];
    else
        o = tracker(:,4);
    end
end
