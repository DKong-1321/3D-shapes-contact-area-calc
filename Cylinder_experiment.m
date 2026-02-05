% contactArea_cylinderExperiment_HYBRID_anim.m
clear; clc; close all;

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = scriptDir;
addpath(genpath(fullfile(projectRoot,'src')));

modelDir = fullfile(projectRoot,'model','cylinder');
dataDir  = fullfile(projectRoot,'data');

stlFemur = fullfile(modelDir,'top-femur.stl');
stlTibia = fullfile(modelDir,'bottom_natural.stl');
matPath  = fullfile(dataDir,'cylinders.mat');

assert(exist(stlFemur,'file')==2, 'Missing: %s', stlFemur);
assert(exist(stlTibia,'file')==2, 'Missing: %s', stlTibia);
assert(exist(matPath,'file')==2,  'Missing: %s', matPath);

tol = 0.20;

wanted = [];              % e.g. ["insertion","rotation"]
animate = true;
animateEvery = 1;         % update every N frames
fps = 30;                 % playback speed (approx)
makeVideo = false;        % true -> saves mp4
videoPath = fullfile(projectRoot,'cylinders_contact_anim.mp4');

S = load(matPath);
trajectories = [];
if isfield(S,'trajectories'); trajectories = S.trajectories; end
if isempty(trajectories) && isfield(S,'Trajectories'); trajectories = S.Trajectories; end
if isempty(trajectories)
    fns = fieldnames(S);
    trajectories = S.(fns{1});
end
if istable(trajectories), trajectories = table2array(trajectories); end
if iscell(trajectories), trajectories = [trajectories{:}]; end
assert(numel(trajectories)>=1, 'No trajectories found.');

[Ff0,Vf0] = loadMeshAny(stlFemur);
[Ft0,Vt0] = loadMeshAny(stlTibia);

tibiaBody0 = buildBodyStruct(Ft0, Vt0);
tibiaBody0.faceNormals = faceNormalsFromFV(tibiaBody0.F, tibiaBody0.V);
tibiaBody0.kdtreeFaces = KDTreeSearcher(tibiaBody0.centroids);

Results = struct('name',{},'A_mm2',{},'Apen_mm2',{},'NcontactTris',{},'NpenTris',{},'tTf',{});

if makeVideo
    vw = VideoWriter(videoPath,'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
end

for it = 1:numel(trajectories)
    tr = trajectories(it);

    lc = "";
    if isfield(tr,'LoadingCondition'); lc = string(tr.LoadingCondition); end
    if lc=="", lc = "traj_"+it; end
    if ~isempty(wanted) && ~any(strcmpi(lc, wanted)), continue; end

    T = tr.Transform;
    hasG = isfield(T,'gTfi') && isfield(T,'gTti') && ~isempty(T.gTfi) && ~isempty(T.gTti);
    hasT = isfield(T,'fTt')  && ~isempty(T.fTt);

    if hasG
        gTf = T.gTfi;
        gTt = T.gTti;
        nFrames = size(gTf,3);
        assert(size(gTt,3)==nFrames, 'Frame mismatch in %s', lc);

        tTf = zeros(4,4,nFrames);
        for k = 1:nFrames
            tTf(:,:,k) = invSE3(gTt(:,:,k)) * gTf(:,:,k);
        end
    elseif hasT
        fTt = T.fTt;
        nFrames = size(fTt,3);
        tTf = zeros(4,4,nFrames);
        for k = 1:nFrames
            tTf(:,:,k) = invSE3(fTt(:,:,k));
        end
    else
        error('No usable transforms in %s', lc);
    end

    A    = nan(nFrames,1);
    Apen = nan(nFrames,1);
    Nct  = zeros(nFrames,1);
    Npn  = zeros(nFrames,1);

    if animate
        fig = figure('Color','w','Name',char(lc)); clf;
        ax = axes(fig); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
        xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
        view(ax,3); camlight(ax); lighting(ax,'gouraud');

        patch('Faces',Ft0,'Vertices',Vt0, ...
            'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.10,'EdgeColor','none','Parent',ax);

        pF = patch('Faces',Ff0,'Vertices',Vf0, ...
            'FaceColor',[0 0.7 0],'FaceAlpha',0.12,'EdgeColor',[0 0 0], ...
            'EdgeAlpha',0.10,'LineWidth',0.25,'Parent',ax);

        pC = patch('Faces',zeros(0,3),'Vertices',Vf0, ...
            'FaceColor',[1 0 1],'FaceAlpha',0.95,'EdgeColor','none','Parent',ax);
    end

    for k = 1:nFrames
        Vk = applySE3(Vf0, tTf(:,:,k));
        femurBody = buildBodyStruct(Ff0, Vk);

        res = computeContactArea_STS_hybrid(tibiaBody0, femurBody, tol, struct());

        A(k)    = res.contactArea;
        Apen(k) = res.penetrationArea;
        Nct(k)  = nnz(res.contactMask);
        Npn(k)  = nnz(res.penMask);

        if animate && mod(k,animateEvery)==0
            set(pF,'Vertices',Vk);

            m = res.contactMask;
            if isempty(m) || ~islogical(m) || numel(m)~=size(Ff0,1), m = false(size(Ff0,1),1); end
            if any(m)
                set(pC,'Vertices',Vk,'Faces',Ff0(m,:));
            else
                set(pC,'Faces',zeros(0,3));
            end

            title(ax, sprintf('%s | frame %d/%d | A=%.3f mm^2 | Apen=%.3f mm^2', ...
                lc, k, nFrames, A(k), Apen(k)), 'Interpreter','none');

            drawnow;

            if makeVideo
                writeVideo(vw, getframe(fig));
            else
                pause(1/fps);
            end
        end
    end

    Results(end+1) = struct('name',lc,'A_mm2',A,'Apen_mm2',Apen, ...
        'NcontactTris',Nct,'NpenTris',Npn,'tTf',tTf); %#ok<SAGROW>

    figure('Color','w');
    plot(A,'-'); hold on; plot(Apen,'-');
    grid on; xlabel('Frame'); ylabel('Area (mm^2)');
    legend('Hybrid (d<=tol)','Penetration (d<=0)','Location','best');
    title(char(lc),'Interpreter','none');
end

if makeVideo
    close(vw);
end

save(fullfile(projectRoot,'results_contactArea_cylinders_HYBRID.mat'), ...
    'Results','tol','stlFemur','stlTibia','matPath');

%% ===== funcs =====
function [F,V] = loadMeshAny(path)
if exist('loadStlMesh','file')
    [F,V] = loadStlMesh(path);
elseif exist('loadAnyStl','file')
    [F,V] = loadAnyStl(path);
else
    error('Need loadStlMesh or loadAnyStl on path.');
end
end

function V2 = applySE3(V, T)
Vh  = [V, ones(size(V,1),1)];
V2h = (T * Vh.').';
V2  = V2h(:,1:3);
end

function Ti = invSE3(T)
R = T(1:3,1:3);
p = T(1:3,4);
Ti = eye(4);
Ti(1:3,1:3) = R.';
Ti(1:3,4)   = -R.'*p;
end

function N = faceNormalsFromFV(F,V)
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);
N = cross(v2 - v1, v3 - v1, 2);
nrm = sqrt(sum(N.^2,2));
nrm(nrm==0) = 1;
N = N ./ nrm;
end
