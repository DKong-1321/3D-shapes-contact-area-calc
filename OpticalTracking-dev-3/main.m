%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             READ ME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program runs every single specimen found in the root folder.
% File structure expected is SpecimenName/State/Loading condition:
% Experiment_3             <====== This is the folder you pick on the popup!
% ├───SN06_XXXX_right
% │   ├───Intact
% │   |    ├───Anterior.csv
% │   |    ├───Neutral.cs
% │   |    └───Sps.csv                   
% │   └───Digitisation                   
% │        ├───tibia_medial.csv          
% |        ...                           
% │        ├───pat med.csv             Notice SN08 appears twice. That is OK
% │        └───FL.csv                They are matched based on the first word (SN08)
% ├───SN08_XXXX_RK            <============╣  
% │   ├───ACL_cut                          ║  
% │   ├───digitise                         ║  
% ├───SN08_yyyy_2_RK          <============╝
% │   ├───ACL_recon
% │   └───Digitise_ACL_recon
% ...

% Make sure each folder is called SpecimenName_whatever_(side)_whatever
% e.g. SN06_right_122024, JJ01 left, XE0_2_left, MA03_333_lk, etc.
% right, left, lk, rk are the acceptable words.
% These can be changed in ./defaults.m`
% 
% lib/configure/clean_specimen_condition.m is a list of regex to help
% clean up inconsistency in your folders. e.g., between cases (ACL vs acl),
% in naming (Reconstruction vs repair vs recon), or mistakes (esp instead of Sps)

% clc; clear; close all;
%% Load default configuration. Check ./lib/configure/defaults.m if you want to modify them.

config = defaults();

disp("Choose the root folder where all the specimens are")
root = uigetdir(".", "Choose the root folder");

%%
module = Module.Knee;
digitisation = Digitisation.new(root, config, module);
jcs = JCS.new(digitisation);
kinematics = jcs.solve();
%%
disp("Loading tension")
kinematics.load_tension("**/*tension.csv");
disp("Done loading data");
%%
% jcs.print_to_file();
% jcs.plot()
kinematics.trajectories.intraspecimen_mean();
kinematics.trajectories.set_flexion_min(-3);
kinematics.trajectories.cp_sensor_to_kinematics();
[flex, ext] = kinematics.trajectories.path.split_flex_ext;

% flex.spm;

normalised = flex.normalise("Neutral", "Intact", "ACL_recon");
spm = normalised.spm();
dunnett = spm.inference(0.05).dunnett(normalised, "ACL_recon");
spm.between_subject.tibiofemoral.anterior.inference(0.05).plot('plot_threshold_label',true, 'plot_p_values',true, 'autoset_ylim',true);

norm_avg = normalised.average();
norm_avg.plot();


% toc


% %% Visualise stl
% neutral = strcmpi({transforms.Intact.name}, "neutral");
% visualise_stl(config, transforms.Intact(neutral).in_femur.tibia, transforms.Intact(neutral).in_femur.femur);
% %% print to file
% print_mean_std_to_file(stats, states, root);
% %% Plot
% truncate_min = -5;
% truncate_max = 90;
% 
% % plot_interspecimen(config, stats, states, truncate_min, truncate_max)
% plot_interspecimen(config, stats, stats_offset, states, truncate_min, truncate_max)
% 
% diary off;
