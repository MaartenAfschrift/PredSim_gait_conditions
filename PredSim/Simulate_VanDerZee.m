%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear; close all; clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);

%% Initialize S
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings);
[S] = initializeSettings();
S.misc.main_path = pathRepo;
addpath([S.misc.main_path '\VariousFunctions'])

%% Required inputs

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
S.subject.IG_selection_gaitCyclePercent = 50;

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

% % S.bounds
S.bounds.a.lower            = 0.01;
S.solver.CasADi_path        = get_casadi_path();
S.subject.mtp_type          = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};

% %S.Cpp2Dll: required inputs to convert .osim to .dll
S.Cpp2Dll.PathCpp2Dll_Exe = fullfile(pathRepo,'Osim2DLL');
S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';

S.solver.N_threads      = 2;
S.solver.N_meshes       = 50;
S.solver.par_cluster_name = 'Cores4';

% model info
S.subject.name = 'Falisse_et_al_2022';
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% copy of the default settings
SDefault = S;

%% Different trials


% our model represents subject 4 approximately, so
s_star = 0.68;
f_star = 1.83./2; % divided by 2 because we implemented stride frequency

S = SDefault;
S = SetVanDerZee_Settings(S, [0.56]*f_star, 0.7, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial1');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [], 0.7, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial2');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 1*f_star, 0.7, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial3');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 0.72*f_star, 0.9, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial4');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, []*f_star, 0.9, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial5');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 1*f_star, 0.9, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial6');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 0.88*f_star, 1.1, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial7');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [], 1.1, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial8');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 1*f_star, 1.1, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial9');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, 1*f_star, 1.6, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial10');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, []*f_star, 1.6, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial11');
add_pred_sim_to_batch(S,osim_path);


S = SDefault;
S = SetVanDerZee_Settings(S, [1.28]*f_star, 1.6, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial12');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.8, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial13');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, []*f_star, 1.8, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial14');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1.44]*f_star, 1.8, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial15');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, []*f_star, 2.0, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial16');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [0.7]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial17');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [0.8]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial18');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [0.9]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial19');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial20');
add_pred_sim_to_batch(S,osim_path);


S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial21');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial22');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1.1]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial23');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1.2]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial24');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1.3]*f_star, 1.25, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial25');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, 0);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial26');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, 0.1);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial27');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, 0.2);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial28');
add_pred_sim_to_batch(S,osim_path);


S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, 0.3);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial29');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.25, 0.4);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial30');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1.12]*f_star, 1.4, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial31');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, []*f_star, 1.4, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial32');
add_pred_sim_to_batch(S,osim_path);

S = SDefault;
S = SetVanDerZee_Settings(S, [1]*f_star, 1.4, []);
S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','VanDerZee','trial33');
add_pred_sim_to_batch(S,osim_path);
