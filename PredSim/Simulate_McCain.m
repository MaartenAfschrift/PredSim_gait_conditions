%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear
close all
clc
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
S.subject.IG_selection = 'quasi-random';

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
S.bounds.a.lower            = 0.01;
S.misc.gaitmotion_type      = 'FullGaitCycle';
S.solver.CasADi_path        = get_casadi_path();
S.subject.mtp_type          = '2022paper';
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};

% %S.Cpp2Dll: required inputs to convert .osim to .dll
S.Cpp2Dll.PathCpp2Dll_Exe = fullfile(pathRepo,'Osim2DLL');
S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';


S.solver.N_threads      = 1;
S.solver.N_meshes       = 75;
S.solver.par_cluster_name = 'Cores3';

% walking speed
S.subject.v_pelvis_x_trgt   = 0.8;

% copy of settings structure
SDefault = S;


%% Loop over different stiffness values

kVect = 0:50:150;
for i= 1:length(kVect)
    kSel =kVect(i);
    kSel_str = num2str(kSel);

    %% brace on right knee
    S = SDefault;
    S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25, ...
        {'knee_angle_r'},kSel};

    S.subject.name = 'Falisse_et_al_2022';
    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
    S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',['MCainKneek' kSel_str]);
    add_pred_sim_to_batch(S,osim_path)


    %% brace on right ankle
    S = SDefault;
    S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25, ...
        {'ankle_angle_r'},kSel};

    S.subject.name = 'Falisse_et_al_2022';
    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
    S.subject.save_folder  = fullfile(pathRepoFolder,'PredSimResults',['MCainAnklek' kSel_str]);
    add_pred_sim_to_batch(S,osim_path)


    %% brace on right ankle and knee
    S = SDefault;
    S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25, ...
        {'ankle_angle_r'},kSel,{'knee_angle_r'},kSel};

    S.subject.name = 'Falisse_et_al_2022';
    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
    S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults',PredSimResults',['MCainAnkleKnee' kSel_str]);
    add_pred_sim_to_batch(S,osim_path)

end