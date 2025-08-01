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
% name of the subject
S.subject.name = 'Falisse2022_7pIncline';

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
S.subject.IG_selection_gaitCyclePercent = 50;

% give the path to the osim model of your subject
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
S.bounds.a.lower            = 0.01;
S.solver.CasADi_path        = get_casadi_path();
S.subject.mtp_type          = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};

% %S.Cpp2Dll: required inputs to convert .osim to .dll
S.Cpp2Dll.PathCpp2Dll_Exe = fullfile(pathRepo,'Osim2DLL');
S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';

S.solver.N_threads      = 1;
S.solver.N_meshes       = 50;
S.solver.par_cluster_name = 'Cores4';

%% Loop over speeds and models and add to batch

% different models for level walking and walking on slopes
SubjNames = {'Falisse_et_al_2022','Falisse2022_7pDecline','Falisse2022_15pDecline',...
    'Falisse2022_7pIncline','Falisse2022_15pIncline',...
    'Fal22_25pMTorso','Fal22_25pMTorso_15pIncline','Fal22_25pMTorso_7pIncline',...
    'Fal22_25pMTorso_7pDecline','Fal22_25pMTorso_15pDecline'};

%  SubjNames = {'Falisse_et_al_2022','Falisse2022_7pIncline','Falisse2022_15pIncline',...
%      'Fal22_25pMTorso','Fal22_25pMTorso_15pIncline','Fal22_25pMTorso_7pIncline'};
% SubjNames = {'Fal22_25pMBackP','Fal22_25pMBackP_7pIncline',...
%     'Fal22_25pMBackP_15pIncline'};


VSpeeds = (2:0.5:7)./3.6;
SpeedNames = {'20kmh','25kmh','30kmh','35kmh','40kmh','45kmh','50kmh','55kmh','60kmh','65kmh','70kmh'};
for s = 1:length(SubjNames)
    % select subject
    S.subject.name = SubjNames{s};
    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
    % loop over walking velocities 
    for i=1:length(VSpeeds)
        OutName = [S.subject.name '_' SpeedNames{i}];
        S.subject.v_pelvis_x_trgt   = VSpeeds(i);
        S.subject.save_folder  = fullfile(pathRepoFolder,'SimResults','PredGom2013_BackP',OutName);
        add_pred_sim_to_batch(S,osim_path)
    end
end

