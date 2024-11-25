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
saveFolderMain = 'C:\Users\u0088756\Documents\Software\Data\ResultsPred';

%% Initialize S
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings);
[S] = initializeSettings();
S.misc.main_path = pathRepo;
addpath([S.misc.main_path '\VariousFunctions'])

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
S.subject.IG_selection_gaitCyclePercent = 50;

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

% % S.bounds
S.bounds.a.lower            = 0.01;
S.solver.CasADi_path        = 'C:\GBW_MyPrograms\casadi_3_5_5';
S.subject.mtp_type          = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};

% %S.Cpp2Dll: required inputs to convert .osim to .dll
S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\GBW_MyPrograms\Osim2Dll_exe_vT\Cpp2Dll_Bin'; %InstallOsim2Dll_Exe('C:\GBW_MyPrograms\Osim2Dll_exe_vT');
S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';


S.solver.N_threads      = 1;
S.solver.N_meshes       = 50;
S.solver.par_cluster_name = 'Cores8';

%% Loop over speeds and add to batch

S.subject.name = 'Falisse_et_al_2022';
osim_path = fullfile(pathRepo,'Subjects',S.subject.name,...
    [S.subject.name '.osim']);
PelvisHeightV = [0.89 0.88 0.87];
S.subject.v_pelvis_x_trgt = 1;

for i=1:length(PelvisHeightV)
    % adapt pelvis height
    S.bounds.coordinates = {{'pelvis_ty'},[],PelvisHeightV(i)};

    % out name
    Str_PelvHeight = num2str(round(PelvisHeightV(i)*100));
    OutName = ['PelviH_' Str_PelvHeight 'ms1'];
    S.subject.save_folder  = fullfile(saveFolderMain,'PredSimResults',OutName);
    add_pred_sim_to_batch(S,osim_path)
end

