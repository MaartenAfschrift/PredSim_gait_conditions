%% Predictive Simulations of Human Gait

% This script starts the predictive simulation of human movement. The
% required inputs are necessary to start the simulations. Optional inputs,
% if left empty, will be taken from getDefaultSettings.m.

clear; close all; clc
% path to the repository folder
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
% path to the folder that contains the repository folder
[pathRepoFolder,~,~] = fileparts(pathRepo);
saveFolderMain = 'E:\u0088756\SimResults';

%% Initialize S
pathDefaultSettings = [pathRepo '\DefaultSettings'];
addpath(pathDefaultSettings);
[S] = initializeSettings();
S.misc.main_path = pathRepo;
addpath([S.misc.main_path '\VariousFunctions'])

%% Required inputs

% either choose "quasi-random" or give the path to a .mot file you want to use as initial guess
S.subject.IG_selection = fullfile(S.misc.main_path,'OCP','IK_Guess_Default.mot');
% S.subject.IG_selection = 'quasi-random';
S.subject.IG_selection_gaitCyclePercent = 50;

% Do you want to run the simulation as a batch job (parallel computing toolbox)
S.solver.run_as_batch_job = 1;

%% Optional inputs
% see README.md in the main folder for information about these optional
% inputs.

% % S.bounds
S.bounds.a.lower            = 0.01;
% S.solver.CasADi_path        = 'C:\Users\Maarten\Documents\Software\downloads\casadi_355';
S.solver.CasADi_path        = 'C:\GBW_MyPrograms\casadi_3_5_5';
S.subject.mtp_type          = '2022paper';
S.subject.set_stiffness_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},25};
S.subject.set_damping_coefficient_selected_dofs = {{'mtp_angle_l','mtp_angle_r'},2};

% %S.Cpp2Dll: required inputs to convert .osim to .dll
% S.Cpp2Dll.PathCpp2Dll_Exe = InstallOsim2Dll_Exe('C:\Osim2Dll_exe');
S.Cpp2Dll.PathCpp2Dll_Exe = 'C:\GBW_MyPrograms\Osim2Dll_exe_vT\Cpp2Dll_Bin'; %InstallOsim2Dll_Exe('C:\GBW_MyPrograms\Osim2Dll_exe_vT');
S.Cpp2Dll.compiler = 'Visual Studio 15 2017 Win64';

% S.solver.N_threads      = 1;
% S.solver.N_meshes       = 50;
% S.solver.par_cluster_name = 'Cores4';
S.solver.N_threads      = 2;
S.solver.N_meshes       = 50;
S.solver.par_cluster_name = 'Cores10';

%% Loop over speeds and models and add to batch

SubjNames = {'Fall22_slope5','Fall22_slope-5','Falisse_et_al_2022'};

VSpeeds = [0.6667;
    0.8611;
    1.0556;
    1.2500;
    1.4444;
    1.6389;
    1.8333;
    2.0278;
    2.4167;
    2.6111;
    2.8056;
    3.0000];

for s = 1:length(SubjNames)
    % select subject
    S.subject.name = SubjNames{s};
    osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);
    % loop over walking velocities
    for i=1:length(VSpeeds)
        SpeedName = num2str(round(VSpeeds(i)*100));
        OutName = [S.subject.name '_' SpeedName];
        S.subject.v_pelvis_x_trgt   = VSpeeds(i);
        S.subject.save_folder  = fullfile(saveFolderMain,'Abe2015',OutName);
        add_pred_sim_to_batch(S,osim_path);
    end
end

