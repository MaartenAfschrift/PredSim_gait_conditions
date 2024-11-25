
%% In this script we want to make variations of the Falisse 2022 model

% the goal is to first build the cpp file
% adapt the cpp file manually if needed
% convert the cpp to dll

clear all; close all; clc;

%% Path information

pathRepo = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\PredSim';
CreateDllExePath = fullfile(pathRepo,'Osim2DLL');
Compiler = 'Visual Studio 15 2017 Win64';

% original model
DefaultModelname = 'Falisse_et_al_2022';
DefaultOsimPath = fullfile(pathRepo,'Subjects',[DefaultModelname '.osim']);

% output path
OutPathModels = fullfile(pathRepo,'Subjects');

%% Create models for walking on slopes

SlopesV = -0.24:0.02:0.24;

for i = 1:length(SlopesV)
    
    % init settings structure
    [S] = initializeSettings();
    S.misc.main_path = pathRepo;

    % Original model settings
    S.subject.name = DefaultModelname;
    S.osim_path_Or1 = DefaultOsimPath;

    % read the model and adapt
    import org.opensim.modeling.*;
    modSel = Model(S.osim_path_Or1);

    % get the gravtiy vector
    g = modSel.getGravity();
    gv = [g.get(0) g.get(1) g.get(2) 0];
    slope = SlopesV(i); % slope of 24%
    fi = asin(slope);
    R = rotz(fi);
    grav_tilt = gv*R;
    gravSlope = Vec3(grav_tilt(1), grav_tilt(2), grav_tilt(3));
    modSel.setGravity(gravSlope);
    slope_str = num2str(round(slope*100));
    SubjName = ['Fall22_slope' slope_str];
    Outpath_model_sel = fullfile(OutPathModels,SubjName);
    OutPathModel = fullfile(Outpath_model_sel,[SubjName '.osim']);
    modSel.print(OutPathModel);

    % conver the model
    S.subject.name = SubjName;

    S.Cpp2Dll.PathCpp2Dll_Exe = CreateDllExePath;
    S.Cpp2Dll.compiler = Compiler;
    % export information
    S.Cpp2Dll.export3DSegmentOrigins = {'calcn_r', 'calcn_l', 'femur_r', 'femur_l',...
        'hand_r','hand_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l'};
    S.Cpp2Dll.jointsOrder = [];
    S.Cpp2Dll.coordinatesOrder = [];
    S.Cpp2Dll.exportGRFs = true;
    S.Cpp2Dll.exportSeparateGRFs = true;
    S.Cpp2Dll.exportGRMs = true;
    S.Cpp2Dll.exportContactPowers = true;
    S.Cpp2Dll.verbose_mode = 1;

    osim2cpp(S.Cpp2Dll,Outpath_model_sel,OutPathModel)

    % convert cpp
    outputFilename = ['Fall22_slope' slope_str];
    verbose = 1;
    cpp2dll(Outpath_model_sel,outputFilename,[],Compiler,S.Cpp2Dll.PathCpp2Dll_Exe,verbose);
end

