

%% In this script we want to make variations of the Falisse 2022 model

% the goal is to first build the cpp file
% adapt the cpp file manually if needed
% convert the cpp to dll

clear all; close all; clc;

%% Path information


currentfilename = mfilename('fullpath');
currentfolder = fileparts(currentfilename);
pathRepo = currentfolder(1:end-18);


CreateDllExePath = fullfile(pathRepo,'Osim2DLL','Cpp2Dll_Bin');
Compiler = 'Visual Studio 17 2022';

% original model
DefaultModelname = 'Falisse_et_al_2022';
DefaultOsimPath = fullfile(pathRepo,'Subjects',DefaultModelname,[DefaultModelname '.osim']);


% output path
OutPathModels = fullfile(pathRepo,'Subjects');


%% Create models for walking on slopes

SlopesV = [-0.08 0.08];
ListFolders = {'Falisse2022_8pDecline','Falisse2022_8pIncline'};


for i = 1:length(ListFolders)

    % init settings structure
    [S] = initializeSettings();
    S.misc.main_path = pathRepo;

    % default settings
    S.subject.name = DefaultModelname;
    S.osim_path_Or1 = DefaultOsimPath;
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

    % read the model and adapt
    import org.opensim.modeling.*;
    modSel = Model(S.osim_path_Or1);

    % get and adapt the gravtiy vector
    g = modSel.getGravity();
    gv = [g.get(0) g.get(1) g.get(2) 0];
    slope = SlopesV(i); % slope of 24%
    fi = asin(slope);
    R = rotz(fi);
    grav_tilt = gv*R;
    gravSlope = Vec3(grav_tilt(1), grav_tilt(2), grav_tilt(3));
    modSel.setGravity(gravSlope);
    slope_str = num2str(round(slope*100));
    SubjName = ListFolders{i};
    Outpath_model_sel = fullfile(OutPathModels,SubjName);
    OutPathModel = fullfile(Outpath_model_sel,[SubjName '.osim']);
    if ~isfolder(Outpath_model_sel)
        mkdir(Outpath_model_sel)
    end
    modSel.print(OutPathModel);

    % conver the model to a cpp file
    S.subject.name = SubjName;
    osim2cpp(S.Cpp2Dll,Outpath_model_sel,OutPathModel)

    % convert cpp
    outputFilename = SubjName;
    verbose = 1;
    cpp2dll(Outpath_model_sel,outputFilename,[],Compiler,S.Cpp2Dll.PathCpp2Dll_Exe,verbose);
end

