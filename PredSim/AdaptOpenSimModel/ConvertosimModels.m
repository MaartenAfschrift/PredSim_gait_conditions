
%% In this script we want to make variations of the Falisse 2022 model

% the goal is to first build the cpp file
% adapt the cpp file manually if needed
% convert the cpp to dll

clear all; close all; clc;

currentfilename = mfilename('fullpath');
currentfolder = fileparts(currentfilename);
pathRepo = currentfolder(1:end-18);

PathCpp2Dll_Exe = fullfile(pathRepo,'Osim2DLL');
compiler = 'Visual Studio 17 2022';
BuildCPP = true;
BuildDll = true;


[S] = initializeSettings();
S.misc.main_path = pathRepo;
%% Original model settings
S.subject.name_Or1 = 'Falisse_et_al_2022';
S.osim_path_Or1 = fullfile(pathRepo,'Subjects','ModelRepo',[S.subject.name_Or1 '.osim']);

S.subject.name_Or2 = 'Falisse_2022_25pcMassTorso';
S.osim_path_Or2 = fullfile(pathRepo,'Subjects','ModelRepo',[S.subject.name_Or2 '.osim']);

S.subject.name_Or3 = 'Fal22_25pcMassBackPack';
S.osim_path_Or3 = fullfile(pathRepo,'Subjects','ModelRepo',[S.subject.name_Or3 '.osim']);

%% New folder

ListFolders = {'Falisse_et_al_2022','Falisse2022_7pDecline','Falisse2022_15pDecline',...
    'Falisse2022_7pIncline','Falisse2022_15pIncline',...
    'Fal22_25pMTorso','Fal22_25pMTorso_15pIncline','Fal22_25pMTorso_7pIncline',...
    'Fal22_25pMTorso_7pDecline','Fal22_25pMTorso_15pDecline',...
    'Fal22_25pMBackP','Fal22_25pMBackP_15pIncline','Fal22_25pMBackP_7pIncline'};

for s=1:length(ListFolders)
    S.subject.name = ListFolders{s};
    S.osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

    %% Copy the model to a new folder

    % copy the opensim model to a new folder
    if ~isfolder(fullfile(pathRepo,'Subjects',S.subject.name))
        mkdir(fullfile(pathRepo,'Subjects',S.subject.name));
    end
    if  s < 6
        copyfile(S.osim_path_Or1,S.osim_path);
    elseif s < 11
        copyfile(S.osim_path_Or2,S.osim_path);
    else
        copyfile(S.osim_path_Or3,S.osim_path);
    end

    %% Create a cpp file from a model
    % name of the subject
    if BuildCPP
        % give the path to the osim model of your subject
        CppDir = fullfile(pathRepo,'Subjects',S.subject.name);
        osim_path = S.osim_path;

        % export information
        Cpp2Dll.export3DSegmentOrigins = {'calcn_r', 'calcn_l', 'femur_r', 'femur_l',...
            'hand_r','hand_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l'};
        Cpp2Dll.jointsOrder = [];
        Cpp2Dll.coordinatesOrder = [];
        Cpp2Dll.exportGRFs = true;
        Cpp2Dll.exportSeparateGRFs = true;
        Cpp2Dll.exportGRMs = true;
        Cpp2Dll.exportContactPowers = true;
        Cpp2Dll.verbose_mode = 0;
        Cpp2Dll.PathCpp2Dll_Exe = PathCpp2Dll_Exe;

        % osim to cpp
        osim2cpp(Cpp2Dll,CppDir,osim_path)
    end


    %% Convert cpp file to .dll
    if BuildDll
        % convert cpp
        CppDir = fullfile(pathRepo,'Subjects',S.subject.name);
        outputFilename = S.subject.name;
        verbose = 1;
        cpp2dll(CppDir,outputFilename,[],compiler,PathCpp2Dll_Exe,verbose);
    end
end



