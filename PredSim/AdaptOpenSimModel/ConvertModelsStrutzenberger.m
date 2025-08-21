
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
S.osim_path_Or1 = fullfile(pathRepo,'Subjects',S.subject.name_Or1,[S.subject.name_Or1 '.osim']);

%% New folder

ListFolders = {'Falisse2022_12pDecline','Falisse2022_6pDecline','Falisse2022_Level',...
    'Falisse2022_6pIncline','Falisse2022_12pIncline',};

for s= 1:length(ListFolders)
    S.subject.name = ListFolders{s};
    S.osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

    %% Copy the model to a new folder

    % copy the opensim model to a new folder
    if ~isfolder(fullfile(pathRepo,'Subjects',S.subject.name))
        mkdir(fullfile(pathRepo,'Subjects',S.subject.name));
    end
    copyfile(S.osim_path_Or1,S.osim_path);


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

%% adapt gravivty vector for slope


grav = [0 -9.80664999999999942304 0 0];
VecSlopes = [0.15 0.08 0.07 -0.07 -0.08 -0.15];
for s= 1:length(VecSlopes)
    slope = VecSlopes(s); % slope in %
    fi = asin(slope);
    R = rotz(fi);
    grav_tilt = grav*R;
    disp(['tilted gravity vector ', num2str(grav_tilt,12) ]);
end


