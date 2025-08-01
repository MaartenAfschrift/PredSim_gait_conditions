
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

%% Create models for walking on slopes

Slopes = [0 0.07 0.15 0 0.07 0.15];
MassTorso = [zeros(1,3), zeros(1,3)+0.25];
ModelNamesOut = {'Fal_level','Fal_level7','Fal_level15',...
    'Fal_level_mtors','Fal_level7_mtors','Fal_level15_mtors'};

for i = 4:length(Slopes)
    
    % init settings structure
    [S] = initializeSettings();
    S.misc.main_path = pathRepo;

    % Original model settings
    S.subject.name = DefaultModelname;
    S.osim_path_Or1 = DefaultOsimPath;

    % read the model and adapt
    import org.opensim.modeling.*;
    modSel = Model(S.osim_path_Or1);

    % adapt the gravity vector of the opensim model
    g = modSel.getGravity();
    gv = [g.get(0) g.get(1) g.get(2) 0];
    slope = Slopes(i); % slope of 24%
    fi = asin(slope);
    R = rotz(fi);
    grav_tilt = gv*R;
    gravSlope = Vec3(grav_tilt(1), grav_tilt(2), grav_tilt(3));
    modSel.setGravity(gravSlope);
    
    % get the total mass of the subject
    bodies = modSel.getBodySet();
    mtot = 0;
    for ij=1:bodies.getSize()
        mtot = mtot + bodies.get(ij-1).getMass();
    end
    % x% added mass
    madd = MassTorso(i)*mtot;    
    if madd ~= 0 
        BodySel = modSel.getBodySet().get('torso');
        mBodySel = BodySel.getMass();
        BodySel.setMass(mBodySel + madd);
    end
        
    % print the model
    SubjName = ModelNamesOut{i};
    PathOsimMod = fullfile(pathRepo,'Subjects',SubjName);
    if ~isfolder(PathOsimMod)
        mkdir(PathOsimMod);
    end
    OutPathModel = fullfile(PathOsimMod,[SubjName,'.osim']);  
    modSel.initSystem();
    modSel.print(OutPathModel);
    
    % convert the model
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

    osim2cpp(S.Cpp2Dll,PathOsimMod,OutPathModel)

    % convert cpp
    outputFilename = S.subject.name;
    verbose = 1;
    cpp2dll(PathOsimMod,outputFilename,[],Compiler,S.Cpp2Dll.PathCpp2Dll_Exe,verbose);
end

