
%% In this script we want to make variations of the Falisse 2022 model

% the goal is to first build the cpp file
% adapt the cpp file manually if needed
% convert the cpp to dll

clear all; close all; clc;

currentfilename = mfilename('fullpath');
currentfolder = fileparts(currentfilename);
pathRepo = currentfolder(1:end-18);


PathCpp2Dll_Exe = fullfile(pathRepo,'Osim2DLL','Cpp2Dll_Bin');
compiler = 'Visual Studio 17 2022';


%% add mass to models

import org.opensim.modeling.*;
% persons walked at 1.25 m/s

Addedmass = 4:1:20;
SegmentAdded = {'torso',[]};
COMlocation = {[0 -0.1 0], []}; % 10cm below COM torso, this is about 20cm above pelvis
mBackPack = 3.8;
% default model
GenModel = fullfile(pathRepo,'Subjects','ModelRepo','Falisse_et_al_2022.osim');

% output model
for i=1:length(Addedmass)
    ModelOut{i} = ['Huang_' num2str(Addedmass(i)) 'kg_BackP'];
end

for i=1:length(Addedmass)    
    % open model
    mSel = Model(GenModel);
    for j =1:2
        segSel = SegmentAdded{1,j};
        if ~isempty(segSel)
            % added mass, COM lcation
            madd = Addedmass(i) + mBackPack;
            madd_COM = COMlocation{1,j};
            madd_Inertia = [0, 0, 0];

            % adapt the segment mass
            BodySel = mSel.getBodySet().get(segSel);
            mBodySel = BodySel.getMass();
            BodySel.setMass(mBodySel + madd);
            % get COM location and inertia
            I = BodySel.getInertia;
            I_diag = I.getMoments;
            I_diag_mat = [I_diag.get(0) I_diag.get(1) I_diag.get(2)];
            COM = BodySel.getMassCenter();
            COM_mat = [COM.get(0) COM.get(1) COM.get(2)];
            % compute new COM location and inertia
            [COM_new,I_new] = AdaptSegmentCOMAndInetia(COM_mat,I_diag_mat,mBodySel,...
                COM_mat+madd_COM,madd_Inertia, madd);
            % update the model
            for ii=1:3
                I_diag.set(ii-1,I_new(ii));
                COM.set(ii-1,COM_new(ii))
            end
            BodySel.setInertia(I);
            BodySel.setMassCenter(COM);
        end        
    end
    % print the model
    mSel.initSystem();
    OutDir = fullfile(pathRepo,'Subjects',ModelOut{i});
    OutModel = fullfile(OutDir,[ModelOut{i} '.osim']);
    if ~isfolder(OutDir)
        mkdir(OutDir);
    end
    mSel.print(OutModel);
    clear mSel;
    % copy the muscle analysis file
    copyfile(fullfile(pathRepo,'Subjects','ModelRepo','f_lMT_vMT_dM_poly_3_9'),...
        fullfile(OutDir,'f_lMT_vMT_dM_poly_3_9'));
end

%% Convert all models to cpp files
for s= 1:length(Addedmass) 
    

    % get osim path and model name
    S.subject.name = ModelOut{s};
    S.osim_path = fullfile(pathRepo,'Subjects',S.subject.name,[S.subject.name '.osim']);

    % test if dll file exists
    dllfile = fullfile(pathRepo,'Subjects',S.subject.name,['F_' S.subject.name '.dll']);
    if ~exist(dllfile,'file')

    % Create a cpp file from a model
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

    % convert cpp
    CppDir = fullfile(pathRepo,'Subjects',S.subject.name);
    outputFilename = S.subject.name;
    verbose = 0;
    cpp2dll(CppDir,outputFilename,[],compiler,PathCpp2Dll_Exe,verbose);
    else
        disp([dllfile ' already converted'])
    end
end

