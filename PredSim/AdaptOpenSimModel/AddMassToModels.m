%% Build osim model with added mass
%-----------------------------------


currentfilename = mfilename('fullpath');
currentfolder = fileparts(currentfilename);
pathRepo = currentfolder(1:end-18);

modelpath = fullfile(pathRepo,'Subjects','Falisse_et_al_2022');
modelname = 'Falisse_et_al_2022.osim';

% load the model
import org.opensim.modeling.*;


%% export model based on Gomeñuka2013

% 25% of body mass is added (assumed located in backpack)
% load model
m = Model(fullfile(modelpath,modelname));

% get the total mass of the subject
bodies = m.getBodySet();
mtot = 0;
for i=1:bodies.getSize()
    mtot = mtot + bodies.get(i-1).getMass();
end
% 25% added mass
addedmass = 0.25*mtot;
% add mass to the COM of the torso
torso = bodies.get('torso');
torso.setMass(torso.getMass() + addedmass);
% print the model
modelnameOut = 'Falisse_2022_25pcMassTorso.osim';
m.print(fullfile(modelpath,modelnameOut));

%% export model based on Gomeñuka2013: mass behind COM torso

% 25% of body mass is added (assumed located in backpack)
% load model
m = Model(fullfile(modelpath,modelname));

% get the total mass of the subject
bodies = m.getBodySet();
mtot = 0;
for i=1:bodies.getSize()
    mtot = mtot + bodies.get(i-1).getMass();
end
% 25% added mass
addedmass = 0.25*mtot;
% add mass to the COM of the torso
torso = bodies.get('torso');
% adapt the COM location
mtorso = torso.getMass();
COM = torso.getMassCenter();
COM_mat = [COM.get(0) COM.get(1) COM.get(2)];
COM_mass = COM_mat; 
COM_mass(1) = COM_mass(1)-0.15;
COM_new = (COM_mat*mtorso + COM_mass*addedmass)./(mtorso + addedmass);
for i=0:2
    COM.set(i, COM_new(i+1));
end

% adapt inertia (we should adapt moment around z-axis as we change the x location)
I = torso.getInertia;
I_diag = I.getMoments;
dBackPack = sqrt((COM_mat(1)-COM_new(1)).^2 +(COM_mat(2)-COM_new(2)).^2);
dMass = sqrt((COM_mass(1)-COM_new(1)).^2 +(COM_mass(2)-COM_new(2)).^2);
I_new_z = I_diag.get(2) + mtorso * dBackPack.^2 + 0 + addedmass * dMass.^2;
I_diag.set(2, I_new_z);
torso.setInertia(I);

% adapt mass
torso.setMass(torso.getMass() + addedmass);

% init model
m.initSystem();

% print the model
modelnameOut = 'Fal22_25pcMassBackPack.osim';
m.print(fullfile(modelpath,modelnameOut));







