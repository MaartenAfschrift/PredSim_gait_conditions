% Detail berekening mechanische arbeid spiervezels in verschillende
% simulaties

datapath = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\SimResults';

% adapat the datafiles

% Default model falisse
Fal13 = load(fullfile(datapath,'Fal22Ref_13','Falisse_et_al_2022_job25.mat'));
Fal13.Mech = getMechPower(Fal13.R, Fal13.model_info);

% model lars
Dhondt13 = load(fullfile(datapath,'DhondtRef_13','DHondt_2023_2seg_job158.mat')); 
Dhondt13.Mech = getMechPower(Dhondt13.R, Dhondt13.model_info);

% model min jointwork
Fal13MinW = load(fullfile(datapath,'Fal22Ref_13_minWork_v3','Falisse_et_al_2022_job240.mat'));
Fal13MinW.Mech = getMechPower(Fal13MinW.R, Fal13MinW.model_info);

% display COT in all models
disp('COT margaria models [J/kg/m]')
disp(['  Falisse default ', num2str(Fal13.Mech.COT_Marg)]);
disp(['  Dhondt default ', num2str(Dhondt13.Mech.COT_Marg)]);
disp(['  Falisse min work ', num2str(Fal13MinW.Mech.COT_Marg)]);
disp(' ')

% display positive work muscle fibers per stride
disp('Positive muscle fibers per stride [J]')
disp(['  Falisse default ', num2str(Fal13.Mech.MusclePosWorkTotal)]);
disp(['  Dhondt default ', num2str(Dhondt13.Mech.MusclePosWorkTotal)]);
disp(['  Falisse min work ', num2str(Fal13MinW.Mech.MusclePosWorkTotal)]);
disp(' ')

% display negative work muscle fibers per stride
disp('Negative muscle fibers per stride [J]')
disp(['  Falisse default ', num2str(Fal13.Mech.MuscleNegWorkTotal)]);
disp(['  Dhondt default ', num2str(Dhondt13.Mech.MuscleNegWorkTotal)]);
disp(['  Falisse min work ', num2str(Fal13MinW.Mech.MuscleNegWorkTotal)]);
disp(' ')

% display energy lost in ground contact
disp('Energy lost in ground contact [J]')
disp(['  Falisse default ', num2str(-Fal13.Mech.TotalWork)]);
disp(['  Dhondt default ', num2str(-Dhondt13.Mech.TotalWork)]);
disp(['  Falisse min work ', num2str(-Fal13MinW.Mech.TotalWork)]);
disp(' ')

% display joint damping
disp('Joint damping [J]')
disp(['  Falisse default ', num2str(Fal13.Mech.A_damping_tot)]);
disp(['  Dhondt default ', num2str(Dhondt13.Mech.A_damping_tot)]);
disp(['  Falisse min work ', num2str(Fal13MinW.Mech.A_damping_tot)]);
disp(' ')

% efficiency 
figure();
histogram(Fal13.Mech.Bhargava_efficiency)

%% effiency uphill

% Default model falisse
McDonald = load(fullfile(datapath,'McDonald_wAct_Slope_24_ms1','Fall22_slope24_job116.mat'));
McDonald.Mech = getMechPower(McDonald.R, McDonald.model_info);

% efficiency 
figure();
histogram(McDonald.Mech.Bhargava_efficiency)
