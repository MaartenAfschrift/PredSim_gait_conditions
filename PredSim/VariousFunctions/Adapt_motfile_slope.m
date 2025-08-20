% adapt solution for slope walking.
%   we have to rotate the pelvis tilt and change the x and y coordinates
%   depending on the slope

% load solution
res_path = 'C:\Users\mat950\OneDrive - Vrije Universiteit Amsterdam\Onderzoek\SimResults\PredSimResults\Fal_Slope_24_ms1';
mat_files = dir(fullfile(res_path,'*.mat'));
res = load(fullfile(mat_files(1).folder,mat_files(1).name));

% slope is 8 pct in this example
slope = 24;
ang = atan2(slope,100);
e_slope = [cos(ang), sin(ang), 0];

% read .mot file
mot_filename = fullfile(mat_files(1).folder,[mat_files(1).name(1:end-4) '.mot']);
mot = ReadMotFile(mot_filename);

% adapt mot file
data = mot.data;
R = rotz(ang)';
R = R(1:3,1:3);

i_pelvis_tilt = strcmp(mot.names,'pelvis_tilt');
data(:,i_pelvis_tilt) = mot.data(:,i_pelvis_tilt)+ang*180/pi;

i_pelvis_tx = find(strcmp(mot.names,'pelvis_tx'));
data_tpelvis =data(:,i_pelvis_tx:i_pelvis_tx+2);

data(:,i_pelvis_tx) = data_tpelvis*R(:,1);
data(:,i_pelvis_tx+1) = data_tpelvis*R(:,2);
data(:,i_pelvis_tx+2) = data_tpelvis*R(:,3);
data_pelvis_out = data(:,i_pelvis_tx:i_pelvis_tx+2); 

% export mot file
mot_filename_out = fullfile(mat_files(1).folder,[mat_files(1).name(1:end-4) '_slope.mot']);
generateMotFile(data,mot.names,mot_filename_out)







