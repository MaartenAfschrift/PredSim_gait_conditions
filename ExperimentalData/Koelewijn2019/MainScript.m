%% Path information

clear all; close all; clc;
addpath(genpath('C:\Users\mat950\Documents\Software\DataAnalysis\NeuromechanicsToolkit'));
Datapath = 'E:\Data\Koelewijn2019';

load(fullfile(Datapath,'KinData'),'KinData');
load(fullfile(Datapath,'MetabolicCostData'),'MetCost');

% Anne did an excellent job. We can extract all the data from this
% structure

%Trial 1	Level, 1.3 m/s
%Trial 2	Level, 0.8 m/s
%Trial 3	Uphill, 1.3 m/s
%Trial 4	Uphill, 0.8 m/s
%Trial 5	Downhill, 1.3 m/s
%Trial 6	Downhill, 0.8 m/s

DatKoel = nan(6,5,12);
HeaderKoel = {'speed','slope','stride time','stepfreq','Pmetab'};
SpeedVect = [1.3, 0.8, 1.3, 0.8, 1.3, 0.8];
grade = [0, 0, 8, 8, -8, -8];
ct = 1;
for i=1:6
    for s = 1:12
        dSel = KinData.(['Subject' num2str(s)]).(['Trial' num2str(i)]);
        DatKoel(ct,1,s) = SpeedVect(i);
        DatKoel(ct,2,s) = grade(i);
        DatKoel(ct,3,s) = dSel.dur;
        DatKoel(ct,4,s) = SpeedVect(i)./dSel.dur;
        if i == 1
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.walk.level;
        elseif  i == 2
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.slow.level;
        elseif  i == 3
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.walk.positive;
        elseif  i == 4
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.slow.positive;
        elseif  i ==5
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.walk.negative;
        elseif  i == 6
            Pmetab = MetCost.(['Subject' num2str(s)]).Experiment.slow.negative;
        end
        DatKoel(ct,5,s) = Pmetab;        
    end
    ct = ct+1;
end

save('DataTableKoelewijn.mat','DatKoel','HeaderKoel');






