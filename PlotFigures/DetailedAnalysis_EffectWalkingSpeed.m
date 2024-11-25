%% Detailed compared kinematics and kinetics various walking speeds
%--------------------------------------------------------------------

% files exported iwht Kinematics_Kinetics_VanDerZee on ubuntu laptop
MainDPath = 'C:\Users\mat950\OneDrive - Vrije Universiteit Amsterdam\Onderzoek\SimResults';
dPathExp = fullfile(MainDPath,'literatureSimpleInt\ExtractData\VanDerZee');

% load experimental dat
WalkSpeed = [0.7 0.9 1.1 1.4 1.6 1.8 2 ];
for i=1:length(WalkSpeed)
    IKfile = fullfile(dPathExp,['mean_' num2str(WalkSpeed(i)*10) '_IK.csv']);
    Dat(i).IK = readtable(IKfile);
    IDfile = fullfile(dPathExp,['mean_' num2str(WalkSpeed(i)*10) '_ID.csv']);
    Dat(i).ID = readtable(IDfile);
    GRFfile = fullfile(dPathExp,['mean_' num2str(WalkSpeed(i)*10) '_GRF.csv']);
    Dat(i).GRF = readtable(GRFfile);
end

PlotFalisse = true;

%% Load simulation results
if PlotFalisse
    ResPathVanDerZee = fullfile(MainDPath,'VanDerZee');
else
    ResPathVanDerZee = fullfile(MainDPath,'VanDerZee_Dhondt');
end
trialnames = [2, 5, 8, 32, 11, 14, 16, ];
for i = 1:length(trialnames)
    simfolder= fullfile(ResPathVanDerZee,['trial' num2str(trialnames(i))]);
    matfiles = dir(fullfile(simfolder,'*.mat'));
    Dat(i).Sim = load(fullfile(matfiles(1).folder,matfiles(1).name));
end

%% Plot figure
SetFigureDefaults
ColorExp = copper(length(WalkSpeed));
ColNames = {'ankle_angle_r','knee_angle_r','hip_flexion_r'};
nr = 3;
nc = 6;
figure()
for i =1:length(WalkSpeed)
    for j = 1:length(ColNames)
        % experimental kinematics
        subplot(nr,nc,j)
        dsel = Dat(i).IK.(ColNames{j});
        if strcmp(ColNames{j},'knee_angle_r')
            dsel = -dsel;
        end
        plot(dsel,'Color',ColorExp(i,:)); hold on;

        % experimental kinetics
        subplot(nr,nc,j+nc)
        dsel = Dat(i).ID.(ColNames{j});
        if strcmp(ColNames{j},'knee_angle_r')
            dsel = -dsel;
        end
        plot(dsel,'Color',ColorExp(i,:)); hold on;

        % simulated kinematics
        dsel = Dat(i).Sim.R.kinematics.Qs(:,strcmp(ColNames{j},Dat(1).Sim.R.colheaders.coordinates));
        subplot(nr,nc,j+3)
        plot(dsel,'Color',ColorExp(i,:)); hold on;

        dsel = Dat(i).Sim.R.kinetics.T_ID(:,strcmp(ColNames{j},Dat(1).Sim.R.colheaders.coordinates));
        subplot(nr,nc,j+nc+3)
        plot(dsel,'Color',ColorExp(i,:)); hold on;



    end


    subplot(nr,nc,nc*2+1)
    dsel = Dat(i).GRF.Flx;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,nc*2+2)
    dsel = Dat(i).GRF.Fly;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,nc*2+3)
    dsel = Dat(i).GRF.Flz;
    plot(dsel,'Color',ColorExp(i,:)); hold on;


    subplot(nr,nc,nc*2+4)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,1)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,nc*2+5)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,2)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,nc*2+6)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,3)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;


end


for i =1:nc*nr
    subplot(nr,nc,i)
    set(gca,'box','off');
    set(gca,'FontSize',10);
end

%% Plot separate figures for kinematics, kinetics and GRF
ColNamesHeader = {'ankle','knee','Hip'};
figure('Name','kinematics');
ColorExp = sky(length(WalkSpeed)+2);
for i =1:length(WalkSpeed)
    for j = 1:length(ColNames)
        % experimental kinematics
        subplot(2,3,j)
        dsel = Dat(i).IK.(ColNames{j});
        if strcmp(ColNames{j},'knee_angle_r')
            dsel = -dsel;
        end
        plot(dsel*180/pi,'Color',ColorExp(i+1,:)); hold on;
        title(ColNamesHeader{j})
        % simulated kinematics
        dsel = Dat(i).Sim.R.kinematics.Qs(:,strcmp(ColNames{j},Dat(1).Sim.R.colheaders.coordinates));
        subplot(2,3,j+3)
        plot(dsel,'Color',ColorExp(i+1,:)); hold on;
    end
end
for i =1:6
    subplot(2,3,i)
    set(gca,'box','off');
    set(gca,'FontSize',10);
    if i > 4
        xlabel('% gait cycle');
    end
    if i ==1
        ylabel({'Experiment','joint angle [deg]'});
    end
    if i ==4
        ylabel({'Simulation','joint angle [deg]'});
    end
    if i == 1 || i == 4
        set(gca,'YLim',[-20 30]);
    end
    if i == 2 || i ==5
        set(gca,'YLim',[-80 10]);
    end
    if i == 3 || i == 6
        set(gca,'YLim',[-30 60]);
    end
end

figure('Name','kinetics');
ColorExp = sky(length(WalkSpeed)+2);
% ColorExp = linspecer(length(WalkSpeed)+2);

for i =1:length(WalkSpeed)
    for j = 1:length(ColNames)
        % experimental kinetics
        subplot(2,3,j)
        dsel = Dat(i).ID.(ColNames{j});
        if strcmp(ColNames{j},'knee_angle_r')
            dsel = -dsel;
        end
        plot(dsel,'Color',ColorExp(i,:)); hold on;
        plot(dsel,'Color',ColorExp(i+1,:)); hold on;
        title(ColNamesHeader{j})
        % simulated kinematics
        subplot(2,3,j+3)
        dsel = Dat(i).Sim.R.kinetics.T_ID(:,strcmp(ColNames{j},Dat(1).Sim.R.colheaders.coordinates));
        plot(dsel./Dat(i).Sim.model_info.mass,'Color',ColorExp(i,:)); hold on;
    end
end
for i =1:6
    subplot(2,3,i)
    set(gca,'box','off');
    set(gca,'FontSize',10);
    if i > 4
        xlabel('% gait cycle');
    end
    if i ==1
        ylabel({'Experiment','joint moment [Nm/kg]'});
    end
    if i ==4
        ylabel({'Simulation','joint moment [Nm/kg]'});
    end
    if i == 1 || i == 4
        set(gca,'YLim',[-2 0.5]);
    end
    if i == 2 || i ==5
        set(gca,'YLim',[-1 1]);
    end
    if i == 3 || i == 6
        set(gca,'YLim',[-2 2]);
    end
end

% detailed figure ground reaction forces
figure('Name','GRF');
nr =2; nc=3;
for i =1:length(WalkSpeed)
    subplot(nr,nc,1)
    dsel = Dat(i).GRF.Flx;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,2)
    dsel = Dat(i).GRF.Fly;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,3)
    dsel = Dat(i).GRF.Flz;
    plot(dsel,'Color',ColorExp(i,:)); hold on;


    subplot(nr,nc,4)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,1)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,5)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,2)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
    subplot(nr,nc,6)
    dsel = Dat(i).Sim.R.ground_reaction.GRF_r(:,3)./Dat(i).Sim.model_info.mass;
    plot(dsel,'Color',ColorExp(i,:)); hold on;
end
for i =1:6
    subplot(2,3,i)
    set(gca,'box','off');
    set(gca,'FontSize',10);
    if i > 4
        xlabel('% gait cycle');
    end
    if i ==1
        ylabel({'Experiment','ground reactin force [N/kg]'});
    end
    if i ==4
        ylabel({'Simulation','ground reactin force [N/kg]'});
    end
    if i == 1 || i == 4
        set(gca,'YLim',[-5 6]);
    end
    if i == 2 || i ==5
        set(gca,'YLim',[0 15]);
    end
    if i == 3 || i == 6
        set(gca,'YLim',[-4 4]);
    end
end
subplot(2,3,1)
title('anterior-posterior')
subplot(2,3,2)
title('vertical')
subplot(2,3,3)
title('medio-lateral')