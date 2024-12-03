% Detailed figure Mconald
%-------------------------


% path information

% path information other software and datapath
MainPath = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions';
ResPathMcDonald = 'C:\Users\mat950\Documents\Data\SimResults_Afschrift2025\PredSimResults\';
ExpData = fullfile(MainPath, 'ExperimentalData');
dExpMcDonald = fullfile(ExpData,'McDonald\EnergyCost_FatigueAvoidance_MasterData.mat');
addpath(genpath(fullfile(MainPath,'NeuromechanicsToolkit')));
addpath(genpath(fullfile(MainPath,'PredSim')));

%% read data

% table information
nSimOutput = 14+7*2 + 2; iSimOut = 1:nSimOutput; % see SimResults2Table for now
nExpOutput = 5; iExpOut = nSimOutput+1:nSimOutput+nExpOutput;
Table = nan(1000,nSimOutput+nExpOutput);
SimInfo =cell(1000,3); % model label, location added mass
ctTable = 1;
Header_Sim = {'COT','StrideFreq','StepWidth','speed','PNetMetab','BodyMass','Slope',...
    'PNetMetab_BM','AddedMass','COT_marg','PNetMetab_Marg_BM','COT_b','PNetMetab_b','PNetMetab_b_BM',...
    'act1','act2','act3','act4','act5','act6','act7','act1sq','act2sq','act3sq','act4sq','act5sq',...
    'act6sq','act7sq','Avgact-Mass','Avgact_sq-Mass'};
Header_Exp = {'Exp_PNetMetab_BM','Exp_COT','Exp_StrideFreq','Exp_StepWidth','Exp_CAct'};
Headers = [Header_Sim, Header_Exp];

% simulation data
Slopes= [0 6 12 18 24];
% PelvisHeightV = [0.88 0.86 0.84];
PelvisHeightV = 0.88;
SlopesV = [Slopes, 0 0 0];
SubjNames = {'Fal22Ref_10'}; % reference simulation walking 1m/s
ctSubjName = 2;
CondHeader{1} = 'None';
for i=2:length(Slopes)
    slope_str = num2str(round(Slopes(i)));
    SubjNames{ctSubjName} = ['McDonald_Slope_' slope_str '_ms1'];
    CondHeader{ctSubjName} = 'None';
    ctSubjName = ctSubjName+1;
end
for i=1:length(PelvisHeightV)
    Str_PelvHeight = num2str(round(PelvisHeightV(i)*100));
    SubjNames{ctSubjName} =  ['PelviH_' Str_PelvHeight 'ms1'];
    CondHeader{ctSubjName} = 'Crouch';
    ctSubjName = ctSubjName+1;
end

% experimental data
Exp = load(dExpMcDonald);
ObjectiveValue = nan(length(SubjNames),8);
for i= 1:length(SubjNames)
    % simulation data
    ResFolder  = fullfile(ResPathMcDonald,SubjNames{i});
    [outputSim] = SimResults2Table_vMcDonald(ResFolder,'slope',SlopesV(i),'RefModelMass',...
        62, 'AddedMass',0);
    Table(ctTable,iSimOut) = outputSim;
    SimInfo{ctTable,1} = 'McDonald';
    SimInfo{ctTable,2} = CondHeader{i};

    % experimental data
    if i <6
        Pmet = nanmean(Exp.C_metP(:,i+1));
        CAct = nanmean(Exp.C_asq(:,i+1));
    elseif i > 5
        Pmet = nanmean(Exp.C_metP(:,1));
        CAct = nanmean(Exp.C_asq(:,1));
    else
        Pmet = NaN;
    end
    dExp = nan(1,nExpOutput);
    dExp(1) = Pmet; % net power
    dExp(2) = Pmet/1; % cost of transport
    dExp(3) = NaN; % stride frequency
    dExp(4) = NaN; % stride width in unknown
    dExp(5) = CAct;
    Table(ctTable,iExpOut) = dExp;
    ctTable = ctTable + 1;

    % store value objective function
    % ObjectiveValue
    [Res] = LoadSimFile(ResFolder);
    % get objective function values
    ObjectiveValue(i,:) = Res.R.objective.absoluteValues;

end
ObjectiveValue_total = sum(ObjectiveValue,2);
Table(ctTable:end,:) = [];

% compute activation cost similar as in
iRef = find(Table(:,strcmp(Headers,'Slope')) == 0);
iRef = iRef(1);

% cost compute cost activation as in McDonald
ActRef = Table(iRef,15:21);
Act = Table(:,15:21)./repmat(ActRef,ctTable-1,1);
CostAct = sum(Act.^2,2)./7;
Table(:,nSimOutput+nExpOutput+1) = CostAct;
ActRef = Table(iRef,22:28);
Act = Table(:,22:28)./repmat(ActRef,ctTable-1,1);
Table(:,nSimOutput+nExpOutput+2) = sum(Act,2)./7;
Headers = [Headers, {'CostActivationSq','CostActivationSq2'}];

%% plot figure McDonald

nc = 2;
nr = 1;
figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic');
mk = 7;

% plot Pmetab exp and sim
subplot(nc,nr,1)
plot([0 5], [0 5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_PNetMetab_BM';
yHeader = 'PNetMetab_b_BM';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);

% plot activation cost
subplot(nc,nr,2)
plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_CAct';
yHeader = 'CostActivationSq';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta ActCost - Exp','Interpreter','tex')
ylabel('\Delta ActCost - Sim','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);

%% plot figure McDonald: No Delta

nc = 2;
nr = 1;
figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic');
mk = 7;

% plot Pmetab exp and sim
subplot(nc,nr,1)
plot([0 10], [0 10],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_PNetMetab_BM';
yHeader = 'PNetMetab_b_BM';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(Table(iSel,strcmp(Headers,xHeader)),...
    Table(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(Table(iCrouch,strcmp(Headers,xHeader)),...
    Table(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel(' Pnet Exp. (W/kg)','Interpreter','tex')
ylabel(' Pnet Sim. (W/kg)','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);

% plot activation cost
subplot(nc,nr,2)
plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_CAct';
yHeader = 'CostActivationSq';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(Table(iSel,strcmp(Headers,xHeader)),...
    Table(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('ActCost - Exp','Interpreter','tex')
ylabel('ActCost - Sim','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);


%% Plot cost function in simulation and compare it with measured costs

figure('Color',[1 1 1],'Name','Figure Objective');
subplot(2,2,1)
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
bar(1,ObjectiveValue_total(1)); hold on;
bar(2,ObjectiveValue_total(2));
bar(3,ObjectiveValue_total(3));
bar(4,ObjectiveValue_total(4));
bar(5,ObjectiveValue_total(5));
bar(6,ObjectiveValue_total(6));
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'Level','6%','12%','18%','24%','crouch'});
title('total');

subplot(2,2,2)
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
bar(1,ObjectiveValue(1,1)); hold on;
bar(2,ObjectiveValue(2,1));
bar(3,ObjectiveValue(3,1));
bar(4,ObjectiveValue(4,1));
bar(5,ObjectiveValue(5,1));
bar(6,ObjectiveValue(6,1));
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'Level','6%','12%','18%','24%','crouch'});
title('metab');

subplot(2,2,3)
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
bar(1,ObjectiveValue(1,2)); hold on;
bar(2,ObjectiveValue(2,2));
bar(3,ObjectiveValue(3,2));
bar(4,ObjectiveValue(4,2));
bar(5,ObjectiveValue(5,2));
bar(6,ObjectiveValue(6,2));
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'Level','6%','12%','18%','24%','crouch'});
title('act');

subplot(2,2,4)
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
bar(1,ObjectiveValue(1,4)); hold on;
bar(2,ObjectiveValue(2,4));
bar(3,ObjectiveValue(3,4));
bar(4,ObjectiveValue(4,4));
bar(5,ObjectiveValue(5,4));
bar(6,ObjectiveValue(6,4));
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'Level','6%','12%','18%','24%','crouch'});
title('qdd');

for i=1:4
    subplot(2,2,i)
    set(gca,'Box','off');
end

%% Same plot but different layout

SetFigureDefaults();
figure('Color',[1 1 1],'Name','Objective Func v2');
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
slope = [0 6 12 18 24];
nr = 2; nc = 3;

mk = 7;

subplot(nr,nc,1);
plot(slope, ObjectiveValue_total(1:5),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue_total(6),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue_total(6) ObjectiveValue_total(6)],'Color',Cs_crouch,'LineStyle','--')

subplot(nr,nc,2);
plot(slope, ObjectiveValue(1:5,1),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,1),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,1) ObjectiveValue(6,1)],'Color',Cs_crouch,'LineStyle','--')

subplot(nr,nc,3);
plot(slope, ObjectiveValue(1:5,2),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,2),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,2) ObjectiveValue(6,2)],'Color',Cs_crouch,'LineStyle','--')

subplot(nr,nc,4);
plot(slope, ObjectiveValue(1:5,4),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,4),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,4) ObjectiveValue(6,4)],'Color',Cs_crouch,'LineStyle','--')

subplot(nr,nc,5);
yHeader = 'PNetMetab_b_BM';
plot(slope, DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0,DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);
line([0 24],[DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader))],'Color',Cs_crouch,'LineStyle','--')

subplot(nr,nc,6);
yHeader = 'CostActivationSq';
plot(slope, DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0,DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);
line([0 24],[DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader))],'Color',Cs_crouch,'LineStyle','--')

for i=1:nr*nc
    subplot(nr,nc,i)
    set(gca,'box','off');
    set(gca,'FontSize',11);
end

%% Plot figure for paper

nc = 3;
nr = 2;
figure('Color',[1 1 1],'Name','Figure Paper');
set(gcf,'Position',[468         690        1092         548]);
mk = 7;

% plot Pmetab exp and sim
subplot(nr,nc,2)
plot([0 5], [0 5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_PNetMetab_BM';
yHeader = 'PNetMetab_b_BM';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

% legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')


% plot activation cost
subplot(nr,nc,3)
plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_CAct';
yHeader = 'CostActivationSq';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

% legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta ActCost - Exp','Interpreter','tex')
ylabel('\Delta ActCost - Sim','Interpreter','tex')


subplot(nr,nc,4);
plot(slope, ObjectiveValue_total(1:5),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue_total(6),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue_total(6) ObjectiveValue_total(6)],'Color',Cs_crouch,'LineStyle','--')
ylabel('total objective');
xlabel('slope [%]');

subplot(nr,nc,5);
plot(slope, ObjectiveValue(1:5,1),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,1),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,1) ObjectiveValue(6,1)],'Color',Cs_crouch,'LineStyle','--')
ylabel('objective metabolic power');
xlabel('slope [%]');

subplot(nr,nc,6);
plot(slope, ObjectiveValue(1:5,2),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,2),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,2) ObjectiveValue(6,2)],'Color',Cs_crouch,'LineStyle','--')
ylabel('objective muscle activation');
xlabel('slope [%]');

for i =1:nr*nc
    subplot(nr,nc,i)
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);
    set(gca, 'FontName', 'Arial')
end

%% plot figure McDonald -- squared

nc = 2;
nr = 1;
figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic');
mk = 7;

% plot Pmetab exp and sim
subplot(nc,nr,1)
plot([0 5], [0 5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_PNetMetab_BM';
yHeader = 'PNetMetab_b_BM';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);

% plot activation cost
subplot(nc,nr,2)
% plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_CAct';
yHeader = 'CostActivationSq2';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

% legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta ActCost - Exp','Interpreter','tex')
ylabel('\Delta ActCost - Sim','Interpreter','tex')
set(gca,'box','off')
set(gca,'FontSize',14);
set(gca,'LineWidth',1.5);

%% Plot figure with objective function values as a function of the slope

% wait on results for multiple slopes ? (currently running)

% temporary figure
%
for i= 1:length(SubjNames)-1
    % simulation data
    ResFolder  = fullfile(ResPathMcDonald,SubjNames{i});
end

%% Plot mass based cost


nc = 3;
nr = 2;
figure('Color',[1 1 1],'Name','Figure Paper');
set(gcf,'Position',[468         690        1092         548]);
mk = 7;

% plot Pmetab exp and sim
subplot(nr,nc,1)
plot([0 5], [0 5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_PNetMetab_BM';
yHeader = 'PNetMetab_b_BM';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

% legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')


% plot activation cost
subplot(nr,nc,2)
plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
DeltaTable = Table - Table(iRef,:);
xHeader = 'Exp_CAct';
yHeader = 'CostActivationSq';
Cs = [0.1 0.1 0.9];
Cs_crouch = [0.9 0.1 0.1];
iSel = 1:5;
iCrouch = 6:length(Table(:,1));
l1 = plot(DeltaTable(iSel,strcmp(Headers,xHeader)),...
    DeltaTable(iSel,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

l2 = plot(DeltaTable(iCrouch,strcmp(Headers,xHeader)),...
    DeltaTable(iCrouch,strcmp(Headers,yHeader)),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk);

% legend([l1(1), l2(1)],{'Slope','Crouch'})
xlabel('\Delta ActCost - Exp','Interpreter','tex')
ylabel('\Delta ActCost - Sim','Interpreter','tex')


subplot(nr,nc,3);
plot(slope, ObjectiveValue_total(1:5),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue_total(6),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue_total(6) ObjectiveValue_total(6)],'Color',Cs_crouch,'LineStyle','--')
ylabel('total objective');
xlabel('slope [%]');

subplot(nr,nc,4);
plot(slope, ObjectiveValue(1:5,1),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,1),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,1) ObjectiveValue(6,1)],'Color',Cs_crouch,'LineStyle','--')
ylabel('objective metabolic power');
xlabel('slope [%]');

subplot(nr,nc,5);
plot(slope, ObjectiveValue(1:5,2),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
plot(0, ObjectiveValue(6,2),...
    'o','Color',Cs_crouch,'MarkerFaceColor',Cs_crouch,'MarkerSize',mk); hold on;
line([0 24],[ObjectiveValue(6,2) ObjectiveValue(6,2)],'Color',Cs_crouch,'LineStyle','--')
ylabel('objective muscle activation');
xlabel('slope [%]');

subplot(nr,nc,6);
iSlope = 1:5;
Cs = [0.1 0.1 0.9];
yHeader = 'Avgact-Mass';
% yHeader = 'Avgact_sq-Mass';
l2 = plot(slope,Table(iSlope,strcmp(Headers,yHeader)),...
    'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
line([0 24],[Table(iCrouch,strcmp(Headers,yHeader)), Table(iCrouch,strcmp(Headers,yHeader))],'Color',Cs_crouch,'LineStyle','--')
ylabel('muscle activation x mass');
xlabel('slope [%]');

for i =1:nr*nc
    subplot(nr,nc,i)
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);
    set(gca, 'FontName', 'Arial')
end
