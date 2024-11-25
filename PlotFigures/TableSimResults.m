%% Create table with all simulation results
%-------------------------------------------
% main script used to plot figures in the publucation

clear all; clc;

Settings.CreateDataTable = false;
Settings.PlotFiguresPaper = true;
Settings.PlotFiguresExplore = true;
Settings.DeltaValues = true;
Settings.vAbe2015_lim = 4; % limit to walking conditions, or keep everything ?

% path information other software
MainDPath = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\SimResults';
ExpData = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\ExperimentalData';
addpath(genpath('C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\NeuromechanicsToolkit'));
addpath(genpath('C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\PredSim'));

% path all simulation and experiments

% relative paths simulation results
RespathBrowning = fullfile(MainDPath,'Browning');
ResPathSchertzer  = fullfile(MainDPath,'Schertzer');
ResPathRef = fullfile(MainDPath,'PredSimResults');
ResPathHuang = fullfile(MainDPath,'Huang2014');
ResPathSchertzer_wObj1 = fullfile(MainDPath,'Schertzer__wMetab60pc'); % 60% weight on metabolic cost
ResPathSchertzer_wObj2 = fullfile(MainDPath,'Schertzer__wqdd30pc'); % 60% weight on metabolic cost
ResPathSchertzer_QR = fullfile(MainDPath,'Schertzer_QR'); % quasi random initial guess
ResPathSchertzer_QRm72 = fullfile(MainDPath,'Schertzer_QR_m72'); % quasi random initial guessm model mass 72 kg
ResPathSchertzer_QRm82 = fullfile(MainDPath,'Schertzer_QR_m82'); % quasi random initial guessm model mass 82 kg
ResPathGom = fullfile(MainDPath,'PredGom2013');
ResPathGom_wObj1 = fullfile(MainDPath,'PredGom2013_wMetab60pc');
ResPathKoel = fullfile(MainDPath,'Koelewijn');
ResPathMcDonald = fullfile(MainDPath,'PredSimResults');
RespathUmberger = fullfile(MainDPath,'PredSimResults');
ResPathAbe = fullfile(MainDPath,'Abe2015');
ResPathVanDerZee = fullfile(MainDPath,'VanDerZee');

% relative path to experimental data (extracted from data using digitizer)
dPathBrowning= fullfile(ExpData,'Browning2008');
dPathSchertzer = fullfile(ExpData,'schertzer2014');
ExpDataFileGom = fullfile(ExpData,'Gomenuka2013\Gomenuka2013_b.csv');
dExpDatKoel= fullfile(ExpData,'Koelewijn2019\COT_Koelewijn2019.csv');
dExpDatKoel2= fullfile(ExpData,'Koelewijn2019\DataTableKoelewijn.mat');
dExpMcDonald = fullfile(ExpData,'McDonald\EnergyCost_FatigueAvoidance_MasterData.mat');
dExpUmberger = fullfile(ExpData,'Umberger2007','Umberger2007.csv');
dExpJordan = fullfile(ExpData,'Jordan2006','Jordan2006.csv');
dExpVanDerZee = fullfile(ExpData,'VanDerZee','stridefreq.csv');

% subject properties to nondim data
Exp_SubjProp_header = {'body_mass','height'};
Studyname = {'Gomenuka','Browning','Schertzer','Huang','Koelewijn',...
    'VanDerZee2022','McDonald','Abe2015','Strutzenberger'};
prop_leg_length = 0.55;
Mass = [71.6, 74.16, 74.88, 71.1, 70, 73.5, 69.6, 58.9, 73.1];
Height = [1.78, 1.82, 1.78, 0.99/prop_leg_length, 1.73, 1.76, 1.70, 1.70, 1.77];
LegLength = Height.*prop_leg_length;


Mass_Sim = 62;
Height_Sim = 1.70;
LegLength_Sim = Height_Sim.*prop_leg_length;
Exp_SubjProp = table(Studyname',Mass',Height',LegLength','VariableNames',...
    {'Study','mass','height','LegLength'});


%% create data table
if Settings.CreateDataTable

    % pre allocate table
    nSimOutput = 27; iSimOut = 1:nSimOutput; % see SimResults2Table for now
    nExpOutput = 4; iExpOut = nSimOutput+1:nSimOutput+nExpOutput;
    Table = nan(1000,nSimOutput+nExpOutput);
    SimInfo =cell(1000,3); % model label, location added mass
    ctTable = 1;
    Header_Sim = {'COT','StrideFreq','StepWidth','speed','PNetMetab','BodyMass','Slope',...
        'PNetMetab_BM','AddedMass','COT_marg','PNetMetab_Marg_BM','COT_b','PNetMetab_b',...
        'PNetMetab_b_BM','COT_b_hr','PNetMetab_b_hr','PNetMetab_b_BM_hr','NetMechWork',...
        'PosMuscleWork','NegMuscleWork','mean_act','mean_actsq','NetMechPower','NetPosMusclePower',...
        'NetNegMusclePower','NetMetabPower_nobasal','margaria_negworkbased'};
    Header_Exp = {'Exp_PNetMetab_BM','Exp_COT','Exp_StrideFreq','Exp_StepWidth'};
    Headers = [Header_Sim Header_Exp];
    Headers_Info = {'ExpModelInfo','LocationAddedMass'};


    %% Create Table - Schertzer 2013

    % Unfortunatly no data related to spatio-temporal outcomes

    % To Do: currently we don't use the correct speed in the reference
    % conditions. Run these simulations again for 4 5 6 kmh
    % loaded conditions
    SubjNames = {'Ref';'Ankle_05';'Ankle_1';'Ankle_2';...
        'Torso_2';'Torso_7';'Torso_10'; 'Torso_16'; 'Torso_22';...
        'Knee_05';'Knee_1';'Knee_2'};
    AddedMassVect = [0 1 2 4 2 7 10 16 22 1 2 4]; % mass on each leg
    LocationAdded = {'None';'Ankle';'Ankle';'Ankle';...clc
        'Torso';'Torso';'Torso'; 'Torso'; 'Torso';...
        'Knee';'Knee';'Knee'};

    % speeds
    VSpeeds = [4 5 6]./3.6;
    SpeedNames = {'40kmh','50kmh','60kmh'};
    SpeedNamesRef = {'11','14','17'};

    % multiple simulations models
    ModelSettings = {ResPathSchertzer,ResPathSchertzer_wObj1,ResPathSchertzer_wObj2,...
        ResPathSchertzer_QR,ResPathSchertzer_QRm72,ResPathSchertzer_QRm82};
    ModelSettings_label = {'Schertzer','Schertzer_wMetab','Schertzer_wqdd','Schertzer_qr',...
        'Schertzer_qrm72','Schertzer_qrm82'};

    % specific names for walking speed conditions without added mass
    PostFixRef = {'','_metab','_qdd30pc','_QR','_m72','_m82'};


    % load experimental data
    % Schertzer reports net metabolic rate divided by (bodymass + addedmass)
    mRef = 62;
    Schert.Exp.Ankle.Tab = importdata(fullfile(dPathSchertzer,'AddedMassAnkle.csv'));
    Schert.Exp.Ankle.speed = Schert.Exp.Ankle.Tab(:,1)./3.6;
    Schert.Exp.Ankle.Pnet = Schert.Exp.Ankle.Tab(:,2);
    Schert.Exp.Ankle.COT = Schert.Exp.Ankle.Tab(:,2)./Schert.Exp.Ankle.speed;
    Schert.Exp.Ankle.mass = [zeros(3,1); zeros(3,1)+0.5; zeros(3,1)+1; zeros(3,1)+2]*2; % times 2 because on both legs

    Schert.Exp.Knee.Tab = importdata(fullfile(dPathSchertzer,'AddedMassKnee.csv'));
    Schert.Exp.Knee.speed = Schert.Exp.Knee.Tab(:,1)./3.6;
    Schert.Exp.Knee.Pnet = Schert.Exp.Knee.Tab(:,2);
    Schert.Exp.Knee.COT = Schert.Exp.Knee.Tab(:,2)./Schert.Exp.Knee.speed;
    Schert.Exp.Knee.mass = [zeros(3,1); zeros(3,1)+0.5; zeros(3,1)+1; zeros(3,1)+2]*2; % times 2 because on both legs

    Schert.Exp.Torso.Tab = importdata(fullfile(dPathSchertzer,'AddeMassTrunk.csv'));
    Schert.Exp.Torso.speed = Schert.Exp.Torso.Tab(:,1)./3.6;
    Schert.Exp.Torso.Pnet = Schert.Exp.Torso.Tab(:,2);
    Schert.Exp.Torso.COT = Schert.Exp.Torso.Tab(:,2)./Schert.Exp.Torso.speed;
    Schert.Exp.Torso.mass = [zeros(3,1)+0; zeros(3,1)+2; zeros(3,1)+7; zeros(3,1)+10;...
        zeros(3,1)+16; zeros(3,1)+22];
    % The COT curve is independent of the added mass for the torso, so we copy the same COT
    % multiple times.
    Schert.Exp.Torso.speed = repmat(Schert.Exp.Torso.speed,6,1);
    Schert.Exp.Torso.COT = repmat(Schert.Exp.Torso.COT,6,1);
    Schert.Exp.Torso.Pnet = repmat(Schert.Exp.Torso.Pnet,6,1);

    % Pnet normalized by body mass (and not by body mass + added mass)
    Schert.Exp.Torso.PnetBM = Schert.Exp.Torso.Pnet.*(mRef + Schert.Exp.Torso.mass)./mRef;
    Schert.Exp.Knee.PnetBM = Schert.Exp.Knee.Pnet.*(mRef + Schert.Exp.Knee.mass)./mRef;
    Schert.Exp.Ankle.PnetBM = Schert.Exp.Ankle.Pnet.*(mRef + Schert.Exp.Ankle.mass)./mRef;


    % loop over all simulations
    for s = 1:length(SubjNames)
        % loop over specific settings
        mAdded = AddedMassVect(s);
        for ms = 1:length(ModelSettings)
            % loop over walking velocities
            for i=1:length(VSpeeds)
                if strcmp(SubjNames{s},'Ref')
                    % simulation data
                    OutName = ['Fal22Ref_' SpeedNamesRef{i} PostFixRef{ms}];
                    if ms == 1 ||  ms == 2 || ms == 3 || ms == 4
                        mRefModel = 62;
                    elseif ms == 5
                        mRefModel = 72;
                    elseif ms == 6
                        mRefModel = 82;
                    end
                    OutPath = fullfile(ResPathRef,OutName);
                else
                    % simulation data
                    if ms == 1 ||  ms == 2 || ms == 3 || ms == 4
                        OutName = [SubjNames{s} '_' SpeedNames{i}];
                        mRefModel = 62;
                    elseif ms == 5
                        OutName = [SubjNames{s} 'm72_' SpeedNames{i}];
                        mRefModel = 72;
                    elseif ms == 6
                        OutName = [SubjNames{s} 'm82_' SpeedNames{i}];
                        mRefModel = 82;
                    end
                    OutPath = fullfile(ModelSettings{ms},OutName);
                end

                [outputSim] = SimResults2Table(OutPath,'slope',0,'RefModelMass',mRefModel, ...
                    'AddedMass',mAdded);
                Table(ctTable,iSimOut) = outputSim;
                SimInfo{ctTable,1} = ModelSettings_label{ms};
                SimInfo{ctTable,2} = LocationAdded{s};

                % find the experimental data that corresponds with this
                % simulation
                dExp = nan(1,nExpOutput);
                if strcmp(LocationAdded{s},'Ankle') || strcmp(LocationAdded{s},'None')
                    % find exp data with correct mass
                    iSel = find(Schert.Exp.Ankle.mass == mAdded & ...
                        (abs(Schert.Exp.Ankle.speed - VSpeeds(i)) < 0.05));
                    if ~isempty(iSel)
                        dExp(1) = Schert.Exp.Ankle.PnetBM(iSel); % net power
                        dExp(2) = Schert.Exp.Ankle.COT(iSel); % cost of transport
                        dExp(3) = NaN; % stride length in unknown
                        dExp(4) = NaN; % stride width in unknown
                    else
                        disp(['Schertzer: cannot find condition in exp data for mass: ' num2str(mAdded), ...
                            ' at Ankle, speed: ' num2str(VSpeeds(i))]);
                    end
                elseif strcmp(LocationAdded{s},'Knee')
                    iSel = find(Schert.Exp.Knee.mass == mAdded & ...
                        (abs(Schert.Exp.Knee.speed - VSpeeds(i)) < 0.05));
                    if ~isempty(iSel)
                        dExp(1) = Schert.Exp.Knee.PnetBM(iSel); % net power
                        dExp(2) = Schert.Exp.Knee.COT(iSel); % cost of transport
                        dExp(3) = NaN; % stride length in unknown
                        dExp(4) = NaN; % stride width in unknown
                    else
                        disp(['Schertzer: cannot find condition in exp data for mass: ' num2str(mAdded), ...
                            ' at Knee, speed: ' num2str(VSpeeds(i))]);
                    end
                elseif strcmp(LocationAdded{s},'Torso')
                    iSel = find(Schert.Exp.Torso.mass == mAdded & ...
                        (abs(Schert.Exp.Torso.speed - VSpeeds(i)) < 0.05));
                    if ~isempty(iSel)
                        dExp(1) = Schert.Exp.Torso.PnetBM(iSel); % net power
                        dExp(2) = Schert.Exp.Torso.COT(iSel); % cost of transport
                        dExp(3) = NaN; % stride frequency
                        dExp(4) = NaN; % stride width in unknown
                    else
                        disp(['Schertzer: cannot find condition in exp data for mass: ' num2str(mAdded), ...
                            ' at Torso, speed: ' num2str(VSpeeds(i))]);
                    end
                end
                Table(ctTable,iExpOut) = dExp;
                ctTable = ctTable + 1;
            end
        end
    end

    %% Create Table - Gomenuka

    % Unfortunatly no data related to spatio-temporal outcomes

    % simulation info
    SubjNames = {'Falisse_et_al_2022','Falisse2022_7pIncline','Falisse2022_15pIncline',...
        'Fal22_25pMTorso','Fal22_25pMTorso_7pIncline','Fal22_25pMTorso_15pIncline'};
    LocationAdded = {'None','None','None','Torso','Torso','Torso'};
    grade = [0 7 15 0 7 15];
    mRefModel = 62;
    addedmass = [0 0 0 0.25 0.25 0.25].*mRefModel;
    SpeedNames = {'20kmh','30kmh','40kmh','50kmh','60kmh'};
    WalkSpeed = 2:1:6;
    ModelSettings = {ResPathGom, ResPathGom_wObj1};
    ModelSettings_label = {'Gomenuka','Gomenuka_wMetab'};

    % load experimental data
    Dat = importdata(ExpDataFileGom);
    Gom2013.Fal22_25pMTorso_15pIncline = Dat(2:6,:).*1.25; % times 1.25 because we want W/kg body mass
    Gom2013.Falisse2022_15pIncline = Dat(7:11,:);
    Gom2013.Fal22_25pMTorso_7pIncline = Dat(12:16,:).*1.25;
    Gom2013.Falisse2022_7pIncline = Dat(17:21,:);
    Gom2013.Fal22_25pMTorso = Dat(22:26,:).*1.25;
    Gom2013.Falisse_et_al_2022 = Dat(27:31,:);


    for s = 1:length(SubjNames)
        if ~isempty(SubjNames{s})
            for i=1:length(SpeedNames)
                for ms = 1:length(ModelSettings_label)
                    % load the results file (there should be one mat file with the results in this folder)
                    OutName = [SubjNames{s} '_' SpeedNames{i}];
                    OutPath = fullfile(ModelSettings{ms},OutName);
                    [outputSim] = SimResults2Table(OutPath,'slope',grade(s),'RefModelMass',mRefModel, ...
                        'AddedMass',addedmass(s));
                    Table(ctTable,iSimOut) = outputSim;
                    SimInfo{ctTable,1} = ModelSettings_label{ms};
                    SimInfo{ctTable,2} = LocationAdded{s};

                    % add experimental data
                    dExp = nan(1,nExpOutput);
                    dsel = Gom2013.(SubjNames{s});
                    dExp(1) = dsel(i,2).*(WalkSpeed(i)./3.6); % net power
                    dExp(2) = dsel(i,2); % cost of transport
                    dExp(3) = NaN; % stride length in unknown
                    dExp(4) = NaN; % stride width in unknown
                    Table(ctTable,iExpOut) = dExp;
                    ctTable = ctTable + 1;
                end
            end
        end
    end

    %% Create Table - a

    % Raw data availlable here: https://zenodo.org/record/1973799.
    % To Do: Processs the data with the AddBiomechanics tool

    % simulation data
    SubjNames = {'Falisse2022_8pDecline','LevelWalk','Falisse2022_8pIncline'};
    Grade = [-8 0 8];
    VSpeeds = [0.8 1.3];
    SpeedNames = {'0c8ms','1c3ms'};
    mRefModel = 62;

    % experimental data
    Koel.Exp.Table = importdata(dExpDatKoel);
    Koel.Exp.grade = [-8 -8 0 0 8 8];
    Koel.Exp.speed = [0.8 1.3 0.8 1.3 0.8 1.3];
    %     Koel.Exp.Tab(:,1) = Koel.Exp.grade;
    %     Koel.Exp.Tab(:,2) = Koel.Exp.speed;
    %     Koel.Exp.Tab(:,3) = Koel.Exp.Table.data;
    %     Koel.Exp.Tab(:,4) = Koel.Exp.Tab(:,3).*Koel.Exp.Tab(:,2);
    %     Koel.Exp.Tab_header = {'grade','speed','COT','MetabPower'};
    % also based on .mat file (raw data)
    Koel.Exp2 = load(dExpDatKoel2);
    DKoel = Koel.Exp2.DatKoel;
    DMean = nanmean(DKoel,3);
    Koel.Exp.Tab = DMean;
    Koel.Exp.Header = Koel.Exp2.HeaderKoel;
    %     disp('ToDo: controlleer exp data Koelewijn')

    % add results to table
    for s = 1:length(SubjNames)
        for i=1:length(VSpeeds)
            % simulation data
            if s == 2
                ResFolderSel  = fullfile(ResPathRef,['Fal22Ref_' num2str(round(VSpeeds(i)*10))]);
            else
                OutName = [SubjNames{s} '_' SpeedNames{i}];
                ResFolderSel  = fullfile(ResPathKoel,OutName);
            end
            [outputSim] = SimResults2Table(ResFolderSel,'slope',Grade(s),...
                'RefModelMass',mRefModel,'AddedMass',0);
            Table(ctTable,iSimOut) = outputSim;
            SimInfo{ctTable,1} = 'Koelewijn';
            SimInfo{ctTable,2} = 'None';
            % experimental data
            iSel = Koel.Exp.Tab(:,strcmp(Koel.Exp.Header,'slope')) == Grade(s) &...
                round(Koel.Exp.Tab(:,strcmp(Koel.Exp.Header,'speed')),2) == VSpeeds(i);
            dExp = nan(1,nExpOutput);
            Pmetab = Koel.Exp.Tab(iSel,strcmp(Koel.Exp.Header,'Pmetab'));
            strideFreq = Koel.Exp.Tab(iSel,strcmp(Koel.Exp.Header,'stepfreq'));
            dExp(1) = Pmetab.*VSpeeds(i); % cost of transport
            dExp(2) = Pmetab; % net power
            dExp(3) = strideFreq; % stride length in unknown
            dExp(4) = NaN; % stride width in unknown
            Table(ctTable,iExpOut) = dExp;
            ctTable = ctTable + 1;
        end
    end


    %% Create Table - Huang

    Addedmass = [0 4:1:20];
    mRef = 62;
    walkspeed = 1.25;
    BodyMass = 62;
    LegLength = 0.9;
    g = 9.81;
    slope_metab = 0.249;
    offset_metab = 0.0859;
    slope_steplength = -0.0029;
    offset_steplength = 0.666;
    slope_stepwidth =0.0308;
    offset_stepwidth = 0.143;

    for i=1:length(Addedmass)
        %   -------  simulation data -------
        %------------------------------------
        if Addedmass(i) == 0
            OutPath  = fullfile(ResPathRef,'Fal22Ref_12');
        else
            OutPath  = fullfile(ResPathHuang,['Huang_' num2str(Addedmass(i)) 'kg']);
        end
        [outputSim] = SimResults2Table(OutPath,'slope',0,'RefModelMass',...
            mRefModel, 'AddedMass',Addedmass(i));
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'Huang';
        SimInfo{ctTable,2} = 'Torso';

        %   -------  experimental data -------
        %-------------------------------------
        % using the regression equations from table 1
        madd = Addedmass(i);
        % metabolic power
        ScalePower =  BodyMass * g ^1.5 * LegLength^0.5; % all power is scaled by this factor
        NormPower = offset_metab+ slope_metab * (madd/BodyMass); % mass is scaled by body mass
        DeltaMetab = NormPower * ScalePower;
        DeltaMetab_BM = DeltaMetab /mRef;
        dExp = nan(1,nExpOutput);
        dExp(1) = DeltaMetab_BM; % net power
        dExp(2) = DeltaMetab_BM./walkspeed; % cost of transport
        % step length
        NormSL = offset_steplength + slope_steplength * (madd/BodyMass);
        SL = NormSL*LegLength;
        stridefrequency = walkspeed/SL;
        dExp(3) = stridefrequency; % stride length in unknown
        % step width
        NormSW = offset_stepwidth + slope_stepwidth * (madd/BodyMass);
        SW = NormSW*LegLength;
        dExp(4) = SW; % stride width in unknown
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;
    end

    % add reference walking to dataset

    %% Create Table - Browning

    % read the simulation results
    SubjNames = {'Ref';'Brown_femur2'; 'Brown_femur4'; 'Brown_femur8';...
        'Brown_foot2';'Brown_foot4';'pelvis_4';'pelvis_8';'pelvis_12';'pelvis_16';...
        'tibia_2';'tibia_4'};
    addedmass = [0 4 8 16 ,...
        4 8 4 8 12 16,...
        4 8];
    LocationAdded = {'None','Femur','Femur','Femur', ...
        'Foot','Foot','Pelvis','Pelvis','Pelvis','Pelvis',...
        'Tibia','Tibia'};

    % read the experimental data
    Browning.Exp.Tab = importdata(fullfile(dPathBrowning,'Browning_COT.csv'));
    % Browning.Exp.Fr = [0.91,0; ... % no added mass
    %     0.9 16; % added mass to pelvis
    %     0.92, 16;... % added mass to femur
    %     0.87, 8;... % added mass to tibia
    %     0.85, 4;... % added mass to foot
    %     0.81, 8]; % added mass to foot
    Browning.Exp.Fr = [0.91, NaN, NaN, 0.92, ...
        0.85, 0.81, NaN, NaN, NaN, 0.9, ...
        NaN, 0.87];

    mRefModel = 62;
    IndexExpTable = [1 NaN, 6, 7, ...
        2, 3, 8 9 NaN 10,...
        4 5]; % index SubjName in exp data table
    walkspeed = 1.25;

    % extract simulation results
    for s = 1:length(SubjNames)

        if strcmp(SubjNames{s},'Ref')
            ResFolder = fullfile(ResPathRef,'Fal22Ref_125');
        else
            % add simulation results
            ResFolder = fullfile(RespathBrowning,SubjNames{s});
        end
        [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
            mRefModel, 'AddedMass',addedmass(s));
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'Browning';
        SimInfo{ctTable,2} = LocationAdded{s};

        % add experimental data
        indexSel = IndexExpTable(s);
        dExp = nan(1,nExpOutput);
        if ~isnan(indexSel)
            COT = Browning.Exp.Tab(indexSel,2);
            m = Browning.Exp.Tab(indexSel,1);
            if abs(m - addedmass(s)) < 0.2
                dExp(1) = COT*walkspeed; % net power
                dExp(2) = COT; % cost of transport
                dExp(3) = Browning.Exp.Fr(s); % stride frequency
                dExp(4) = NaN; % stride width in unknown
            else
                disp(['Problem in indexing experimental data in Browning ' SubjNames{s}]);
            end
        end
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;
    end

    %% Create Table - Umberger2007

    % preferred step frequency is 0.89 in the simulation model at 1.3 m/s

    % get experimental data
    Dat = importdata(dExpUmberger);
    Freq.ExpUmb.freqDev = Dat(:,1);
    Freq.ExpUmb.Power = Dat(:,2); % net metabolic power, norm to body mass
    Freq.ExpUmb.freq = (1+Freq.ExpUmb.freqDev/100).*0.89;

    % get simulation results
    % loop over alls simulations
    StrideFreq = 0.7:0.05:1.3;
    for i=1:length(StrideFreq)
        % filename
        freqName = num2str(round(StrideFreq(i)*100));
        ResFolder  = fullfile(RespathUmberger,['freq_' freqName]);
        [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
            62, 'AddedMass',0);
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'Umberger';
        SimInfo{ctTable,2} = 'None';
        % experimental data
        DevStrideFreq = abs(StrideFreq(i)-Freq.ExpUmb.freq);
        dExp = nan(1,nExpOutput);
        if any(DevStrideFreq<0.025)
            iSel = find(DevStrideFreq<0.025);
            dExp(1) = Freq.ExpUmb.Power(iSel); % net power
            dExp(2) = Freq.ExpUmb.Power(iSel)/1.3; % cost of transport
            dExp(3) = Freq.ExpUmb.freq(iSel); % stride frequency
            dExp(4) = NaN; % stride width in unknown
        end
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;
    end

    %% McDonald Data

    % simulation data
    Slopes= [0 6 12 18 24];
    PelvisHeightV = [0.88 0.86 0.84];
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
    for i= 1:length(SubjNames)
        % simulation data
        ResFolder  = fullfile(ResPathMcDonald,SubjNames{i});
        [outputSim] = SimResults2Table(ResFolder,'slope',SlopesV(i),'RefModelMass',...
            62, 'AddedMass',0);
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'McDonald';
        SimInfo{ctTable,2} = CondHeader{i};

        % experimental data
        if i <6
            Pmet = nanmean(Exp.C_metP(:,i+1));
        elseif i == 6
            Pmet = nanmean(Exp.C_metP(:,1));
        else
            Pmet = NaN;
        end
        dExp = nan(1,nExpOutput);
        dExp(1) = Pmet; % net power
        dExp(2) = Pmet/1; % cost of transport
        dExp(3) = NaN; % stride frequency
        dExp(4) = NaN; % stride width in unknown
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;

    end

    %% Abe 2015

    % experimental data
    [Abe_Dat,Abe_header] = GetDataAbe2015();

    % walking at various speeds and inclines
    SubjNames = {'Fall22_slope5','Fall22_slope-5','Falisse_et_al_2022'};
    Slopes = [5, -5, 0];
    VSpeeds = [0.6667;
        0.8611;
        1.0556;
        1.2500;
        1.4444;
        1.6389;
        1.8333;
        2.0278;
        2.4167;
        2.6111;
        2.8056;
        3.0000];
    for s = 1:length(SubjNames)
        for i=1:length(VSpeeds)
            SpeedName = num2str(round(VSpeeds(i)*100));
            OutName = [SubjNames{s} '_' SpeedName];
            ResFolder  = fullfile(ResPathAbe,OutName);
            [outputSim] = SimResults2Table(ResFolder,'slope',Slopes(s),'RefModelMass',...
                62, 'AddedMass',0);
            Table(ctTable,iSimOut) = outputSim;
            SimInfo{ctTable,1} = 'Abe2015';
            SimInfo{ctTable,2} = 'None';

            % experimental data
            dExp = nan(1,nExpOutput);
            % find matching speed
            ExpSpeeds = Abe_Dat(:,strcmp(Abe_header,'speed'));
            slope = Slopes(s);
            DeltaSpeeds = abs(ExpSpeeds- VSpeeds(i));
            if any(DeltaSpeeds<0.05)
                iSel = find(DeltaSpeeds<0.05 & Abe_Dat(:,strcmp(Abe_header,'grade')) == slope);
                dExp(1) = Abe_Dat(iSel,strcmp(Abe_header,'MetabPower'));
                dExp(2) = Abe_Dat(iSel,strcmp(Abe_header,'Metab_COT'));
                dExp(3) = NaN; % stride frequency
                dExp(4) = NaN; % stride width in unknown
            else
                disp('test');
            end
            Table(ctTable,iExpOut) = dExp;
            ctTable = ctTable + 1;
        end
    end

    %% VanDerZee
    DExp = importdata(dExpVanDerZee);
    walkspeed = [0.7 0.9 1.1 1.6 1.8 2 1.4];
    iTrials = [2 5 8 11 14 16 32];
    for i =1:length(walkspeed)
        % exp data
        stridefreqAv = nanmean(DExp.data(DExp.data(:,2) == walkspeed(i),3));
        sf_noNorm = stridefreqAv*(1.72*0.5);
        dExp(1) = NaN;
        dExp(2) = NaN;
        dExp(3) = sf_noNorm; % stride frequency
        dExp(4) = NaN; % stride width in unknown
        Table(ctTable,iExpOut) = dExp;
        % sim data
        ResFolder = fullfile(ResPathVanDerZee,['trial' num2str(iTrials(i))]);
        [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
            62, 'AddedMass',0);
        Table(ctTable,iSimOut) = outputSim;
        % sim info
        SimInfo{ctTable,1} = 'VanDerZee2022';
        SimInfo{ctTable,2} = 'None';
        ctTable = ctTable + 1;
    end


    %% Jordan 2006

    % experimental data
    DExp = importdata(dExpJordan);
    PreferedSpeed = 1.3;
    ExpSpeeds = DExp(:,1).*PreferedSpeed;
    ExpStrideFreq = ExpSpeeds./DExp(:,2);

    % walking at various speeds - COT
    VSpeeds = 0.5:0.1:2;
    SpeedNames = round(VSpeeds*10);
    for i=1:length(VSpeeds)
        OutName = ['Fal22Ref_' num2str(SpeedNames(i))];
        ResFolder  = fullfile(ResPathRef,OutName);
        [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
            62, 'AddedMass',0);
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'Jordan';
        SimInfo{ctTable,2} = 'None';

        % experimental data
        dExp = nan(1,nExpOutput);
        % find matching speed
        slope = 0;
        DeltaSpeeds = abs(ExpSpeeds- VSpeeds(i));
        if any(DeltaSpeeds<0.05)
            [~,iSel] = min(DeltaSpeeds);
            dExp(1) = NaN;
            dExp(2) = NaN;
            dExp(3) = ExpStrideFreq(iSel); % stride frequency
            dExp(4) = NaN; % stride width in unknown
        end
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;
    end

    %% Effect ignoring muscle dynamics

    VSpeeds = 0.5:0.1:2;
    SpeedNames = round(VSpeeds*10);
    ModelTypeNames = {'RigidTendon','RigidTendon_IgnoreFLV','Default'};
    PrefixNames = {'Fal22_RigTendon_','Fal22_RigTendon_NoFLV_','Fal22Ref_'};
    for j = 1:length(ModelTypeNames)
        for i=1:length(VSpeeds)
            % rigid tendon
            OutName = [PrefixNames{j} num2str(SpeedNames(i))];
            ResFolder  = fullfile(ResPathRef,OutName);
            [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
                62, 'AddedMass',0);
            Table(ctTable,iSimOut) = outputSim;
            SimInfo{ctTable,1} = 'WalkSpeedRef';
            SimInfo{ctTable,2} = ModelTypeNames{j};
            ctTable = ctTable + 1;
        end
    end

    %% Effect weights objective function: different walking speeds

    VSpeeds = 0.5:0.1:2;
    SpeedNames = round(VSpeeds*10);
    ModelTypeNames = {'lower_w_qdd','lower_w_metab','Default','Default_QR'};
    PostfixNames = {'_qdd30pc','_metab','','_QR'};
    for j = 1:length(ModelTypeNames)
        for i=1:length(VSpeeds)
            % rigid tendon
            OutName = ['Fal22Ref_' num2str(SpeedNames(i)) PostfixNames{j}];
            ResFolder  = fullfile(ResPathRef,OutName);
            [outputSim] = SimResults2Table(ResFolder,'slope',0,'RefModelMass',...
                62, 'AddedMass',0);
            Table(ctTable,iSimOut) = outputSim;
            SimInfo{ctTable,1} = 'WalkSpeedRef_Obj';
            SimInfo{ctTable,2} = ModelTypeNames{j};
            ctTable = ctTable + 1;
        end
    end

    %% Strutzenberger

    Slopes_exp = [-12 -6 0 6 12];
    Freq_Exp = [106 104 100 97 97]./60./2;
    Slopes_Sim = Slopes_exp;
    for i = 1:length(Slopes_Sim)
        OutName = ['Strutzenberger_', num2str(Slopes_Sim(i)) ,'_ms1c1'];
        ResFolder  = fullfile(ResPathRef,OutName);
        [outputSim] = SimResults2Table(ResFolder,'slope',Slopes_Sim(i),'RefModelMass',...
            62, 'AddedMass',0);
        Table(ctTable,iSimOut) = outputSim;
        SimInfo{ctTable,1} = 'Strutzenberger';
        SimInfo{ctTable,2} = 'none';

        %{'Exp_PNetMetab_BM','Exp_COT','Exp_StrideFreq','Exp_StepWidth'};
        % experimental data
        dExp = nan(1,nExpOutput);
        dExp(1) = NaN;
        dExp(2) = NaN;
        dExp(3) = Freq_Exp(i); % stride frequency
        dExp(4) = NaN; % stride width in unknown
        Table(ctTable,iExpOut) = dExp;
        ctTable = ctTable + 1;
    end


    %% Remove unwanted lines
    Table(ctTable:end,:) = [];
    SimInfo(ctTable:end,:) = [];

    %% save results

    save('ExpSimDataTable.mat','Table','SimInfo','Headers','Headers_Info');
else
    load('ExpSimDataTable.mat','Table','SimInfo','Headers','Headers_Info');
end

if Settings.PlotFiguresExplore
    %% Test Plot metabolic power in absolute terms

    % test plot Schertzer
    figure();
    Cs = linspecer(1);
    mk = 6;
    iSel = strcmp(SimInfo(:,1),'Schertzer')';
    plot([-1 8], [-1 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    plot(Table(iSel,strcmp(Headers,'PNetMetab_BM')),...
        Table(iSel,strcmp(Headers,'Exp_PNetMetab_BM')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
    ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    set(gca,'box','off');
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    title('Schertzer (loaded walking at various speeds');

    % test plot McDonald
    figure();
    Cs = linspecer(1);
    mk = 6;
    iSel = strcmp(SimInfo(:,1),'McDonald')';
    plot([-1 8], [-1 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    plot(Table(iSel,strcmp(Headers,'PNetMetab_BM')),...
        Table(iSel,strcmp(Headers,'Exp_PNetMetab_BM')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
    ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    set(gca,'box','off');
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    title('McDonald slope and crouch');


    % test plotAbe 2015
    figure();
    Cs = linspecer(1);
    mk = 6;
    iSel = strcmp(SimInfo(:,1),'Abe2015')' ;
    plot([-1 8], [-1 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    plot(Table(iSel,strcmp(Headers,'PNetMetab_BM')),...
        Table(iSel,strcmp(Headers,'Exp_PNetMetab_BM')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
    ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    set(gca,'box','off');
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    title('Abe2015 walking speed variations');

    % test stride frequency figure
    figure();
    Cs = linspecer(1);
    mk = 6;
    iSel = strcmp(SimInfo(:,1),'Jordan')';
    plot([0.75 1.1], [0.75 1.1],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    plot(Table(iSel,strcmp(Headers,'StrideFreq')),...
        Table(iSel,strcmp(Headers,'Exp_StrideFreq')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    xlabel('stride frequency','Interpreter','tex')
    ylabel('stride frequency','Interpreter','tex')
    set(gca,'box','off');
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    title('Jordan 2006 step frequency in various walking speeds');

end


%% Figures paper --

% select a color for each unique experiment we try to reproduce
Cols = linspecer(12);
ColBrown = Cols(1,:);
ColSchertzer = Cols(2,:);
ColGom = Cols(3,:);
ColHuang = Cols(5,:);
ColKoelewijn = Cols(7,:);
ColJordan= Cols(6,:);
ColAbe= Cols(8,:);
ColMcDonald= Cols(9,:);
ColUmberger = Cols(10,:);
ColVanDerZee = Cols(11,:);
ColStrutzenberger = Cols(12,:);


if Settings.PlotFiguresPaper
    %% Figure 1: Plot change in metabolic energy in various conditions
    %   (always w.r.t. reference condition)
    figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic [defaultV]');
    mk = 7;
    nc = 6;
    nr = 2;
    subplot(2,6,4:6);
    if ~Settings.DeltaValues
        plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    %     plot([-1 8], [-1 8]*0.75,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.25,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*0.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %xHeader = 'Exp_PNetMetab_BM';
    %yHeader = 'PNetMetab_b_BM';
    xHeader = 'Exp_PNetMetab_BM';
    yHeader = 'PNetMetab_b_BM';

    %     yHeader = 'PNetMetab_Marg_BM';

    % ------ Browning ------
    iSel = strcmp(SimInfo(:,1),'Browning');
    iRef = strcmp(SimInfo(:,1),'Browning') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0;
    if Settings.DeltaValues
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
    else
        DeltaBrowning = Table(iSel,:);
    end
    Cs = ColBrown;
    l = plot(DeltaBrowning(:,strcmp(Headers,xHeader)),...
        DeltaBrowning(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(1) = l(1); legN{1} = SimInfo{iRef};

    % ---- Schertzer -----
    iSel = strcmp(SimInfo(:,1),'Schertzer');
    iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) >= 1.1 &...
        Table(:,strcmp(Headers,'speed')) <= 1.2);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColSchertzer;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(2) = l(1); legN{2} = SimInfo{iRef};

    % ---- Gomenuka -----
    iSel = strcmp(SimInfo(:,1),'Gomenuka');
    iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.1 &...
        Table(:,strcmp(Headers,'speed')) < 1.2) ;
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColGom;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(3) = l(1); legN{3} = SimInfo{iRef};

    % ---- Huang -----
    iSel = strcmp(SimInfo(:,1),'Huang');
    iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.3);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColHuang;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(4) = l(1); legN{4} = SimInfo{iRef};

    %  ------  Koelewijn ------
    iSel = strcmp(SimInfo(:,1),'Koelewijn');
    iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.4);
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColKoelewijn;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(5) = l(1); legN{5} = SimInfo{iRef};

    %  ------  Umberger ------
    iSel = strcmp(SimInfo(:,1),'Umberger');
    iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
        Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColUmberger;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(6) = l(1); legN{6} = SimInfo{iRef};

    %  ------  Jordan ------ (only spatio-temporal data)

    % ---- Abe -----
    iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
        Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
    iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.05 &...
        Table(:,strcmp(Headers,'speed')) < 1.15) ;
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColAbe;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(7) = l(1); legN{7} = SimInfo{iRef};

    % ---- McDonald -----
    iSel = strcmp(SimInfo(:,1),'McDonald');
    iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColMcDonald;
    l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        DeltaDat(:,strcmp(Headers,yHeader)),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(8) = l(1); legN{8} = SimInfo{iRef};


    % legend and layout
    legend(leg,legN);
    if ~Settings.DeltaValues
        xlabel('Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('Pnet Sim. (W/kg)','Interpreter','tex')
    else
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    % part 2 figure 1Plot change in metabolic energy for each type of intervention

    %-------------------------------------------
    %           effect of walking speed
    %-------------------------------------------
    %   Schertzer
    %   Koelewijn
    %   Abe 2015
    %   Gomenuka
    %
    subplot(nr,nc,7:8);
    if ~Settings.DeltaValues
        plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('speed')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('Pnet Sim. (W/kg)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);



    %------------------------------
    % ---- Loaded Walking ---------
    %------------------------------
    subplot(nr,nc,9:10);
    if ~Settings.DeltaValues
        plot([0 10], [0 10],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-1 3], [-1 3],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    % --- Huang -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColHuang,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Huang'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColBrown,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

    % --- Schertzer -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;



    legend(leg,legN);
    title('loaded')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('Pnet Sim. (W/kg)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    %-------------------------------
    %----- Walking on a slope ------
    %-------------------------------
    subplot(nr,nc,11:12);
    if ~Settings.DeltaValues
        plot([0 20], [0 20],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 6], [-5 6],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColMcDonald,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'McDonald'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,Settings.DeltaValues);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('slope')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('Pnet Sim. (W/kg)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    SetDescriptiveStatistics_Fig(gcf,'paper_predint');

    %% Figure 2:  Plot change in step frequency for each type of intervention

    BoolPlotDelta = Settings.DeltaValues;
    figure('Color',[1 1 1],'Name','Figure step frequency type intervention [DefaultV]');
    nc = 3;
    nr = 1;
    xHeader = 'Exp_StrideFreq';
    yHeader = 'StrideFreq';

    %-------------------------------------------
    %           effect of walking speed
    %-------------------------------------------
    %   Schertzer
    %   Koelewijn
    %   Abe 2015
    %   Gomenuka
    %
    subplot(nr,nc,1);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Jordan','speed');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColJordan,mk,xHeader,yHeader,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Jordan'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','speed');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'VanDerZee2022','speed');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColVanDerZee,mk,xHeader,yHeader,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'VanDerZee'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('speed')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);



    %------------------------------
    % ---- Loaded Walking ---------
    %------------------------------
    subplot(nr,nc,2);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColBrown,mk,xHeader,yHeader,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('loaded')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    %-------------------------------
    %----- Walking on a slope ------
    %-------------------------------
    subplot(nr,nc,3);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;
    %
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
    l = PlotIndexPlot(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('slope')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    SetDescriptiveStatistics_Fig(gcf,'paper_predint');


    %% Figure S1: effect energy equation on prediced change in metabolic energy

    % effect of energy equation on results change in metabolic energy:
    yHeaders = {'PNetMetab_b_BM','PNetMetab_Marg_BM'};
    %     yHeaders = {'PNetMetab_BM','PNetMetab_b_BM'};
    Cols = linspecer(2);
    figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic dependent on energy eq [defaultV]');
    mk = 7;
    nc = 6;
    nr = 2;
    for ih = 1:2
        Cs = Cols(ih,:);

        subplot(2,6,4:6);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        xHeader = 'Exp_PNetMetab_BM';
        yHeader = yHeaders{ih};
        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaBrowning(:,strcmp(Headers,xHeader)),...
            DeltaBrowning(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Umberger ------
        iSel = strcmp(SimInfo(:,1),'Umberger');
        iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
            Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % legend and layout
        if ih ==2
            xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
            set(gca,'box','off')
            set(gca,'FontSize',14);
            set(gca,'LineWidth',1.5);
        end

        % part 2 figure 1Plot change in metabolic energy for each type of intervention

        %-------------------------------------------
        %           effect of walking speed
        %-------------------------------------------
        subplot(nr,nc,7:8);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        title('speed')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %------------------------------
        % ---- Loaded Walking ---------
        %------------------------------
        subplot(nr,nc,9:10);
        plot([-1 3], [-1 3],'--','Color',[0 0 0],'LineWidth',1.3); hold on;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        % --- Schertzer -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        % --- Huang -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);
        title('loaded')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %-------------------------------
        %----- Walking on a slope ------
        %-------------------------------
        subplot(nr,nc,11:12);
        plot([-5 6], [-5 6],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader);
        title('slope')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);
    end
    SetDescriptiveStatistics_Fig(gcf,'paper_predint');

    %% Figure S1: only margaria energy model

    % effect of energy equation on results change in metabolic energy:
    yHeaders = {'PNetMetab_Marg_BM'};
    %     yHeaders = {'PNetMetab_BM','PNetMetab_b_BM'};
    Cols = linspecer(2);
    figure('Color',[1 1 1],'Name','Figure Delate Pmetab Marg [defaultV]');
    mk = 7;
    nc = 6;
    nr = 2;
    for ih = 1
        Cs = Cols(ih,:);

        subplot(2,6,4:6);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        xHeader = 'Exp_PNetMetab_BM';
        yHeader = yHeaders{ih};
        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        if Settings.DeltaValues
            DeltaBrowning = Table(iSel,:) - Table(iRef,:);
        else
            DeltaBrowning = Table(iSel,:);
        end
        l = plot(DeltaBrowning(:,strcmp(Headers,xHeader)),...
            DeltaBrowning(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Umberger ------
        iSel = strcmp(SimInfo(:,1),'Umberger');
        iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
            Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerSize',mk);

        % legend and layout
        if ih ==2
            xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
            set(gca,'box','off')
            set(gca,'FontSize',14);
            set(gca,'LineWidth',1.5);
        end

        % part 2 figure 1Plot change in metabolic energy for each type of intervention

        %-------------------------------------------
        %           effect of walking speed
        %-------------------------------------------
        subplot(nr,nc,7:8);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        title('speed')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %------------------------------
        % ---- Loaded Walking ---------
        %------------------------------
        subplot(nr,nc,9:10);
        plot([-1 3], [-1 3],'--','Color',[0 0 0],'LineWidth',1.3); hold on;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        % --- Schertzer -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        % --- Huang -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);
        title('loaded')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %-------------------------------
        %----- Walking on a slope ------
        %-------------------------------
        subplot(nr,nc,11:12);
        plot([-5 6], [-5 6],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
        l = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,Settings.DeltaValues);
        title('slope')
        xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);
    end
    SetDescriptiveStatistics_Fig(gcf,'paper_predint');


    %% Plot effect of ignoring tendon compliance and/or force-length-velocity properties

    figure('Color',[1 1 1],'Name','Figure S2: effect ignoring muscle dynamics [defaultV]');
    nr = 1; nc = 2;
    Cols = linspecer(4);

    % plot effect on metabolic power
    subplot(1,2,1)
    yHeader = 'COT_b';
    xHeader = 'speed';
    ModelTypeNames = {'RigidTendon','RigidTendon_IgnoreFLV','Default'};
    leg = [];    legN = []; ctLeg = 1;
    for i = 1:length(ModelTypeNames)
        iSel = strcmp(SimInfo(:,1),'WalkSpeedRef') & ...
            strcmp(SimInfo(:,2),ModelTypeNames{i});
        Cs = Cols(i,:);
        l = plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
        leg(ctLeg) = l(1); legN{ctLeg} = ModelTypeNames{i}; ctLeg = ctLeg + 1;
    end
    [DatAbe, HeaderAbe] = GetDataAbe2015();
    iSel = DatAbe(:,strcmp(HeaderAbe,'grade')) == 0;
    Cs = [0.6 0.6 0.6];
    l = plot(DatAbe(iSel,strcmp(HeaderAbe,'speed')),...
        DatAbe(iSel,strcmp(HeaderAbe,'Metab_COT')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;
    legend(leg,legN,'interpreter','none');
    xlabel('walking speed [m/s]','Interpreter','tex')
    ylabel('COT [J/kg/m]','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);

    % effect on stride frequency
    subplot(1,2,2)
    yHeader = 'StrideFreq';
    xHeader = 'speed';
    ModelTypeNames = {'RigidTendon','RigidTendon_IgnoreFLV','Default'};
    leg = [];    legN = []; ctLeg = 1;
    for i = 1:length(ModelTypeNames)
        iSel = strcmp(SimInfo(:,1),'WalkSpeedRef') & ...
            strcmp(SimInfo(:,2),ModelTypeNames{i});
        Cs = Cols(i,:);
        l = plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
        leg(ctLeg) = l(1); legN{ctLeg} = ModelTypeNames{i}; ctLeg = ctLeg + 1;
    end
    % experimental data Jordan
    DExp = importdata(dExpJordan);
    PreferedSpeed = 1.3;
    ExpSpeeds = DExp(:,1).*PreferedSpeed;
    ExpStrideFreq = ExpSpeeds./DExp(:,2);
    Cs = [0.6 0.6 0.6];
    l = plot(ExpSpeeds,ExpStrideFreq,'o','Color',Cs,'MarkerFaceColor',...
        Cs,'MarkerSize',mk); hold on;
    leg(ctLeg) = l(1); legN{ctLeg} = 'Jordan'; ctLeg = ctLeg + 1;
    legend(leg,legN,'interpreter','none');
    xlabel('walking speed [m/s]','Interpreter','tex')
    ylabel('stride frequency','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);


    %% Plot effect weight objective function : as a function of walking speed
    figure('Color',[1 1 1],'Name','Figure S3: effect weight objective function [defaultV]');
    nr = 1; nc = 2;
    Cols = linspecer(4);

    % plot effect on metabolic power
    subplot(1,2,1)
    yHeader = 'COT_b';
    xHeader = 'speed';
    ModelTypeNames = {'lower_w_metab','lower_w_qdd','Default','Default_QR'};
    leg = [];    legN = []; ctLeg = 1;
    for i = 1:length(ModelTypeNames)
        iSel = strcmp(SimInfo(:,1),'WalkSpeedRef_Obj') & ...
            strcmp(SimInfo(:,2),ModelTypeNames{i});
        Cs = Cols(i,:);
        l = plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
        plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            '--','Color',Cs); hold on;
        leg(ctLeg) = l(1); legN{ctLeg} = ModelTypeNames{i}; ctLeg = ctLeg + 1;
    end
    [DatAbe, HeaderAbe] = GetDataAbe2015();
    iSel = DatAbe(:,strcmp(HeaderAbe,'grade')) == 0;
    Cs = [0.6 0.6 0.6];
    l = plot(DatAbe(iSel,strcmp(HeaderAbe,'speed')),...
        DatAbe(iSel,strcmp(HeaderAbe,'Metab_COT')),...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    plot(DatAbe(iSel,strcmp(HeaderAbe,'speed')),...
        DatAbe(iSel,strcmp(HeaderAbe,'Metab_COT')),...
        '--','Color',Cs); hold on;
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;
    legend(leg,legN,'interpreter','none');
    xlabel('walking speed [m/s]','Interpreter','tex')
    ylabel('COT [J/kg/m]','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);

    % effect on stride frequency
    subplot(1,2,2)
    yHeader = 'StrideFreq';
    xHeader = 'speed';
    ModelTypeNames = {'lower_w_metab','lower_w_qdd','Default','Default_QR'};
    leg = [];    legN = []; ctLeg = 1;
    for i = 1:length(ModelTypeNames)
        iSel = strcmp(SimInfo(:,1),'WalkSpeedRef_Obj') & ...
            strcmp(SimInfo(:,2),ModelTypeNames{i});
        Cs = Cols(i,:);
        l = plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
        plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            '--','Color',Cs); hold on;
        leg(ctLeg) = l(1); legN{ctLeg} = ModelTypeNames{i}; ctLeg = ctLeg + 1;
    end
    % experimental data Jordan
    DExp = importdata(dExpJordan);
    PreferedSpeed = 1.3;
    ExpSpeeds = DExp(:,1).*PreferedSpeed;
    ExpStrideFreq = ExpSpeeds./DExp(:,2);
    Cs = [0.6 0.6 0.6];
    l = plot(ExpSpeeds,ExpStrideFreq,'o','Color',Cs,'MarkerFaceColor',...
        Cs,'MarkerSize',mk); hold on;
    plot(ExpSpeeds,ExpStrideFreq,'--','Color',Cs); hold on;
    leg(ctLeg) = l(1); legN{ctLeg} = 'Jordan'; ctLeg = ctLeg + 1;
    legend(leg,legN,'interpreter','none');
    xlabel('walking speed [m/s]','Interpreter','tex')
    ylabel('stride frequency','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1.5);


    %% Plot error in metabolic energy prediction as a function of mechanical
    % non-dim units and converted to osim model

    %   (always w.r.t. reference condition)
    figure('Color',[1 1 1],'Name','Figure Error Pmetab [nondim]');
    mk = 7;
    nc = 5;
    nr = 2;
    xHeader = 'Exp_PNetMetab_BM';
    yHeader = 'PNetMetab_b_BM';
    %horHeaders = {'NetMechWork','PosMuscleWork','NegMuscleWork','mean_act'};
    horHeaders = {'NetMechPower','NetPosMusclePower','NetNegMusclePower','mean_act'};
    horHeadersLabels = {'Net. mech. power','pos. muscle power','neg. muscle power','mean muscle act'};
    for ih = 1:length(horHeaders)
        horHeader = horHeaders{ih};
        subplot(1,4,ih);

        %     yHeader = 'PNetMetab_Marg_BM';

        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColBrown;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;
        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(1) = l(1); legN{1} = SimInfo{iRef};

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColSchertzer;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(2) = l(1); legN{2} = SimInfo{iRef};

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColGom;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(3) = l(1); legN{3} = SimInfo{iRef};

        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColHuang;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));

        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(4) = l(1); legN{4} = SimInfo{iRef};

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColKoelewijn;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(5) = l(1); legN{5} = SimInfo{iRef};

        % %  ------  Umberger ------
        % iSel = strcmp(SimInfo(:,1),'Umberger');
        % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
        %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
        %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        % if  Settings.DeltaValues
        %     DeltaDat = Table(iSel,:) - Table(iRef,:);
        % else
        %     DeltaDat = Table(iSel,:);
        % end
        % Cs = ColUmberger;
        % xdat = DeltaDat(:,strcmp(Headers,xHeader));
        % ydat = DeltaDat(:,strcmp(Headers,yHeader));
        % Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Umberger'));
        % mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Umberger'));
        % scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        % scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        % xdat = xdat./scale_exp.*scale_sim;
        % l = plot(Table(iSel,strcmp(Headers,horHeader)),(ydat - xdat)*Mass_Sim,...
        %     'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        % leg(6) = l(1); legN{6} = SimInfo{iRef};

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColAbe;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(7) = l(1); legN{7} = SimInfo{iRef};

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColMcDonald;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        dffPlot= (ydat - xdat)*Mass_Sim;

        l = plot(DeltaDat(:,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(8) = l(1); legN{8} = SimInfo{iRef};


        % legend and layout

        % legend(leg,legN);
        if ih == 1
            ylabel('Error relative metab prediction (W/kg)','Interpreter','tex')
        end
        xlabel(horHeadersLabels{ih},'Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',12);
        set(gca,'LineWidth',1.5);
        % compute correlation cofficients
        ah = gca;
        XData = [ah.Children.XData];
        YData = [ah.Children.YData];
        iNan = find(isnan(XData) | isnan(YData));
        XData(iNan) = [];
        YData(iNan)= [];
        stats = regstats(XData',YData','linear');
        title('Rsq ', num2str(stats.rsquare));


    end
    subplot(1,4,3)
    % set(gca,'XLim',[-140 0])

    %% Mechanical efficiency of muscle fibers

    figure('Color',[1 1 1],'Name','Figure mechanical efficiency');
    mk = 7;
    nc = 6;
    nr = 2;
    SimConditions = {'Abe2015','Browning','Gomenuka','Huang','Jordan','Koelewijn','McDonald','Schertzer','VanDerZee2022'};
    VarSpeed = [1 ]

    % ------ Browning ------
    iSel = strcmp(SimInfo(:,1),'Browning');
    Cs = ColBrown;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    subplot(2,3,2)
    plot(AddedMass, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,5)
    plot(AddedMass, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    iSel = strcmp(SimInfo(:,1),'Huang');
    Cs = ColHuang;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    subplot(2,3,2)
    plot(AddedMass, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,5)
    plot(AddedMass, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    iSel = strcmp(SimInfo(:,1),'Schertzer') & round(Table(:,strcmp(Headers,'speed')),2) == 1.11;
    Cs = ColSchertzer
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    subplot(2,3,2)
    plot(AddedMass, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,5)
    plot(AddedMass, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    % walking speed
    iSel = strcmp(SimInfo(:,1),'Jordan');
    Cs = ColJordan;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    walkspeed = Table(iSel,strcmp(Headers,'speed'));
    subplot(2,3,1)
    plot(walkspeed, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,4)
    plot(walkspeed, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    % walking slope
    iSel = strcmp(SimInfo(:,1),'Gomenuka');
    Cs = ColGom;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    walkspeed = Table(iSel,strcmp(Headers,'speed'));
    slope = Table(iSel,strcmp(Headers,'Slope'));
    subplot(2,3,3)
    plot(slope, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,6)
    plot(slope, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    % walking slope
    iSel = strcmp(SimInfo(:,1),'Abe2015') & Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
    Cs = ColAbe;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    walkspeed = Table(iSel,strcmp(Headers,'speed'));
    slope = Table(iSel,strcmp(Headers,'Slope'));
    subplot(2,3,3)
    plot(slope, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,6)
    plot(slope, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    % walking slope
    iSel = strcmp(SimInfo(:,1),'McDonald');
    Cs = ColMcDonald;
    NetMechWork = Table(iSel,strcmp(Headers,'NetMechPower'));
    PosMechWork = Table(iSel,strcmp(Headers,'NetPosMusclePower'));
    PNetMetab = Table(iSel,strcmp(Headers,'NetMetabPower_nobasal'));
    AddedMass = Table(iSel,strcmp(Headers,'AddedMass'));
    walkspeed = Table(iSel,strcmp(Headers,'speed'));
    slope = Table(iSel,strcmp(Headers,'Slope'));
    subplot(2,3,3)
    plot(slope, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
    subplot(2,3,6)
    plot(slope, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;

    for i = 1:6
        subplot(2,3,i)
        if i < 4
            set(gca,'YLim',[0 0.5]);
        else
            set(gca,'YLim',[0 0.4]);
        end
        set(gca,'box','off')
        set(gca,'FontSize',12);
        set(gca,'LineWidth',1.5);
    end

    %% alternative for figure above
    figure('Color',[1 1 1],'Name','Figure mechanical efficiency 2');
    nSim = length(Table);
    for i =1:nSim
        SimCond = SimInfo{i,1};
        BoolPlot = true;
        if strcmp(SimCond,'Abe2015')
            Cs = ColAbe;
        elseif strcmp(SimCond,'Browning')
            Cs = ColBrown;
        elseif strcmp(SimCond,'Gomenuka')
            Cs = ColGom;
        elseif strcmp(SimCond,'Huang')
            Cs = ColHuang;
        elseif strcmp(SimCond,'Jordan')
            Cs = ColJordan;
        elseif strcmp(SimCond,'Koelewijn')
            Cs = ColKoelewijn;
        elseif strcmp(SimCond,'McDonald')
            Cs = ColMcDonald;
        elseif strcmp(SimCond,'Schertzer')
            Cs = ColSchertzer;
        elseif strcmp(SimCond,'VanDerZee2022')
            Cs = ColVanDerZee;
        else
            BoolPlot = false;
        end
        % ,'NetMechPower','NetPosMusclePower',...
        % 'NetNegMusclePower','NetMetabPower_nobasal'};
        if BoolPlot
            NetMechWork = Table(i,strcmp(Headers,'NetMechPower'));
            PosMechWork = Table(i,strcmp(Headers,'NetPosMusclePower'));
            PNetMetab = Table(i,strcmp(Headers,'NetMetabPower_nobasal'));
            subplot(1,2,1)
            plot(PosMechWork, PosMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
            subplot(1,2,2)
            plot(NetMechWork, NetMechWork./PNetMetab,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk); hold on;
        end
    end
    for i = 1:2
        subplot(1,2,i)
        set(gca,'box','off')
        set(gca,'FontSize',12);
        set(gca,'LineWidth',1.5);
        if i ==1
            set(gca,'YLim',[0.1 0.5]);
        end
        if i == 2
            set(gca,'YLim',[-0.2 0.5]);
        end
    end


    %% Figure 2:  Plot change in step frequency for each type of intervention -- nondim

    BoolPlotDelta = Settings.DeltaValues;
    figure('Color',[1 1 1],'Name','Figure step frequency type intervention [nondim]');
    nc = 3;
    nr = 1;
    xHeader = 'Exp_StrideFreq';
    yHeader = 'StrideFreq';

    %-------------------------------------------
    %           effect of walking speed
    %-------------------------------------------
    %   Schertzer
    %   Koelewijn
    %   Abe 2015
    %   Gomenuka
    %
    subplot(nr,nc,1);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    % IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Jordan','speed');
    % Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Jordan'));
    % scale_exp =  sqrt(9.81/Lexp);
    % scale_sim = sqrt(9.81/LegLength_Sim);
    % l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColJordan,mk,xHeader,yHeader,...
    %     scale_exp, scale_sim,BoolPlotDelta);
    % leg(ctLeg) = l(1); legN{ctLeg} = 'Jordan'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  sqrt(9.81/Lexp);
    scale_sim = sqrt(9.81/LegLength_Sim);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,...
        scale_exp, scale_sim,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'VanDerZee2022','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'VanDerZee2022'));
    scale_exp =  sqrt(9.81/Lexp);
    scale_sim = sqrt(9.81/LegLength_Sim);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColVanDerZee,mk,xHeader,yHeader,...
        scale_exp, scale_sim,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'VanDerZee'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('speed')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);



    %------------------------------
    % ---- Loaded Walking ---------
    %------------------------------
    subplot(nr,nc,2);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
    scale_exp =  sqrt(9.81/Lexp);
    scale_sim = sqrt(9.81/LegLength_Sim);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColBrown,mk,xHeader,yHeader,...
        scale_exp, scale_sim,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('loaded')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    %-------------------------------
    %----- Walking on a slope ------
    %-------------------------------
    subplot(nr,nc,3);
    if BoolPlotDelta
        plot([-0.5 0.5], [-0.5 0.5],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([0.5 1.2], [0.5 1.2],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;
    %
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  sqrt(9.81/Lexp);
    scale_sim = sqrt(9.81/LegLength_Sim);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,...
        scale_exp,scale_sim,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Strutzenberger','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Strutzenberger'));
    scale_exp =  sqrt(9.81/Lexp);
    scale_sim = sqrt(9.81/LegLength_Sim);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColStrutzenberger,mk,xHeader,yHeader,...
        scale_exp,scale_sim,BoolPlotDelta);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Struzenberger'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('slope')
    xlabel('\Delta step frequency (Hz)','Interpreter','tex')
    ylabel('\Delta step frequency (Hz)','Interpreter','tex')
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);
    SetDescriptiveStatistics_Fig(gcf,'paper_predint');

    %% Figure metabolic energy consumption, nondim

    %   (always w.r.t. reference condition)
    figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic [nondim]',...
        'Position',[0.2250    0.2490    1.1167    0.6580]*10^3);

    mk = 7;
    nc = 6;
    nr = 2;
    subplot(2,6,4:6);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 11]*Mass_Sim, [-5 11]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    %     plot([-1 8], [-1 8]*0.75,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.25,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*0.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %xHeader = 'Exp_PNetMetab_BM';
    %yHeader = 'PNetMetab_b_BM';
    xHeader = 'Exp_PNetMetab_BM';
    yHeader = 'PNetMetab_b_BM';

    %     yHeader = 'PNetMetab_Marg_BM';

    % ------ Browning ------
    iSel = strcmp(SimInfo(:,1),'Browning');
    iRef = strcmp(SimInfo(:,1),'Browning') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0;
    if Settings.DeltaValues
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
    else
        DeltaBrowning = Table(iSel,:);
    end
    Cs = ColBrown;

    xdat = DeltaBrowning(:,strcmp(Headers,xHeader));
    ydat = DeltaBrowning(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;

    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(1) = l(1); legN{1} = SimInfo{iRef};

    % ---- Schertzer -----
    iSel = strcmp(SimInfo(:,1),'Schertzer');
    iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) >= 1.1 &...
        Table(:,strcmp(Headers,'speed')) <= 1.2);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColSchertzer;

    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(2) = l(1); legN{2} = SimInfo{iRef};

    % ---- Gomenuka -----
    iSel = strcmp(SimInfo(:,1),'Gomenuka');
    iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.1 &...
        Table(:,strcmp(Headers,'speed')) < 1.2) ;
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColGom;

    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(3) = l(1); legN{3} = SimInfo{iRef};

    % ---- Huang -----
    iSel = strcmp(SimInfo(:,1),'Huang');
    iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.3);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColHuang;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(4) = l(1); legN{4} = SimInfo{iRef};

    %  ------  Koelewijn ------
    iSel = strcmp(SimInfo(:,1),'Koelewijn');
    iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.4);
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColKoelewijn;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(5) = l(1); legN{5} = SimInfo{iRef};

    % %  ------  Umberger ------
    % iSel = strcmp(SimInfo(:,1),'Umberger');
    % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
    %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
    %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
    %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
    %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
    % if  Settings.DeltaValues
    %     DeltaDat = Table(iSel,:) - Table(iRef,:);
    % else
    %     DeltaDat = Table(iSel,:);
    % end
    % Cs = ColUmberger;
    % xdat = DeltaDat(:,strcmp(Headers,xHeader));
    % ydat = DeltaDat(:,strcmp(Headers,yHeader));
    % Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Umberger'));
    % mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Umberger'));
    % scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    % scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    % xdat = xdat./scale_exp.*scale_sim;
    % l = plot(xdat,ydat,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    % leg(6) = l(1); legN{6} = SimInfo{iRef};

    %  ------  Jordan ------ (only spatio-temporal data)

    % ---- Abe -----
    iSel = strcmp(SimInfo(:,1),'Abe2015') &...
        Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
    iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.05 &...
        Table(:,strcmp(Headers,'speed')) < 1.15) ;
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColAbe;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(6) = l(1); legN{6} = SimInfo{iRef};

    % ---- McDonald -----
    iSel = strcmp(SimInfo(:,1),'McDonald');
    iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColMcDonald;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(7) = l(1); legN{7} = SimInfo{iRef};


    % legend and layout
    legend(leg,legN);
    if ~Settings.DeltaValues
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    % part 2 figure 1Plot change in metabolic energy for each type of intervention

    %-------------------------------------------
    %           effect of walking speed
    %-------------------------------------------
    %   Schertzer
    %   Koelewijn
    %   Abe 2015
    %   Gomenuka
    %
    subplot(nr,nc,7:8);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 8]*Mass_Sim, [-5 8]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('speed')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);



    %------------------------------
    % ---- Loaded Walking ---------
    %------------------------------
    subplot(nr,nc,9:10);
    if ~Settings.DeltaValues
        plot([0 10]*Mass_Sim, [0 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-1 3]*Mass_Sim, [-1 3]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColBrown,mk,xHeader,yHeader,...
        scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

    % --- Schertzer -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,...
        scale_exp, scale_sim, Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

    % --- Huang -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColHuang,mk,xHeader,yHeader,...
        scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Huang'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('loaded')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    %-------------------------------
    %----- Walking on a slope ------
    %-------------------------------
    subplot(nr,nc,11:12);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 6]*Mass_Sim, [-5 6]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColMcDonald,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'McDonald'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    legend(leg,legN);
    title('slope')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    SetDescriptiveStatistics_Fig(gcf,'paper_predint');

    %% dimensieloos margaria model


    %   (always w.r.t. reference condition)
    figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic - Marg [nondim]');
    mk = 7;
    nc = 6;
    nr = 2;
    subplot(2,6,4:6);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-6 18]*Mass_Sim, [-6 18]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    %     plot([-1 8], [-1 8]*0.75,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.25,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*0.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %     plot([-1 8], [-1 8]*1.5,'--','Color',[0.2 0.2 0.2],'LineWidth',0.8); hold on;
    %xHeader = 'Exp_PNetMetab_BM';
    %yHeader = 'PNetMetab_b_BM';
    xHeader = 'Exp_PNetMetab_BM';
    yHeader = 'PNetMetab_Marg_BM';

    % ------ Browning ------
    iSel = strcmp(SimInfo(:,1),'Browning');
    iRef = strcmp(SimInfo(:,1),'Browning') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0;
    if Settings.DeltaValues
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
    else
        DeltaBrowning = Table(iSel,:);
    end
    Cs = ColBrown;

    xdat = DeltaBrowning(:,strcmp(Headers,xHeader));
    ydat = DeltaBrowning(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;

    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,...
        'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(1) = l(1); legN{1} = SimInfo{iRef};

    % ---- Schertzer -----
    iSel = strcmp(SimInfo(:,1),'Schertzer');
    iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) >= 1.1 &...
        Table(:,strcmp(Headers,'speed')) <= 1.2);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColSchertzer;

    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(2) = l(1); legN{2} = SimInfo{iRef};

    % ---- Gomenuka -----
    iSel = strcmp(SimInfo(:,1),'Gomenuka');
    iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.1 &...
        Table(:,strcmp(Headers,'speed')) < 1.2) ;
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColGom;

    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(3) = l(1); legN{3} = SimInfo{iRef};

    % ---- Huang -----
    iSel = strcmp(SimInfo(:,1),'Huang');
    iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.3);
    if Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColHuang;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(4) = l(1); legN{4} = SimInfo{iRef};

    %  ------  Koelewijn ------
    iSel = strcmp(SimInfo(:,1),'Koelewijn');
    iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.2 &...
        Table(:,strcmp(Headers,'speed')) < 1.4);
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColKoelewijn;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(5) = l(1); legN{5} = SimInfo{iRef};

    % %  ------  Umberger ------
    % iSel = strcmp(SimInfo(:,1),'Umberger');
    % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
    %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
    %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
    %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
    %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
    % if  Settings.DeltaValues
    %     DeltaDat = Table(iSel,:) - Table(iRef,:);
    % else
    %     DeltaDat = Table(iSel,:);
    % end
    % Cs = ColUmberger;
    % xdat = DeltaDat(:,strcmp(Headers,xHeader));
    % ydat = DeltaDat(:,strcmp(Headers,yHeader));
    % Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Umberger'));
    % mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Umberger'));
    % scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    % scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    % xdat = xdat./scale_exp.*scale_sim;
    % l = plot(xdat,ydat,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    % leg(6) = l(1); legN{6} = SimInfo{iRef};

    %  ------  Jordan ------ (only spatio-temporal data)

    % ---- Abe -----
    iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
        Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
    iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,'speed')) > 1.05 &...
        Table(:,strcmp(Headers,'speed')) < 1.15) ;
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColAbe;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(6) = l(1); legN{6} = SimInfo{iRef};

    % ---- McDonald -----
    iSel = strcmp(SimInfo(:,1),'McDonald');
    iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
        Table(:,strcmp(Headers,'Slope')) == 0 & ...
        Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
    if  Settings.DeltaValues
        DeltaDat = Table(iSel,:) - Table(iRef,:);
    else
        DeltaDat = Table(iSel,:);
    end
    Cs = ColMcDonald;
    xdat = DeltaDat(:,strcmp(Headers,xHeader));
    ydat = DeltaDat(:,strcmp(Headers,yHeader));
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    xdat = xdat./scale_exp.*scale_sim;
    l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    leg(7) = l(1); legN{7} = SimInfo{iRef};


    % legend and layout
    % legend(leg,legN);
    if ~Settings.DeltaValues
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    % part 2 figure 1Plot change in metabolic energy for each type of intervention

    %-------------------------------------------
    %           effect of walking speed
    %-------------------------------------------
    %   Schertzer
    %   Koelewijn
    %   Abe 2015
    %   Gomenuka
    %
    subplot(nr,nc,7:8);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-6 18]*Mass_Sim, [-6 18]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

    % legend(leg,legN);
    title('speed')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);



    %------------------------------
    % ---- Loaded Walking ---------
    %------------------------------
    subplot(nr,nc,9:10);
    if ~Settings.DeltaValues
        plot([0 10]*Mass_Sim, [0 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-1 3]*Mass_Sim, [-1 3]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColBrown,mk,xHeader,yHeader,...
        scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

    % --- Schertzer -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColSchertzer,mk,xHeader,yHeader,...
        scale_exp, scale_sim, Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

    % --- Huang -----
    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColHuang,mk,xHeader,yHeader,...
        scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Huang'; ctLeg = ctLeg + 1;

    % legend(leg,legN);
    title('loaded')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    %-------------------------------
    %----- Walking on a slope ------
    %-------------------------------
    subplot(nr,nc,11:12);
    if ~Settings.DeltaValues
        plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    else
        plot([-5 10]*Mass_Sim, [-5 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
    end
    leg = [];    legN = []; ctLeg = 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColKoelewijn,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColGom,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColMcDonald,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'McDonald'; ctLeg = ctLeg + 1;

    IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
    Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
    mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
    scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
    scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
    l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,ColAbe,mk,xHeader,yHeader,...
        scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
    leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

    % legend(leg,legN);
    title('slope')
    if Settings.DeltaValues
        xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
        ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
    else
        xlabel('Pnet Exp. (W)','Interpreter','tex')
        ylabel('Pnet Sim. (W)','Interpreter','tex')
    end
    set(gca,'box','off')
    set(gca,'FontSize',14);
    set(gca,'LineWidth',1.5);

    SetDescriptiveStatistics_Fig(gcf,'paper_predint');

    %% Plot error in metabolic energy prediction as a function of mechanical
    % dim units and converted to osim model

    %   (always w.r.t. reference condition)
    figure('Color',[1 1 1],'Name','Figure Error Pmetabolic [defaultV]');
    mk = 7;
    nc = 5;
    nr = 2;
    xHeader = 'Exp_PNetMetab_BM';
    yHeader = 'PNetMetab_b_BM';
    % horHeaders = {'NetMechWork','PosMuscleWork','NegMuscleWork','mean_act'};
    horHeadersLabels = {'Net. mech. power','pos. muscle power','neg. muscle power','mean muscle act'};
    horHeaders = {'NetMechPower','NetPosMusclePower','NetNegMusclePower','mean_act'};
    for ih = 1:length(horHeaders)
        horHeader = horHeaders{ih};
        subplot(1,4,ih);

        %     yHeader = 'PNetMetab_Marg_BM';

        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColBrown;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(1) = l(1); legN{1} = SimInfo{iRef};

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColSchertzer;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(2) = l(1); legN{2} = SimInfo{iRef};

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColGom;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(3) = l(1); legN{3} = SimInfo{iRef};

        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        if Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColHuang;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(4) = l(1); legN{4} = SimInfo{iRef};

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColKoelewijn;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(5) = l(1); legN{5} = SimInfo{iRef};

        % %  ------  Umberger ------
        % iSel = strcmp(SimInfo(:,1),'Umberger');
        % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
        %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
        %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        % if  Settings.DeltaValues
        %     DeltaDat = Table(iSel,:) - Table(iRef,:);
        % else
        %     DeltaDat = Table(iSel,:);
        % end
        % Cs = ColUmberger;
        % xdat = DeltaDat(:,strcmp(Headers,xHeader));
        % ydat = DeltaDat(:,strcmp(Headers,yHeader));
        % Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Umberger'));
        % mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Umberger'));
        % scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        % scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        % xdat = xdat./scale_exp.*scale_sim;
        % l = plot(Table(iSel,strcmp(Headers,horHeader)),(ydat - xdat)*Mass_Sim,...
        %     'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        % leg(6) = l(1); legN{6} = SimInfo{iRef};

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColAbe;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(7) = l(1); legN{7} = SimInfo{iRef};

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        if  Settings.DeltaValues
            DeltaDat = Table(iSel,:) - Table(iRef,:);
        else
            DeltaDat = Table(iSel,:);
        end
        Cs = ColMcDonald;
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        dffPlot= ydat - xdat;
        l = plot(Table(iSel,strcmp(Headers,horHeader)),dffPlot,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);  hold on;
        leg(8) = l(1); legN{8} = SimInfo{iRef};


        % legend and layout

        % legend(leg,legN);
        if ih == 1
            ylabel('Error relative metab prediction (W/kg)','Interpreter','tex')
        end
        xlabel(horHeadersLabels{ih},'Interpreter','tex')
        set(gca,'box','off')
        set(gca,'FontSize',12);
        set(gca,'LineWidth',1.5);
        % compute correlation cofficients
        ah = gca;
        XData = [ah.Children.XData];
        YData = [ah.Children.YData];
        iNan = find(isnan(XData) | isnan(YData));
        XData(iNan) = [];
        YData(iNan)= [];
        stats = regstats(XData',YData','linear');
        title('Rsq ', num2str(stats.rsquare));


    end
    subplot(1,4,3)
    % set(gca,'XLim',[-140 0])


    %% Effect energy equation on delta metabolic [nondim]
    % effect of energy equation on results change in metabolic energy:
    yHeaders = {'PNetMetab_b_BM','PNetMetab_Marg_BM'};
    %     yHeaders = {'PNetMetab_BM','PNetMetab_b_BM'};
    Cols = [88, 101, 218; 178, 22, 22]./255;
    figure('Color',[1 1 1],'Name','Figure Delta Pmetabolic dependent on energy eq [nondim]');
    mk = 7;
    nc = 6;
    nr = 2;
    for ih = 1:2
        Cs = Cols(ih,:);

        subplot(2,6,4:6);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        xHeader = 'Exp_PNetMetab_BM';
        yHeader = yHeaders{ih};
        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaBrowning(:,strcmp(Headers,xHeader));
        ydat = DeltaBrowning(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        DeltaDat = Table(iSel,:) - Table(iRef,:);

        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(2) = l(1); legN{2} = SimInfo{iRef};

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);

        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(3) = l(1); legN{3} = SimInfo{iRef};


        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(4) = l(1); legN{4} = SimInfo{iRef};

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(5) = l(1); legN{5} = SimInfo{iRef};

        % %  ------  Umberger ------
        % iSel = strcmp(SimInfo(:,1),'Umberger');
        % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
        %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
        %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        % DeltaDat = Table(iSel,:) - Table(iRef,:);
        % l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        %     DeltaDat(:,strcmp(Headers,yHeader)),...
        %     'o','Color',Cs,'MarkerSize',mk);

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(6) = l(1); legN{6} = SimInfo{iRef};

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(7) = l(1); legN{7} = SimInfo{iRef};

        % legend and layout
        if ih ==2
            xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
            set(gca,'box','off')
            set(gca,'FontSize',14);
            set(gca,'LineWidth',1.5);
        end

        % part 2 figure 1Plot change in metabolic energy for each type of intervention

        %-------------------------------------------
        %           effect of walking speed
        %-------------------------------------------
        %   Schertzer
        %   Koelewijn
        %   Abe 2015
        %   Gomenuka
        %
        subplot(nr,nc,7:8);
        if ~Settings.DeltaValues
            plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-6 18]*Mass_Sim, [-6 18]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('speed')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);



        %------------------------------
        % ---- Loaded Walking ---------
        %------------------------------
        subplot(nr,nc,9:10);
        if ~Settings.DeltaValues
            plot([0 10]*Mass_Sim, [0 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-1 3]*Mass_Sim, [-1 3]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

        % --- Schertzer -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim, Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

        % --- Huang -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Huang'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('loaded')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %-------------------------------
        %----- Walking on a slope ------
        %-------------------------------
        subplot(nr,nc,11:12);
        if ~Settings.DeltaValues
            plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-5 10]*Mass_Sim, [-5 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'McDonald'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('slope')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

    end


    %% Effect energy equation on delta metabolic [nondim]

    % same figure without delta values ?

    % effect of energy equation on results change in metabolic energy:
    yHeaders = {'PNetMetab_b_BM','margaria_negworkbased'};
    %     yHeaders = {'PNetMetab_BM','PNetMetab_b_BM'};
    Cols = [88, 101, 218; 178, 22, 22]./255;
    figure('Color',[1 1 1],'Name','Delta Metab energy, neglecting energy dissipation [nondim]');
    mk = 7;
    nc = 6;
    nr = 2;
    for ih = 1:2
        Cs = Cols(ih,:);

        subplot(2,6,4:6);
        plot([-5 8], [-5 8],'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        xHeader = 'Exp_PNetMetab_BM';
        yHeader = yHeaders{ih};
        % ------ Browning ------
        iSel = strcmp(SimInfo(:,1),'Browning');
        iRef = strcmp(SimInfo(:,1),'Browning') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0;
        DeltaBrowning = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaBrowning(:,strcmp(Headers,xHeader));
        ydat = DeltaBrowning(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

        % ---- Schertzer -----
        iSel = strcmp(SimInfo(:,1),'Schertzer');
        iRef = find(strcmp(SimInfo(:,1),'Schertzer') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) >= 1.1 &...
            Table(:,strcmp(Headers,'speed')) <= 1.2);
        DeltaDat = Table(iSel,:) - Table(iRef,:);

        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(2) = l(1); legN{2} = SimInfo{iRef};

        % ---- Gomenuka -----
        iSel = strcmp(SimInfo(:,1),'Gomenuka');
        iRef = find(strcmp(SimInfo(:,1),'Gomenuka') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.1 &...
            Table(:,strcmp(Headers,'speed')) < 1.2) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);

        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(3) = l(1); legN{3} = SimInfo{iRef};


        % ---- Huang -----
        iSel = strcmp(SimInfo(:,1),'Huang');
        iRef = find(strcmp(SimInfo(:,1),'Huang') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.3);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(4) = l(1); legN{4} = SimInfo{iRef};

        %  ------  Koelewijn ------
        iSel = strcmp(SimInfo(:,1),'Koelewijn');
        iRef = find(strcmp(SimInfo(:,1),'Koelewijn') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.2 &...
            Table(:,strcmp(Headers,'speed')) < 1.4);
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(5) = l(1); legN{5} = SimInfo{iRef};

        % %  ------  Umberger ------
        % iSel = strcmp(SimInfo(:,1),'Umberger');
        % iRef = find(strcmp(SimInfo(:,1),'Umberger') & ...
        %     Table(:,strcmp(Headers,'Slope')) == 0 & ...
        %     Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) > 0.87 & ...
        %     Table(:,strcmp(Headers,'Exp_StrideFreq')) < 0.9);
        % DeltaDat = Table(iSel,:) - Table(iRef,:);
        % l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
        %     DeltaDat(:,strcmp(Headers,yHeader)),...
        %     'o','Color',Cs,'MarkerSize',mk);

        %  ------  Jordan ------ (only spatio-temporal data)

        % ---- Abe -----
        iSel = strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'speed')) < Settings.vAbe2015_lim;
        iRef = find(strcmp(SimInfo(:,1),'Abe2015') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,'speed')) > 1.05 &...
            Table(:,strcmp(Headers,'speed')) < 1.15) ;
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(6) = l(1); legN{6} = SimInfo{iRef};

        % ---- McDonald -----
        iSel = strcmp(SimInfo(:,1),'McDonald');
        iRef = find(strcmp(SimInfo(:,1),'McDonald') & ...
            Table(:,strcmp(Headers,'Slope')) == 0 & ...
            Table(:,strcmp(Headers,'AddedMass')) == 0 & ...
            Table(:,strcmp(Headers,xHeader)) < 6); % the other one is the crouch walking condition
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*Mass_Sim,ydat*Mass_Sim,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
        leg(7) = l(1); legN{7} = SimInfo{iRef};

        % legend and layout
        if ih ==2
            xlabel('\Delta Pnet Exp. (W/kg)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W/kg)','Interpreter','tex')
            set(gca,'box','off')
            set(gca,'FontSize',14);
            set(gca,'LineWidth',1.5);
        end

        % part 2 figure 1Plot change in metabolic energy for each type of intervention

        %-------------------------------------------
        %           effect of walking speed
        %-------------------------------------------
        %   Schertzer
        %   Koelewijn
        %   Abe 2015
        %   Gomenuka
        %
        subplot(nr,nc,7:8);
        if ~Settings.DeltaValues
            plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-6 18]*Mass_Sim, [-6 18]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','speed','MaxSpeedSel',Settings.vAbe2015_lim);
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','speed');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','speed');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('speed')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);



        %------------------------------
        % ---- Loaded Walking ---------
        %------------------------------
        subplot(nr,nc,9:10);
        if ~Settings.DeltaValues
            plot([0 10]*Mass_Sim, [0 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-1 3]*Mass_Sim, [-1 3]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Browning','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Browning'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Browning'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Browning'; ctLeg = ctLeg + 1;

        % --- Schertzer -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Schertzer','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Schertzer'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Schertzer'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim, Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Schertzer'; ctLeg = ctLeg + 1;

        % --- Huang -----
        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Huang','addedmass');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Huang'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Huang'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp, scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Huang'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('loaded')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

        %-------------------------------
        %----- Walking on a slope ------
        %-------------------------------
        subplot(nr,nc,11:12);
        if ~Settings.DeltaValues
            plot([0 20]*Mass_Sim, [0 20]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        else
            plot([-5 10]*Mass_Sim, [-5 10]*Mass_Sim,'--','Color',[0 0 0],'LineWidth',1.3); hold on;
        end
        leg = [];    legN = []; ctLeg = 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Koelewijn','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Koelewijn'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Koelewijn'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Gomenuka','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Gomenuka'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Gomenuka'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'McDonald','slope');
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'McDonald'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'McDonald'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'McDonald'; ctLeg = ctLeg + 1;

        IndexPlot = SelectIndices_PlotTable(Table,SimInfo,Headers,'Abe2015','slope','MaxSpeedSel',Settings.vAbe2015_lim);
        Lexp = Exp_SubjProp.LegLength(strcmp(Exp_SubjProp.Study,'Abe2015'));
        mexp = Exp_SubjProp.mass(strcmp(Exp_SubjProp.Study,'Abe2015'));
        scale_exp =  mexp*9.81^(3/2)*sqrt(Lexp);
        scale_sim = Mass_Sim*9.81^(3/2)*sqrt(Lexp);
        l = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,...
            scale_exp,scale_sim,Settings.DeltaValues,Mass_Sim);
        leg(ctLeg) = l(1); legN{ctLeg} = 'Abe2015'; ctLeg = ctLeg + 1;

        % legend(leg,legN);
        title('slope')
        if Settings.DeltaValues
            xlabel('\Delta Pnet Exp. (W)','Interpreter','tex')
            ylabel('\Delta Pnet Sim. (W)','Interpreter','tex')
        else
            xlabel('Pnet Exp. (W)','Interpreter','tex')
            ylabel('Pnet Sim. (W)','Interpreter','tex')
        end
        set(gca,'box','off')
        set(gca,'FontSize',14);
        set(gca,'LineWidth',1.5);

    end




end

