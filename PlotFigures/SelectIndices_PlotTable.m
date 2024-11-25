function [Index] = SelectIndices_PlotTable(Table,SimInfo,Headers,ExpName,Effect,varargin)
%SelectIndices_PlotTable Selects the correct index in the data table
% (ToDo: proper description how this works)

% defaults slope

% default speed bounds depending on the type of experiment
SpeedBoundRef = [1.1 1.2];
if strcmp(ExpName,'Browning')
    SpeedBoundRef = [0 99]; % constant speed experiment
elseif strcmp(ExpName,'Schertzer')
    SpeedBoundRef = [1.1 1.2];
elseif strcmp(ExpName,'Gomenuka')
    SpeedBoundRef = [1.1 1.2];
elseif strcmp(ExpName,'Huang')
    SpeedBoundRef = [1.2 1.3];
elseif strcmp(ExpName,'Koelewijn')
    SpeedBoundRef = [1.2 1.4];
elseif strcmp(ExpName,'Umberger')
    SpeedBoundRef = [0 99]; % constant speed exp
elseif strcmp(ExpName,'Abe2015')
    SpeedBoundRef = [1.05 1.15];
elseif strcmp(ExpName,'McDonald')
    SpeedBoundRef = [0 99];
elseif strcmp(ExpName,'Jordan')
    SpeedBoundRef = [0.95 1.05];
elseif strcmp(ExpName,'VanDerZee2022')
    SpeedBoundRef = [1.05 1.15];
elseif strcmp(ExpName,'Strutzenberger')
    SpeedBoundRef = [1.05 1.15];
end

% adapt bounds on reference speed based on user input
maxSpeedSel = 999;
if ~isempty(varargin)
    nInputs = floor(length(varargin)/2);
    for i=1:nInputs
        if strcmp(varargin{i*2-1},'SpeedBoundRef')
            SpeedBoundRef = varargin{i*2};
            %         elseif strcmp(varargin{i*2-1},'RefModelMass')
            %             RefModelMass = varargin{i*2};
            %         elseif strcmp(varargin{i*2-1},'AddedMass')
            %             AddedMass = varargin{i*2};
        elseif strcmp(varargin{i*2-1}, 'MaxSpeedSel')
            maxSpeedSel = varargin{i*2};
        end
    end
end

% DeltaSpeed
% The simulation is not always exactly at the same speed as in the
% experimental (especially in an early phase of the research project). To
% solve this issue I search for the reference condition between a minimum
% and maximum bound of the target speed
DeltaSpeed = 0.05;


% all possible slopes, masslocations, added mass and walking speeds in the
% experiment
iExp = strcmp(SimInfo(:,1),ExpName);
Slopes = unique(Table(iExp,strcmp(Headers,'Slope')));
MassLoc = unique(SimInfo(iExp,2));
mAdded = unique(Table(iExp,strcmp(Headers,'AddedMass')));
speeds = unique(Table(iExp,strcmp(Headers,'speed')));
Index = [];


if strcmp(Effect,'speed')
    % select indicies for various walking speeds in each unique condition
    ct = 1;
    for i = 1:length(Slopes)
        for j = 1:length(MassLoc)
            for k = 1:length(mAdded)
                iSel = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'Slope')) == Slopes(i) & ...
                    Table(:,strcmp(Headers,'AddedMass')) == mAdded(k) & ...
                    strcmp(SimInfo(:,2), MassLoc{j}) & ...
                    Table(:,strcmp(Headers,'speed')) < maxSpeedSel);

                iRef = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'Slope')) == Slopes(i) & ...
                    Table(:,strcmp(Headers,'AddedMass')) == mAdded(k) & ...
                    strcmp(SimInfo(:,2), MassLoc{j}) & ...
                    Table(:,strcmp(Headers,'speed')) > SpeedBoundRef(1) &...
                    Table(:,strcmp(Headers,'speed')) < SpeedBoundRef(2)) ;
                if ~isempty(iSel) && ~isempty(iRef)
                    Index(ct).iSel = iSel;
                    Index(ct).iRef = iRef;
                    Index(ct).dat = {Slopes(i),mAdded(k),MassLoc{j},NaN,ExpName};
                    Index(ct).dat_header = {'slope','addedmass','masslocation','speed','ExpName'};
                    ct = ct+1;
                end
            end
        end
    end
elseif strcmp(Effect,'slope')
    % select indicies for various walking speeds in each unique condition
    ct = 1;
    for i = 1:length(speeds)
        for j = 1:length(MassLoc)
            for k = 1:length(mAdded)
                iSel = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'speed')) == speeds(i) & ...
                    Table(:,strcmp(Headers,'AddedMass')) == mAdded(k) & ...
                    strcmp(SimInfo(:,2), MassLoc{j}) & ...
                    Table(:,strcmp(Headers,'speed')) < maxSpeedSel);
                iRef = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'Slope')) == 0 & ...
                    Table(:,strcmp(Headers,'AddedMass')) == mAdded(k) & ...
                    strcmp(SimInfo(:,2), MassLoc{j}) & ...
                    Table(:,strcmp(Headers,'speed')) >= speeds(i)-DeltaSpeed &...
                    Table(:,strcmp(Headers,'speed')) <= speeds(i)+DeltaSpeed);
                if ~isempty(iSel) && ~isempty(iRef)
                    Index(ct).iSel = iSel;
                    Index(ct).iRef = iRef;
                    Index(ct).dat = {NaN,mAdded(k),MassLoc{j},speeds(i),ExpName};
                    Index(ct).dat_header = {'slope','addedmass','masslocation','speed','ExpName'};
                    ct = ct+1;
                end
            end
        end
    end
elseif strcmp(Effect,'addedmass')
    % select indicies for various walking speeds in each unique condition
    ct = 1;
    for i = 1:length(Slopes)
        for j = 1:length(MassLoc)
            for k = 1:length(speeds)
                iSel = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'Slope')) == Slopes(i) & ...
                    Table(:,strcmp(Headers,'speed')) == speeds(k) & ...
                    strcmp(SimInfo(:,2), MassLoc{j}) & ...
                    Table(:,strcmp(Headers,'speed')) < maxSpeedSel);
                iRef = find(strcmp(SimInfo(:,1),ExpName) & ...
                    Table(:,strcmp(Headers,'Slope')) == Slopes(i) & ...
                    Table(:,strcmp(Headers,'speed')) >= speeds(k)-DeltaSpeed &...
                    Table(:,strcmp(Headers,'speed')) <= speeds(k)+DeltaSpeed &...
                    Table(:,strcmp(Headers,'AddedMass')) == 0);
                if ~isempty(iSel) && ~isempty(iRef)
                    Index(ct).iSel = iSel;
                    Index(ct).iRef = iRef;
                    Index(ct).dat = {Slopes(i),NaN,MassLoc{j},speeds(k),ExpName};
                    Index(ct).dat_header = {'slope','addedmass','masslocation','speed','ExpName'};
                    ct = ct+1;
                end
            end
        end
    end

end



end