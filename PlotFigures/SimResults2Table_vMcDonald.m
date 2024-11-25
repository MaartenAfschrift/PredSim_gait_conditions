function [output] = SimResults2Table_vMcDonald(SimFolder,varargin)
%SimResults2Table Loads the simulation results and selects a subset of the
%outcomes

% defaults slope
slope = 0;
RefModelMass = 62; % model mass without added mass
AddedMass = 0;
if ~isempty(varargin)
    nInputs = floor(length(varargin)/2);
    for i=1:nInputs
        if strcmp(varargin{i*2-1},'slope')
            slope = varargin{i*2};
        elseif strcmp(varargin{i*2-1},'RefModelMass')
            RefModelMass = varargin{i*2};
        elseif strcmp(varargin{i*2-1},'AddedMass')
            AddedMass = varargin{i*2};
        end
    end
end

% Load the results file (assumes that this the only .mat file in the folder)
matFiles = dir(fullfile(SimFolder,'*.mat'));
if ~isempty(matFiles)
    if length(matFiles) > 1
        disp([ num2str(length(matFiles))  ' ,mat files in in folder ', SimFolder, '  selected ' matFiles(1).name]);
    end
    Res = load(fullfile(matFiles(1).folder,matFiles(1).name));
else
    disp(['No .mat files in ', SimFolder]);
    Res = [];
end

output = nan(1,14+7);
if ~isempty(Res)
    % extract all outputs
    output(1) = Res.R.metabolics.Bhargava2004.COT;
    output(2) = Res.R.spatiotemp.stride_freq;
    output(3) = Res.R.spatiotemp.step_width_COP;
    output(4) = Res.R.S.subject.v_pelvis_x_trgt;
    output(5) = output(1).*output(4); % metabolic power
    output(6) = Res.R.misc.body_mass;
    output(7) = slope;
    output(8) = output(5).*output(6)./RefModelMass; % W/kg subject mass (and not W/kg total mass)
    output(9) = AddedMass;
    if ~isfield(Res.R.metabolics,'Marg')
        % compute metabolic energy Marg
        Res.R.metabolics.Marg.Power = getMetabolicEnergy_MargariaSmooth(Res.R.muscles.Fce,Res.R.muscles.vM,1000);
        NetPower = sum(Res.R.metabolics.Marg.Power,2);
        Work = trapz(Res.R.time.mesh_GC(1:end-1),NetPower);
        Res.R.metabolics.Marg.COT = Work/Res.R.misc.body_mass/Res.R.spatiotemp.dist_trav;
    end
    output(10) = Res.R.metabolics.Marg.COT;
    output(11) = output(10).*output(4).*output(6)./RefModelMass; % W/kg subject mass (and not W/kg total mass), marg model
    % compute metabolic energy consumption again with other tanh smoothing
    % factor
    exc = Res.R.muscles.e;
    act = Res.R.muscles.a;
    lMtilde = Res.R.muscles.lMtilde;
    vM = Res.R.muscles.vM;
    Fce = Res.R.muscles.Fce;
    Fpass = Res.R.muscles.Fpass;
    pctst = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');
    musclemass = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'muscle_mass');
    Fiso = Res.R.muscles.Fiso;
    modelmass = Res.model_info.mass;
    %     b = Res.R.S.metabolicE.tanh_b;
    b = 10000;
    N = length(exc(:,1));
    FMo = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'FMo');
    Etot = nan(N,1);
    for i=1:N
        [energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
            getMetabolicEnergySmooth2004all(exc(i,:)',act(i,:)',lMtilde(i,:)',vM(i,:)',...
            Fce(i,:)',Fpass(i,:)',musclemass,pctst,Fiso(i,:)',FMo,modelmass,b);
        Etot(i) = full(energy_model);
    end
    E_sum_GC = trapz(Res.R.time.mesh_GC(1:end-1),Etot);
    output(12) = E_sum_GC/Res.R.misc.body_mass/Res.R.spatiotemp.dist_trav;
    output(13) = output(12).*output(4); % metabolic power
    output(14) = output(13).*output(6)./RefModelMass; % W/kg subject mass (and not W/kg total mass)
    [At,Atsq] = GetA_McDonald(Res);
    output(15:21) = At;
    output(22:28) = Atsq;

    % add activation cost dependend on muscle mass.
    MassAct = act.*musclemass';
    output(29) = nanmean(sum(MassAct,2));
    output(30) = nanmean(sum(MassAct.^2,2));
    


end

end