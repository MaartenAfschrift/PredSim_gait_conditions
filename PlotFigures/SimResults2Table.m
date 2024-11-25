function [output] = SimResults2Table(SimFolder,varargin)
%SimResults2Table Loads the simulation results and selects a subset of the
%outcomes
import org.opensim.modeling.*;

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

output = nan(1,27);
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
        Etot_nobasal(i,:) = energy_total;
    end
    E_sum_GC = trapz(Res.R.time.mesh_GC(1:end-1),Etot);
    output(12) = E_sum_GC/Res.R.misc.body_mass/Res.R.spatiotemp.dist_trav;
    output(13) = output(12).*output(4); % metabolic power
    output(14) = output(13).*output(6)./RefModelMass; % W/kg subject mass (and not W/kg total mass)

    % compute metabolic power with slightly adapted heath rate coefficient
    b = 10000;
    scaleRate = 2;
    N = length(exc(:,1));
    FMo = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'FMo');
    Etot = nan(N,1);
    for i=1:N
        [energy_total,Adot,Mdot,Sdot,Wdot,energy_model] = ...
            getMetabolicEnergySmooth2004all_scaleheatrate(exc(i,:)',act(i,:)',lMtilde(i,:)',vM(i,:)',...
            Fce(i,:)',Fpass(i,:)',musclemass,pctst,Fiso(i,:)',FMo,modelmass,b,scaleRate);
        Etot(i) = full(energy_model);
    end
    E_sum_GC = trapz(Res.R.time.mesh_GC(1:end-1),Etot);
    output(15) = E_sum_GC/Res.R.misc.body_mass/Res.R.spatiotemp.dist_trav;
    output(16) = output(15).*output(4); % metabolic power
    output(17) = output(16).*output(6)./RefModelMass; % W/kg subject mass (and not W/kg total mass)

    %% mechanical energy
    % compute positive and negative work done by muscle fibers (and net work done)
    tau =  Res.R.kinetics.T_ID;
    qd = Res.R.kinematics.Qdots_rad;

    % power for each joint
    Power = tau .* qd;

    % total mechanical power
    TotalPower = sum(Power,2);

    % mechanical work at joints
    t = Res.R.time.mesh_GC(1:end-1)';
    JointWork = trapz(t,Power);
    TotalWork = trapz(t,TotalPower);

    % muscle power
    vM = -Res.R.muscles.vM;
    Fce = Res.R.muscles.Fce;
    MusclePower = Fce.*vM;

    % Muscle work
    MuscleWorkInd = trapz(t,MusclePower);
    MuscleWork    = sum(MuscleWorkInd);

    % positive and negative muscle power
    PosPower = MusclePower;
    PosPower(PosPower<0) = 0;
    NegPower = MusclePower;
    NegPower(NegPower>0) = 0;
    MusclePosWork = trapz(t,PosPower);
    MuscleNegWork = trapz(t,NegPower);
    MusclePosWorkTotal = sum(MusclePosWork);
    MuscleNegWorkTotal = sum(MuscleNegWork);
    output(18) = TotalWork;
    output(19) = MusclePosWorkTotal;
    output(20) = MuscleNegWorkTotal;

    %% other metablic energy related parameters
    output(21) = mean(mean(Res.R.muscles.a));
    output(22) = mean(mean(Res.R.muscles.a.^2));

    %% metabolic work

    output(23) = TotalWork./(t(end)-t(1));
    output(24) = MusclePosWorkTotal./(t(end)-t(1));
    output(25) = MuscleNegWorkTotal./(t(end)-t(1));

    % metabolic work
    MetabWork = sum(trapz(t,Etot_nobasal));
    output(26) = MetabWork./(t(end)-t(1));

    %% margaria based metabolic work based on the idea that the energy
    % dissipation in joints and contact model is wrong.

    % we have to do some tricks here because I don't measure energy
    % dissipated in ground contact model
    %   - limit cycle so Delta Ekin = 0
    Delta_Ekin = 0;
    %   - compute work against gravity (a bit hard because we tilted) the
    %   gravity vector
    %       read osim model
    % osim_model = Res.model_info.osim_path;
    % if exist(osim_model,'file')
    %     my_model = Model(osim_model);
    %
    %     g = my_model.getGravity();
    %     g = [g.get(0), g.get(1), g.get(2)];
    fi = asin(slope/100);
    R = rotz(fi);
    gv = [0, -9.81, 0];
    g = gv*R(1:3, 1:3);
    G = Res.model_info.mass * g;
    r = Res.R.kinematics.Qs(:,4:6); % pelvis position
    dr = r(end,:)-r(1,:);
    A_gravity = dr*G'; % conservative force so we can do force times displacement

    % energy equation: A_ground + A_gravity + A_intern = Delta Ekin
    A_ground = Delta_Ekin - TotalWork - A_gravity;
    if A_ground>0
        disp('error');
    end

    % joint damping
    P_damping = Res.R.kinetics.T_damping.*Res.R.kinematics.Qdots_rad;
    A_damping = trapz(t, P_damping);
    A_jointdamping = sum(A_damping);
    if A_jointdamping >0
        disp('error');
    end

    % positive work muscles without energy losses in model
    Energy_dissipated = A_jointdamping + A_ground;
    A_muscles_pos = MusclePosWorkTotal + Energy_dissipated;
    output(27) = (A_muscles_pos/0.25 + MuscleNegWorkTotal/1.2)./RefModelMass ./(t(end)-t(1));
    %
    % else
    %     disp(['cannot find ', osim_model]);
    %     output(27) = NaN;
    % end


end

end