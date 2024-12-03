function [MechE] = getMechPower(R, model_info)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% mechanical work at joints

tau =  R.kinetics.T_ID;
qd = R.kinematics.Qdots_rad;

% power for each joint
Power = tau .* qd;

% total mechanical power
TotalPower = sum(Power,2);

t = R.time.mesh_GC(1:end-1)';
JointWork = trapz(t,Power);
TotalWork = trapz(t,TotalPower);

% muscle power
vM = R.muscles.vM;
Fce = R.muscles.Fce;
MusclePower = Fce.*-vM;

% % we miss the energy dissipated in muscle damper here (or is this included
% in Fce)
Fmo = [model_info.muscle_info.parameters().FMo];
vMtilde = R.muscles.vMtilde;
F_damping = vMtilde*R.S.misc.dampingCoefficient.*Fmo;
P_muscle_damping = vM.*F_damping;
A_muscle_damping =  trapz(t, P_muscle_damping);
A_muscle_dampingtot = sum(A_muscle_damping);

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

% store outcomes
MechE.JointPower = Power;
MechE.TotalPower = TotalPower;
MechE.JointWork = JointWork;
MechE.TotalWork = TotalWork;
MechE.MusclePower = MusclePower;
MechE.MuscleWork = MuscleWork;
MechE.PosPower = PosPower;
MechE.NegPower = NegPower;
MechE.MusclePosWork = MusclePosWork;
MechE.MuscleNegWork = MuscleNegWork;
MechE.MusclePosWorkTotal = MusclePosWorkTotal;
MechE.MuscleNegWorkTotal = MuscleNegWorkTotal;
MechE.A_muscle_dampingtot = A_muscle_dampingtot;


% also compute metabolic energy margaria model
MechE.P_Marg = PosPower./0.25 - NegPower./1.2;
MechE.COT_Marg = sum(trapz(t,MechE.P_Marg),2)./R.spatiotemp.dist_trav./model_info.mass;

% compute net power
dt = t(end)-t(1);
MechE.Pnet_Joints = TotalWork./dt;
MechE.Pnet_PosMuscles = MusclePosWorkTotal./dt;
MechE.Pnet_NegMuscles = MuscleNegWorkTotal./dt;

% energy torque actuqtors
iTorqueAct = [model_info.actuator_info.parameters().coordi];
P_torque = R.torque_actuators.T .*  R.kinematics.Qdots_rad(:,iTorqueAct);
A_torque = trapz(t,P_torque);
MechE.A_Torque = A_torque;
MechE.A_Torque_tot = sum(A_torque);

% joint damping
P_damping = R.kinetics.T_damping.*R.kinematics.Qdots_rad;
A_damping = trapz(t, P_damping);
A_damping_tot = sum(A_damping);
MechE.A_damping_tot =  A_damping_tot;
% check energy equation
% periodic motion so delta Ekin is = 0 
% work done by gravity is also = 0
% assumption that -TotalWork (inverse dynamics) equals energy dissipated in
% ground contact

A_contact = - TotalWork;
MechE.error_energyeq = MechE.MusclePosWorkTotal + MechE.MuscleNegWorkTotal +...
    MechE.A_Torque_tot + A_contact + MechE.A_damping_tot;

% also add metabolic power
Edot = R.metabolics.Bhargava2004.Edot_gait;
metab_work = trapz(t,Edot);
MechE.Barhava_metab = sum(metab_work);
MechE.Bhargava_efficiency = MusclePosWork./metab_work;




end