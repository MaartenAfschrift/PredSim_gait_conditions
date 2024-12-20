function [f_forceEquilibrium_FtildeState_all_tendon,f_FiberLength_TendonForce_tendon,...
    f_FiberVelocity_TendonForce_tendon,f_lT_vT] = createCasadi_ContractDynam(S,model_info)
% --------------------------------------------------------------------------
% createCasadi_ContractDynam
%   Function to create Casadi functions for muscle contraction dynamics.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - f_forceEquilibrium_FtildeState_all_tendon -
%   * function for force equilibrium between muscle fiber and tendon
%
%   - f_FiberLength_TendonForce_tendon -
%   * function to compute fiber length of a muscle
%
%   - f_FiberVelocity_TendonForce_tendon -
%   * function to compute fiber lenghtening velocity of a muscle
%
%   - f_lT_vT -
%   * function to compute tendon length and lengthening velocity of a tendon
% 
% Original author: Ines Vandekerckhove, Tom Buurke & Dhruv Gupta, KU Leuven
% Original date: 30-11-2021 
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

import casadi.*
N_muscles = model_info.muscle_info.NMuscle;



%% Muscle contraction dynamics
% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde',N_muscles); % Normalized tendon forces
a           = SX.sym('a',N_muscles); % Muscle activations
dFTtilde    = SX.sym('dFTtilde',N_muscles); % Time derivative tendon forces
lMT         = SX.sym('lMT',N_muscles); % Muscle-tendon lengths
vMT         = SX.sym('vMT',N_muscles); % Muscle-tendon velocities
tension_SX  = SX.sym('tension',N_muscles); % Tensions
% atendon_SX  = SX.sym('atendon',NMuscle); % Tendon stiffness
% shift_SX    = SX.sym('shift',NMuscle); % shift curve
Hilldiff    = SX(N_muscles,1); % Hill-equilibrium
FT          = SX(N_muscles,1); % Tendon forces
Fce         = SX(N_muscles,1); % Contractile element forces
Fiso        = SX(N_muscles,1); % Normalized forces from force-length curve
vMmax       = SX(N_muscles,1); % Maximum contraction velocities
massM       = SX(N_muscles,1); % Muscle mass
Fpass       = SX(N_muscles,1); % Passive element forces
% Parameters of force-length-velocity curves

% load(fullfile(MainPath,'MuscleModel','Fvparam.mat'),'Fvparam');
% load(fullfile(MainPath,'MuscleModel','Fpparam.mat'),'Fpparam');
% load(fullfile(MainPath,'MuscleModel','Faparam.mat'),'Faparam');
load('Fvparam.mat','Fvparam');
load('Fpparam.mat','Fpparam');
load('Faparam.mat','Faparam');

% Parameters of force-length-velocity curves
for m = 1:N_muscles
    [Hilldiff(m),FT(m),Fce(m),Fpass(m),Fiso(m),vMmax(m),massM(m)] = ...
        ForceEquilibrium_FtildeState_all_tendon(a(m),FTtilde(m),...
        dFTtilde(m),lMT(m),vMT(m),model_info.muscle_info.parameters(m).FMo,...
        model_info.muscle_info.parameters(m).lMo,model_info.muscle_info.parameters(m).lTs,...
        model_info.muscle_info.parameters(m).alphao,model_info.muscle_info.parameters(m).vMmax,...
        Fvparam,Fpparam,Faparam,tension_SX(m),model_info.muscle_info.parameters(m).tendon_stiff,...
        model_info.muscle_info.parameters(m).tendon_stiff_shift,S.misc.constant_pennation_angle,...
        S.misc.dampingCoefficient,model_info.muscle_info.parameters(m).muscle_pass_stiff_shift,...
        model_info.muscle_info.parameters(m).muscle_pass_stiff_scale,...
        model_info.muscle_info.parameters(m).muscle_strength,S);
end
f_forceEquilibrium_FtildeState_all_tendon = ...
    Function('f_forceEquilibrium_FtildeState_all_tendon',{a,FTtilde,...
    dFTtilde,lMT,vMT,tension_SX},{Hilldiff,FT,Fce,Fpass,Fiso,vMmax,massM},...
    {'a','FTtilde','dFTtilde','lMT','vMT','tension_SX'},...
    {'Hilldiff','FT','Fce','Fpass','Fiso','vMmax','massM'});

% Function to get (normalized) muscle fiber lengths
lM      = SX(N_muscles,1);
lMtilde = SX(N_muscles,1);
lT      = SX(N_muscles,1);
for m = 1:N_muscles
    [lM(m),lMtilde(m),lT(m)] = FiberLength_TendonForce_tendon(FTtilde(m),...
        model_info.muscle_info.parameters(m).lMo,model_info.muscle_info.parameters(m).lTs,...
        model_info.muscle_info.parameters(m).alphao,lMT(m),...
        model_info.muscle_info.parameters(m).tendon_stiff,...
        model_info.muscle_info.parameters(m).tendon_stiff_shift,...
        S.misc.constant_pennation_angle);
end
f_FiberLength_TendonForce_tendon = Function(...
    'f_FiberLength_Ftilde_tendon',{FTtilde,lMT},{lM,lMtilde},...
    {'FTtilde','lMT'},{'lM','lMtilde'});

% Function to get (normalized) muscle fiber velocities
vM      = SX(N_muscles,1);
vMtilde = SX(N_muscles,1);
vT      = SX(N_muscles,1);
for m = 1:N_muscles
    [vM(m),vMtilde(m),vT(m)] = FiberVelocity_TendonForce_tendon(FTtilde(m),...
        dFTtilde(m),model_info.muscle_info.parameters(m).lMo,model_info.muscle_info.parameters(m).lTs,...
        model_info.muscle_info.parameters(m).alphao,model_info.muscle_info.parameters(m).vMmax,lMT(m),...
        vMT(m),model_info.muscle_info.parameters(m).tendon_stiff,...
        model_info.muscle_info.parameters(m).tendon_stiff_shift,...
        S.misc.constant_pennation_angle);
end
f_FiberVelocity_TendonForce_tendon = Function(...
    'f_FiberVelocity_Ftilde_tendon',{FTtilde,dFTtilde,lMT,vMT},...
    {vM,vMtilde},{'FTtilde','dFTtilde','lMT','vMT'},{'vM','vMtilde'});

f_lT_vT = Function('f_lT_vT',{FTtilde,dFTtilde,lMT,vMT},...
    {lT,vT},{'FTtilde','dFTtilde','lMT','vMT'},{'lT','vT'});
end