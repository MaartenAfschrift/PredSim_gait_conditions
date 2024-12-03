%% get mechanical efficiency for each muscle in our model

SetFigureDefaults();
Cs = [0.4 0.4 0.4];

% simulation resutls
mainpath = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions';
addpath(genpath(fullfile(mainpath, 'NeuromechanicsToolkit')));
addpath(genpath(fullfile(mainpath, 'PredSim')));
MainDPath = fullfile(mainpath, 'SimResults');

% load simulation results to get model properties (I know a bit weird to
% do but this is more convenient for me)
Res = load(fullfile(MainDPath,'Fal22Ref_13','Falisse_et_al_2022_job25.mat'));

% get some model information
pctst = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'slow_twitch_fiber_ratio');
musclemass = struct_array_to_double_array(Res.model_info.muscle_info.parameters,'muscle_mass');
modelmass = Res.model_info.mass;
% maximal mechanical efficiency when full activated at optimal length
vmax = 10;
vMtilde = linspace(0,-10,200); %[lopt/s]
e = 1;
a = e;
lMtilde = 1;
b = 10000; % smoothing factor
E_all = nan(length(musclemass), length(vMtilde));
E_max = nan(length(musclemass), 1);
for i =1:length(musclemass)
    % model properties
    lceopt = Res.model_info.muscle_info.parameters(i).lMo;
    Fiso =  Res.model_info.muscle_info.parameters(i).FMo;
    vM = vMtilde.*lceopt;
    
    % force-length-velocity properties
    [Fpe_norm,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,vmax);
    Fce = a.*FMltilde*FMvtilde.*Fiso;
    Fpe = Fpe_norm.*Fiso;
    
    % metabolic power
    [Bharg.energy_total,Bharg.Adot,Bharg.Mdot,Bharg.Sdot,Bharg.Wdot,Bharg.energy_model] = ...
        getMetabolicEnergySmooth2004all(e,a,lMtilde,vM,Fce,Fpe,musclemass(i),pctst(i),FMltilde,...
        Fiso,modelmass,b);

    % mechanical work contractile element
    Pmech = -Fce.*vM;

    % max efficiency
    E_all(i,:) = Pmech./Bharg.energy_total;
    E_max(i) = max(E_all(i,:));   
end

figure(); 
subplot(1,4,2)

histogram(E_max,40,'FaceColor', Cs, 'EdgeColor',Cs);
set(gca,'box','off');
xlabel('mechanical efficiency []')
ylabel('number of muscles');

subplot(1,4,1)
plot(vMtilde, E_all,'k','Color',Cs);
xlabel('rel. contraction velocity');
ylabel('mechanical efficiency []');
set(gca,'box','off');

[~, imax] = max(E_max);
[~, imin] = min(E_max);
disp(Res.model_info.muscle_info.muscle_names{imax})
disp(Res.model_info.muscle_info.muscle_names{imin})
%% mechanical efficiency negative contractions

% maximal mechanical efficiency when full activated at optimal length
vmax = 10;
vMtilde = linspace(0,10,200); %[lopt/s]
e = 1;
a = e;
lMtilde = 1.3;
b = 10000; % smoothing factor
E_all = nan(length(musclemass), length(vMtilde));
E_max = nan(length(musclemass), 1);
for i =1:length(musclemass)
    % model properties
    lceopt = Res.model_info.muscle_info.parameters(i).lMo;
    Fiso =  Res.model_info.muscle_info.parameters(i).FMo;
    vM = vMtilde.*lceopt;
    
    % force-length-velocity properties
    [Fpe_norm,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,vmax);
    Fce = a.*FMltilde*FMvtilde.*Fiso;
    Fpe = Fpe_norm.*Fiso;
    
    % metabolic power
    [Bharg.energy_total,Bharg.Adot,Bharg.Mdot,Bharg.Sdot,Bharg.Wdot,Bharg.energy_model] = ...
        getMetabolicEnergySmooth2004all(e,a,lMtilde,vM,Fce,Fpe,musclemass(i),pctst(i),FMltilde,...
        Fiso,modelmass,b);

    % mechanical work contractile element
    Pmech = -Fce.*vM;

    % max efficiency
    E_all(i,:) = Pmech./Bharg.energy_total;
    E_max(i) = min(E_all(i,:));   
end

% figure(); 
subplot(1,4,4)


histogram(E_max,40,'FaceColor', Cs, 'EdgeColor',Cs);
set(gca,'box','off');
xlabel('mechanical efficiency []')
ylabel('number of muscles');

subplot(1,4,3)
plot(vMtilde, E_all,'Color',Cs);
xlabel('rel. contraction velocity');
ylabel('mechanical efficiency []');
set(gca,'box','off');
