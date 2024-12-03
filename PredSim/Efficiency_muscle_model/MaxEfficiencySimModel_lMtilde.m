%% get mechanical efficiency for each muscle in our model

SetFigureDefaults();
Cs = [0.4 0.4 0.4];

% simulation resutls
addpath(genpath('C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\NeuromechanicsToolkit'));
addpath(genpath('C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\PredSim'));
MainDPath = 'C:\Users\mat950\Documents\Software\Publications\PredSim_gait_conditions\SimResults';

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
lMtilde_v = linspace(0.5,1.5,200);
b = 10000; % smoothing factor
E_all = nan(length(musclemass), length(vMtilde), length(lMtilde_v));
E_max = nan(length(musclemass), length(lMtilde_v));
for j = 1:length(lMtilde_v)
    lMtilde = lMtilde_v(j);
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
        E_all(i,:,j) = Pmech./Bharg.energy_total;
        E_max(i,j) = max(E_all(i,:));
    end
end


figure()
subplot(2,2,1:2)
[X, Y] = meshgrid(lMtilde_v, vMtilde);
iM = strcmp(Res.model_info.muscle_info.muscle_names,'soleus_r');
Z = squeeze(E_all(iM,:,:));
surf(X, Y, Z);  % Surface plot

% Add labels and a title
xlabel('lMtilde');
ylabel('vMtilde');
zlabel('efficiency');
title('efficiency ');

% Add a color bar
colorbar;

subplot(2,2,3)
iSel = find(vMtilde <-2,1,'first');
plot(lMtilde_v,Z(iSel,:));
xlabel('lMtilde');
ylabel('efficiency');
title('efficiency at lM_tilde_dot = -2');

subplot(2,2,4)
iSel = find(lMtilde_v > 1,1,'first');
plot(vMtilde,Z(:,iSel));
xlabel('vMtilde');
ylabel('efficiency');
title('efficiency at lMtilde = 0');




% figure(); 
% subplot(1,4,2)
% 
% histogram(E_max,40,'FaceColor', Cs, 'EdgeColor',Cs);
% set(gca,'box','off');
% xlabel('mechanical efficiency []')
% ylabel('number of muscles');
% 
% subplot(1,4,1)
% plot(vMtilde, E_all,'k','Color',Cs);
% xlabel('rel. contraction velocity');
% ylabel('mechanical efficiency []');
% set(gca,'box','off');
% 
% [~, imax] = max(E_max);
% [~, imin] = min(E_max);
% disp(Res.model_info.muscle_info.muscle_names{imax})
% disp(Res.model_info.muscle_info.muscle_names{imin})
