%% get mechanical efficiency for each muscle in our model

clear all; close all; clc;
SetFigureDefaults();
Cs = [0.4 0.4 0.4];

% simulation resutls
addpath(genpath('C:\Users\mat950\Documents\Software\DataAnalysis\NeuromechanicsToolkit'));
addpath(genpath('C:\Users\mat950\Documents\Software\Sim\PredInt\src\PredSim'));
MainDPath = 'C:\Users\mat950\OneDrive - Vrije Universiteit Amsterdam\Onderzoek\SimResults';
ResPath= fullfile(MainDPath,'PredSimResults');

Res = load(fullfile(ResPath,'Fal22Ref_10_v2','Falisse_et_al_2022_job296.mat'));

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

    % Umberger 2003
    for j = 1:length(vMtilde)
        [E_Umb2003(j),energy_am(j),energy_sl(j),energy_mech(j),E_Umb2003_model(j)] = ...
            getMetabolicEnergySmooth2003all(e,a,lMtilde,vMtilde(j),vM(j),Fce(j),...
            musclemass(i),pctst(i),10,FMltilde,modelmass,b);
    end

    % Umberger 2010
    for j = 1:length(vMtilde)
        [E_Umb2010(j),~,~,~,E_Umb2010_model(j)] = ...
            getMetabolicEnergySmooth2010all(e,a,lMtilde,vMtilde(j),vM(j),Fce(j),...
            musclemass(i),pctst(i),10,FMltilde,modelmass,b);
    end

    % Uchida 2016
    for j = 1:length(vMtilde)
        [E_Uch2016(j),~,~,~,E_Uch2016_model(j)] = ...
            getMetabolicEnergySmooth2016all(e,a,lMtilde,vMtilde(j),vM(j),Fce(j),...
            musclemass(i),pctst(i),10,FMltilde,modelmass,b);
    end


    % mechanical work contractile element
    Pmech = -Fce.*vM;

    % max efficiency
    E_all(i,:) = Pmech./Bharg.energy_total;
    E_max(i) = max(E_all(i,:));   

    E_all_Umb2003(i,:) = Pmech./E_Umb2003;
    E_max_Umb2003(i) = max(E_all_Umb2003(i,:));   

    E_all_Umb2010(i,:) = Pmech./E_Umb2010;
    E_max_Umb2010(i) = max(E_all_Umb2010(i,:));   

    E_all_Uch2016(i,:) = Pmech./E_Uch2016;
    E_max_Uch2016(i) = max(E_all_Uch2016(i,:)); 
    E_mech(i,:) = energy_mech;
end

figure(); 

Cols = linspecer(5);

subplot(1,2,1)
l1 = plot(vMtilde, E_all,'Color',[0.4 0.4 0.4]); hold on;
l2 = plot(vMtilde, E_all_Umb2003,'Color',Cols(4,:)); hold on;
% plot(vMtilde, E_all_Umb2010,'Color',Cols(3,:)); hold on;
% plot(vMtilde, E_all_Uch2016,'Color',Cols(4,:)); hold on;
xlabel('rel. contraction velocity');
ylabel('mechanical efficiency []');
set(gca,'box','off');
legend([l1(1), l2(2)],{'Bhargava','Umberger'});


subplot(1,2,2)
histogram(E_max,40,'FaceColor', [0.4 0.4 0.4], 'EdgeColor',[0.4 0.4 0.4]); hold on;
histogram(E_max_Umb2003, 40,'FaceColor', Cols(4,:), 'EdgeColor',Cols(4,:));
set(gca,'box','off');
xlabel('mechanical efficiency []')
ylabel('number of muscles');






% 
% %% mechanical efficiency negative contractions
% 
% % maximal mechanical efficiency when full activated at optimal length
% vmax = 10;
% vMtilde = linspace(0,10,200); %[lopt/s]
% e = 1;
% a = e;
% lMtilde = 1.3;
% b = 10000; % smoothing factor
% E_all = nan(length(musclemass), length(vMtilde));
% E_max = nan(length(musclemass), 1);
% for i =1:length(musclemass)
%     % model properties
%     lceopt = Res.model_info.muscle_info.parameters(i).lMo;
%     Fiso =  Res.model_info.muscle_info.parameters(i).FMo;
%     vM = vMtilde.*lceopt;
% 
%     % force-length-velocity properties
%     [Fpe_norm,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,vmax);
%     Fce = a.*FMltilde*FMvtilde.*Fiso;
%     Fpe = Fpe_norm.*Fiso;
% 
%     % metabolic power
%     [Bharg.energy_total,Bharg.Adot,Bharg.Mdot,Bharg.Sdot,Bharg.Wdot,Bharg.energy_model] = ...
%         getMetabolicEnergySmooth2004all(e,a,lMtilde,vM,Fce,Fpe,musclemass(i),pctst(i),FMltilde,...
%         Fiso,modelmass,b);
% 
%     % mechanical work contractile element
%     Pmech = -Fce.*vM;
% 
%     % max efficiency
%     E_all(i,:) = Pmech./Bharg.energy_total;
%     E_max(i) = min(E_all(i,:));   
% end
% 
% % figure(); 
% subplot(1,4,4)
% 
% 
% histogram(E_max,40,'FaceColor', Cs, 'EdgeColor',Cs);
% set(gca,'box','off');
% xlabel('mechanical efficiency []')
% ylabel('number of muscles');
% 
% subplot(1,4,3)
% plot(vMtilde, E_all,'Color',Cs);
% xlabel('rel. contraction velocity');
% ylabel('mechanical efficiency []');
% set(gca,'box','off');
% 
