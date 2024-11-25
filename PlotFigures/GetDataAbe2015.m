function [Dat,header] = GetDataAbe2015()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Dat = nan(12*3,3);

speeds =[2.4 3.1 3.8 4.5 5.2 5.9 6.6 7.3 8.7 9.4 10.1 10.8]./3.6; % conver t to m/s
grade = [0, -5, 5];
VO2_COT_level = [239.5 200.3 174.1 164.1 163.2 172.7 185.1 203.3 211.1 212.5 209 205.7];
VO2_COT_downhill = [201.2 162.2 142.2 130.6 126.1 131.1 138.6 154.2 171.4 170.6 166.9 167.4];
VO2_COT_uphill = [317.8 270.9 247.4 236.2 236.2 240.3 251.7 267.6 276.8 274.9 272.6 270.2];
VO2_COT = [VO2_COT_level'; VO2_COT_downhill'; VO2_COT_uphill']./1000; % in mL per kg per m

for i=1:3
    iSel = (i-1)*12+1:i*12;
    Dat(iSel,1) = speeds;
    Dat(iSel,2) = grade(i);    
end
Dat(:,3) = VO2_COT;
VO2_rate = VO2_COT.*Dat(:,1); % mL/kg/s
C02_rate = 0.85*VO2_rate; % this is a very simple assumption as VCO2 is not reported
MetabE_rate = 16.89.*VO2_rate + 4.84.*C02_rate; % in J/kg/s
MetabE_COT = MetabE_rate./Dat(:,1);% in J/kg/m
Dat(:,4) = MetabE_rate;
Dat(:,5) = MetabE_COT;
header = {'speed','grade','VO2','MetabPower','Metab_COT'};


% figure(); plot(Dat(:,1),MetabE_COT,'ok')





end