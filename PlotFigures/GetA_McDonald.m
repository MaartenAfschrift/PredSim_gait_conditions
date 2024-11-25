function [At,Atsq] = GetA_McDonald(Res)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% get the activation cost
a = Res.R.muscles.a;
muscles = Res.R.colheaders.muscles;
mList = {'glut_max1_r','glut_max2_r','glut_max3_r','glut_max1_l','glut_max2_l','glut_max3_l',...
    'bifemlh_r','bifemsh_r','bifemlh_l','bifemsh_l',...
    'rect_fem_r','rect_fem_l',...
    'vas_med_r','vas_med_l',...
    'med_gas_r','med_gas_l',...
    'soleus_r','soleus_l',...
    'tib_ant_l','tib_ant_r'};
igroup = [1 1 1 1 1 1, ...
    2 2 2 2,...
    3 3,...
    4 4,...
    5,5,...
    6,6,...
    7,7];
iMuscles = nan(length(mList),1);
for i=1:length(mList)
    iSel = find(strcmp(mList{i},muscles)); % find index muscle
    if ~isempty(iSel)
        iMuscles(i) = iSel;
    else
        disp(['Muscle ' mList{i} ' not in the simulation model']);
    end
end

% duration stride
dt = Res.R.time.mesh(end) - Res.R.time.mesh(1);
% compute activation of muscle groups
Groups = unique(igroup);
nGroups = length(Groups);
At = nan(nGroups,1);
for i=1:nGroups
    % all muscles in this group
    iSel = igroup == Groups(i);
    asel = a(:,iMuscles(iSel));
    AvAct = nanmean(asel); % average avation of each muscle in the gait cycle
    At(i) = nanmean(AvAct).*dt; % numerical integral of average actvation of muscles in group
    % integrate square
    Atsq(i) = nanmean(AvAct.^2).*dt; % numerical integral of average actvation of muscles in group
end

% integrate square


end