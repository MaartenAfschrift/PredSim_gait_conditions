function [Res] = LoadSimFile(SimFolder)
%LoadSimFile Load mat file with sim results
%   Detailed explanation goes here

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



end