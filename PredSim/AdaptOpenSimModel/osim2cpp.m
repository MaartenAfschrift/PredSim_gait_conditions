function [] = osim2cpp(Cpp2Dll,CppDir,osim_path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

% disp
disp('Convert .osim to .dll file started ...');
% list with all input arguments to create the .dll file
pathOpenSimModel = osim_path;
outputFilename = ['F_' osim_file_name];
export3DSegmentOrigins = Cpp2Dll.export3DSegmentOrigins;
jointsOrder = Cpp2Dll.jointsOrder;
coordinatesOrder = Cpp2Dll.coordinatesOrder;
exportGRFs = Cpp2Dll.exportGRFs;
exportSeparateGRFs = Cpp2Dll.exportSeparateGRFs;
exportGRMs = Cpp2Dll.exportGRMs;
exportContactPowers = Cpp2Dll.exportContactPowers;

% Create the .cpp file using an executable
if isempty(Cpp2Dll.PathCpp2Dll_Exe)
    error(['You should specify the path to the executables to convert', ...
        'osim file to an .dll file. You can install the exectubale using the' ...
        'function: InstallOsim2Dll_Exe']);
else
    PathCpp2Dll_Exe = Cpp2Dll.PathCpp2Dll_Exe;
end

% create folders if needed
if ~isfolder(CppDir)
    mkdir(CppDir);
end

% create opensim bin folder if needed (quick fix for issues with use of
% python api in exectutable. (To Do: update this)


% create a mat file with all input information
MatFileExoInfo = fullfile(pwd,'MatFileExoInfo.mat');
% delete old file
if isfile(MatFileExoInfo)
    delete(MatFileExoInfo);
end
pathOut = CppDir;
save(MatFileExoInfo,'pathOpenSimModel','outputFilename','export3DSegmentOrigins',...
    'jointsOrder','coordinatesOrder','exportGRFs','exportSeparateGRFs',...
    'exportGRMs','exportContactPowers','pathOut');

% run the exectuable
MainPath = pwd;
cd(PathCpp2Dll_Exe);
command = ['osimtocppexe.exe "' MatFileExoInfo '"'];
disp('   converting .osim file to .cpp file to solve ID');
if Cpp2Dll.verbose_mode == 0
    [~,~] = system(command);
else
    system(command);
end
cd(MainPath);
delete(MatFileExoInfo);

end