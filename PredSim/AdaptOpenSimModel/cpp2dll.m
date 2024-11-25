function [] = cpp2dll(CppDir,outputFilename,dll_file,compiler,PathCpp2Dll_Exe,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% load .mat file with input and outputs
MainPath = pwd;
verbose_mode = 0;
if ~isempty(verbose_mode)
    verbose_mode = varargin{1};
end

outputFilename = ['F_' outputFilename];

% path to the opensim-AD submodule
pathFunction = mfilename('fullpath');
[pathPreProcessing, ~, ~] = fileparts(pathFunction);
[pathMain, ~] = fileparts(pathPreProcessing);
PathOsimAD = fullfile(pathMain,'opensimAD');

cd(PathCpp2Dll_Exe);
load(fullfile(CppDir,[outputFilename '_IO.mat']),'IO');

% Create foo file
cd(MainPath);
disp('   creating project for .cpp file using cmake');
[fooPath] = buildExternalFunction(outputFilename, CppDir, compiler,...
    PathOsimAD, verbose_mode);

% create foo_jac (temporary solution using executable);
disp('   converting foo.py to foo_jac.c using casadi');
cd(PathCpp2Dll_Exe)
command = ['GenF.exe ' fooPath ' ' num2str(IO.nCoordinates*3)]; % ToDo: Add nInputes to IO structure
if verbose_mode == 0
    [~,~] = system(command);
else
    system(command);
end
cd(MainPath);

% createDll
disp('   creating .dll file using cmake');
CreateDllFromFoo_Jac(outputFilename, CppDir, compiler, PathOsimAD,...
    verbose_mode);
cd(MainPath);

% display message
disp('   convert .osim to .dll file finished ...');
end