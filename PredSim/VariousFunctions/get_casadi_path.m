function [casadi_path] = get_casadi_path()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

try casadi.GlobalOptions.getCasadiPath();
    casadi_path = casadi.GlobalOptions.getCasadiPath();
    disp(['casadi installation found in ', casadi_path]);
catch
    disp('no casadi installation found. please download and install casadi here');
    disp('https://web.casadi.org/');
end
end