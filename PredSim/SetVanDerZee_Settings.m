function [S] = SetVanDerZee_Settings(S, stride_freq, speed, stepwidth)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% set stride frequency
if ~isempty(stride_freq) && ~isnan(stride_freq)
    dtStep = 1./stride_freq;
    S.bounds.t_final.lower = dtStep-0.01;
    S.bounds.t_final.upper = dtStep+0.01;
end

% set walking speed
S.subject.v_pelvis_x_trgt = speed;

% constrain stepwidth (now ugly implementation)
if ~isempty(stepwidth) && ~isnan(stepwidth)
    S.bounds.calcn_dist.lower = stepwidth;
end

end