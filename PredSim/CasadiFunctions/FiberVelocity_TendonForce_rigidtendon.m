function [vM,vMtilde,varargout] = FiberVelocity_TendonForce_rigidtendon(FTtile,...
    lMo_in,lTs_in,alphao_in,vMmax_in,lMT,vMT,MuscMoAsmp)
% --------------------------------------------------------------------------
% FiberVelocity_TendonForce_tendon
%    This function computes muscle fiber velocities from muscle-tendon forces.
%    More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%   
% INPUT:
%
% OUTPUT:
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   Adapted to allow assumption of constant pennation angle, by Lars D'Hondt.
%   Adapted to return tendon lengthening velocity, by Lars D'Hondt
% Last edit by: Lars D'Hondt
% Last edit date: 30/May/2022
% --------------------------------------------------------------------------

lMo = ones(size(FTtile,1),1)*lMo_in;
lTs = ones(size(FTtile,1),1)*lTs_in;
alphao = ones(size(FTtile,1),1)*alphao_in;
vMmax = ones(size(FTtile,1),1)*vMmax_in;

% Inverse tendon force-length characteristic
lTtilde = 1;

% Hill-type muscle model: geometric relationships
vT = 0;
if(MuscMoAsmp == 0) % b = cst
    lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
    cos_alpha = (lMT-lTs.*lTtilde)./lM;
else    % alpha = cst = alphao
    cos_alpha = cos(alphao);
end
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;

if nargout == 3
    varargout{1} = vT;
end

end