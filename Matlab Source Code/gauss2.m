%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 20 March 2020
%
% Function: gauss
% 
% Returns Gaussian with standard deviation s. Input variable x can be an
% array as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = gauss2(eps_k,sigma)
    out = exp(-eps_k*eps_k/(2*sigma*sigma))/sqrt(2*pi*sigma*sigma);
