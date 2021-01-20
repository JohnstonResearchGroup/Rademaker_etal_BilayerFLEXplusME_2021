%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 20 January 2021
%
% Function: nnormgauss
% 
% Returns Non-normalized Gaussian with standard deviation sigma. Input variable x can be an
% array as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = nnormgauss(eps_k,sigma)
    out = exp(-eps_k.*eps_k/(2*sigma*sigma));