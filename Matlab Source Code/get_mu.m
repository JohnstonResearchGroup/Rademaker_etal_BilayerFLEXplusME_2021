%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 26 March 2020
%
% >>Function: get_mu
% Find the chemical potential with bisection method.
% The bisection method is guaranteed to converge if DOS(w) >= 0. The code
% uses the matlab function 'fzero'.
%
%
% >>Input:
% S:        normal self-energy (4-dimensional array: kx, ky, wn, orbital)
% P:        anomalous self-energy (4-dimensional array: kx, ky, wn, orbital)
% WN:       Matsubara Frequency grid
% ek:       Electron dispersion (3-dimensional array: kx, ky, orbital)
% muLim:    Bounds of the chemical potential: (muL, muR)
% xifill:   Determines the high-frequency tail of the bare Greens function,
%           typically T*log(2*Norb/filling0 - 1)
% norb:     Nr. of orbitals (parameter passed on to solve_dyson)
% beta:     Inverse temperature
% useSymmetry: Whether to use symmetry tricks to speed up the calculation
%              (parameter passed on to solve_dyson)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function out = get_mu(S,P,WN,ek,muLim,xifill,norb,beta,useSymmetry)

% Set the threshold for tolerance in the Matlab Fzero method
erreps = 1e-8;

% Initial filling determined by xifill
fill = 2*norb*fermi(xifill,beta);       %for two spins

%set the initial values of the chemical potential.
muL = muLim(1);
muR = muLim(2);
muM = (muR+muL)/2;

% Get the fillings associaed with the left- and right-bound of the chemical
% potenial
fillL = get_filling(S,P,WN,ek,muL,xifill,norb,beta,useSymmetry);
fillR = get_filling(S,P,WN,ek,muR,xifill,norb,beta,useSymmetry);

% In some rare cases, the self-energy X is so large that fill is not in
% [fillL fillR], so we expand the interval [muL muR] below.
iter = 0;
Wband = muR - muL;
while ( (fill-fillL)*(fill-fillR)>0 )
    iter = iter + 1;
    if (fill > fillR)
        muR = muR + Wband*0.5;
        fillR = get_filling(S,P,WN,ek,muR,xifill,norb,beta,useSymmetry);
    elseif (fill < fillL)
        muL = muL - Wband*0.5;
        fillL = get_filling(S,P,WN,ek,muL,xifill,norb,beta,useSymmetry);
    end
    if (iter > 4)
        error('Self-energy too large, cannot find the chemical potential!')
    end
end
muLim = [muL muR];

% Find mu using Matlab fzero method (super-linear convergence) ------------
% fzero uses the algorithm that was originated by T. Dekker, which uses a
% combination of bisection, secant, and inverse quadratic interpolation
% methods. See Press et al., Numerical Recipes, Chapter 9.3 Van
% Wijngaarden-Dekker-Brent Method.
options = optimset;
options.TolX = erreps;
out = fzero(@(xmu) get_filling(S,P,WN,ek,xmu,xifill,norb,beta,useSymmetry)-fill,muLim,options);
