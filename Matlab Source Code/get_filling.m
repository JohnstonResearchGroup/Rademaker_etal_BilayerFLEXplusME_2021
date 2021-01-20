%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 26 March 2020
%
% >>Function: get_filling
% Compute the filling given the electron dispersion and self-energies.
% The description is given in Appendix C of the notes.
%
% >>Description:
% We substract the high-frequency tail using the parameter xifill, so that
% the remainder of the sum over Matsubara frequencies converges faster.
%
% >>Input:
% S:        normal self-energy (4-dimensional array: kx, ky, wn, orbital)
% P:        anomalous self-energy (4-dimensional array: kx, ky, wn, orbital)
% WN:       Matsubara Frequency grid
% ek:       Electron dispersion (3-dimensional array: kx, ky, orbital)
% xifill:   Determines the high-frequency tail of the bare Greens function,
%           typically T*log(2*Norb/filling0 - 1)
% norb:     Nr. of orbitals (parameter passed on to solve_dyson)
% beta:     Inverse temperature
% useSymmetry: Whether to use symmetry tricks to speed up the calculation
%              (parameter passed on to solve_dyson)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = get_filling(S,P,WN,ek,mu,xifill,norb,beta,useSymmetry)

% S and P should be 4-dimensional array and ek 3-dimensional array. 
% G only includes all the diagonal matrix elements in orbital space. 

% Get the size of array S
[nk1 nk2 nw nel] = size(S);
% Check whether the frequency grid size is compatible
if mod(nw,2) ~= 0 || nw ~= numel(WN)
    error('Wrong number of frequency grid!')
end
% Convert grid size numbers
Nk(2) = nk2/2;
Nk(1) = nk1/2;
nkpt = nk1*nk2;
numwi = nw/2;
% Compute the Greens function using solve_dyson
G = solve_dyson(S,P,WN,ek,mu,norb,Nk,numwi,useSymmetry,0);
% Extract only the [1 1] orbital (WHY?)
G = G(:,:,:,[1 1]);     % orb1 is equivalent to orb2

% See Eqn. (C3) from notes, to be updated
fill = 2*norb*fermi(xifill,beta);       %for two spins

% Reshape the Greens function
G = reshape(permute(G,[1 2 4 3]),nkpt*norb,[]);
% Take only the real part of the Matsubara Greens function, because we are
% summing G(-iwn) + G(+iwn)
G = real(G(:,1:numwi));

% Sum over Matsubara frequencies
for nn = 1:numwi
    wn = WN(nn);
    % Add to the filling 
    fill = fill + 4*sum(G(:,nn))/(nkpt*beta) ...
        + (4*norb/beta)*xifill/(wn^2 + xifill^2);
end

% Output
out = fill;
