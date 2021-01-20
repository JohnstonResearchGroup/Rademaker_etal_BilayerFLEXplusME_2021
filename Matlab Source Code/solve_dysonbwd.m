%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% Function: solve_dysonbwd 
% Solve Dyson's equation backward to find normal (S) and anomalous (P) 
% self-energy from the normal (G) and anomalous (F) Green's function.
%
% Usage:
% [S P] = solve_dysonbwd(G,F,WN,ek,mu,norb,Nk,Nc,useSymmetry,varargin)
%
% >>Input:
% G:        normal self-energy in (k,i\omega_n) space.
% F:        anomalous self-energy in (k,i\omega_n) space.
% WN:       fermionic Matsubara frequencies.
% ek:       dispersion on (2*Nk(1))*(2*Nk(2)) k-grid
% mu:       chemical potential
% norb:     number of orbitals (has to be 2)
% Nk:       Nk = [Nk(1) Nk(2)], k-grid and r-grid.
% Nc:       2*Nc = number of Matsubara frequencies.
% useSymmetry: use C_2 symmetry for k-grid. Nk(1) and Nk(2) must both be
%           even number if useSymmetry=1. Set useSymmetry=0 and 
%           Nk = [1/2 1/2] to calculate one k-point (\Gamma point).
% varargin: if default (1), calculate S and P, otherwise only S.
%
% >>Output:
% S:        normal Green's function in (k,i\omega_n) space.
% P:        anomalous Green's function in (k,i\omega_n) space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [S P] = solve_dysonbwd(G,F,WN,ek,mu,norb,Nk,Nc,useSymmetry,varargin)


% useFreqSymmetry: =1 
% Use frequency symmetry G(k,-i*w_n) = conj[G(k,i*w_n)]
% and F(k,-i*w_n) = conj[F(k,i*w_n)].
useFreqSymmetry = 1;


% Only 2-orbital systems supported
if norb~=2
    error('Currently, only two-orbital case is supported!')
end


% Calculate S&P or only S?
if isempty(varargin)
    solveCase = 1;  %default: calculate both S and P
else
    solveCase = 0;  %calculate only S
end


% Check the sizes of the input arrays
siz_G = size(G);
siz_e = size(ek);
if ~isequal(siz_G(1:3),[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array G has a wrong size.')
end
if ~isequal(siz_e(1:2),[2*Nk(1) 2*Nk(2)])
    error(' Input array ek has a wrong size.')
end
if size(WN,1)~=1
    WN = reshape(WN,1,[]);
end
nel = siz_G(4);     %number of elements of the matrix in orbital space stored
if ~isequal(siz_e(3),nel)
    error(' Input array ek or G has a wrong size.')
end

% If using momentum space symmetry, only take one quarter of the BZ
if useSymmetry
    G = G(1:Nk(1)+1,1:Nk(2)+1,:,:);
    F = F(1:Nk(1)+1,1:Nk(2)+1,:,:);
    ek = ek(1:Nk(1)+1,1:Nk(2)+1,:);
    nk1 = Nk(1)+1;
    nk2 = Nk(2)+1;
else
    nk1 = 2*Nk(1);
    nk2 = 2*Nk(2);
end
% If using frequency symmetry, only take one half of the frequency points
if useFreqSymmetry
    G = G(:,:,1:Nc,:);
    F = F(:,:,1:Nc,:);
    WN = WN(1:Nc);
    nwpt = Nc;
else
    nwpt = 2*Nc;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here the actual Dyson's equation starts

% Set the array sizes
nkpt = nk1*nk2;
npt = nkpt*nwpt;
norbsq = norb*norb;

% how to order the orbital degrees of freedom
idxmap = [1 2 2 1];       %orb1 is equivalent to orb2
%idxmap = [1 2 2 3];       %orb1 is not equivalent to orb2

% Select the right orbitals and reshape G, F
G = G(:,:,:,idxmap);
F = F(:,:,:,idxmap);
G = reshape(G,npt,norbsq);   %now G becomes [nk1*nk2*nw norbsq]
F = reshape(F,npt,norbsq);
% F is Hermitian 
F(:,3) = conj(F(:,3));
F(:,[1 4]) = real(F(:,[1 4]));

% Compute the inverse of the noninteracting Greens function (iG0)
ek = ek(:,:,idxmap);
iG0(1:npt,1:norbsq) = 0;
for ss = 1:4
    if ss==1 || ss==4
        iG0(:,ss) = complex(reshape(repmat(mu - ek(:,:,ss), [1, 1, nwpt]),npt,1), ...
            reshape(repmat(WN, [nkpt, 1]),npt,1));
    else
        iG0(:,ss) = reshape(repmat(-ek(:,:,ss), [1, 1, nwpt]),npt,1);
    end
end

% Compute the S
Mtmp = mul2(inv2(G),F);
S = iG0 - inv2(G + mul2(F,conj(Mtmp)));
% Compute P if asked
if solveCase == 1
    % Compute P
    P = mul2(Mtmp,conj(S-iG0));
    % P is Hermitian
    P(:,[1 4]) = real(P(:,[1 4]));
    % Reshape P
    P = reshape(P(:,[1 2]),[nk1, nk2, nwpt, nel]);
else
    P = [];
end
% Reshape S
S = reshape(S(:,[1 2]),[nk1, nk2, nwpt, nel]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change back to (2*Nk(1), 2*Nk(2), 2*Nc, nel)-array, if we used
% frequency symmetry
if useFreqSymmetry
    S = S(:,:,[1:end, end:-1:1],:);
    S(:,:,(Nc+1):end,:) = conj(S(:,:,(Nc+1):end,:));
    if solveCase == 1
        P = P(:,:,[1:end, end:-1:1],:);
        P(:,:,(Nc+1):end,:) = conj(P(:,:,(Nc+1):end,:));
    end
end
% and momentum symmetry
if useSymmetry
    S = S([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    if solveCase == 1
        P = P([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    end
end
