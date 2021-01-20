%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% Function: SOLVE_DYSON 
% Solve Dyson's equation to find normal (G) and anomalous (F)
% Green's function from normal (S) and anomalous (P) self-energy.
%
% >>Usage:
% [G F] = solve_dyson(S,P,WN,ek,mu,Nk,Nc,useSymmetry,varargin)
%
% >>Input:
% S:    normal self-energy in (k,i\omega_n) space.
% P:    anomalous self-energy in (k,i\omega_n) space.
% WN:   fermionic Matsubara frequencies.
% ek:   dispersion on (2*Nk(1))*(2*Nk(2)) k-grid
% mu:   chemical potential
% norb: number of orbitals (has to be 2)
% Nk:   Nk = [Nk(1) Nk(2)], k-grid and r-grid.
% Nc:   2*Nc = number of Matsubara frequencies.
% useSymmetry: use C_2 symmetry for k-grid. Nk(1) and Nk(2) must both be
%       even number if useSymmetry=1. Set useSymmetry=0 and Nk = [1/2 1/2] 
%       to calculate one k-point (\Gamma point).
% varargin: If set to 0, then only calculate the normal Greens function G
%
% >>Output:
% G:    normal Green's function in (k,i\omega_n) space.
% F:       anomalous Green's function in (k,i\omega_n) space.t_filling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [G F] = solve_dyson(S,P,WN,ek,mu,norb,Nk,Nc,useSymmetry,varargin)


% useFreqSymmetry: =1 
% Use frequency symmetry G(k,-i*w_n) = conj[G(k,i*w_n)]
% and F(k,-i*w_n) = conj[F(k,i*w_n)].
useFreqSymmetry = 1;


% Only 2-orbital systems supported
if norb~=2
    error('Currently, only two-orbital case is supported!')
end

% Calculate G and F or only G?
if isempty(varargin)
    solveCase = 1;  %default: calculate both G and F
else
    solveCase = 0;  %calculate only G
end

% Check the sizes of the input arrays
siz_S = size(S);
siz_e = size(ek);
if ~isequal(siz_S(1:3),[2*Nk(1) 2*Nk(2) 2*Nc])
    error(' Input array S has a wrong size.')
end
if ~isequal(siz_e(1:2),[2*Nk(1) 2*Nk(2)])
    error(' Input array ek has a wrong size.')
end
if size(WN,1)~=1
    WN = reshape(WN,1,[]);
end
nel = siz_S(4);     %number of elements of the matrix in orbital space stored
if ~isequal(siz_e(3),nel)
    error(' Input array ek or S has a wrong size.')
end


% If using momentum space symmetry, only take one quarter of the BZ
if useSymmetry
    S = S(1:Nk(1)+1,1:Nk(2)+1,:,:);
    P = P(1:Nk(1)+1,1:Nk(2)+1,:,:);
    ek = ek(1:Nk(1)+1,1:Nk(2)+1,:);
    nk1 = Nk(1)+1;
    nk2 = Nk(2)+1;
else
    nk1 = 2*Nk(1);
    nk2 = 2*Nk(2);
end

% If using frequency symmetry, only take one half of the frequency points
if useFreqSymmetry
    S = S(:,:,1:Nc,:);
    P = P(:,:,1:Nc,:);
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

% Select the right orbitals and reshape S, P
S = S(:,:,:,idxmap);
P = P(:,:,:,idxmap);
S = reshape(S,npt,norbsq);   %now S becomes [nk1*nk2*nw norbsq]
P = reshape(P,npt,norbsq);

% P is Hermitian, so use that symmetry
P(:,3) = conj(P(:,3));
P(:,[1 4]) = real(P(:,[1 4]));

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

% Add self-energy, invert, multiply by the anomalous Greens function
% See equation (58) of notes01_FLEX.pdf
Mtmp = mul2(inv2(S-iG0),P);
% Eq. 58a
G = inv2(iG0 - S - mul2(P,conj(Mtmp)));
% If we want to find F as well
if solveCase == 1
    % Eq. 58c
    F = mul2(Mtmp,conj(G));
    % Note that F is Hermitian
    F(:,[1 4]) = real(F(:,[1 4]));   
    % Reshape F
    F = reshape(F(:,[1 2]),[nk1, nk2, nwpt, nel]);
else
    F = [];
end
% Reshape G
G = reshape(G(:,[1 2]),[nk1, nk2, nwpt, nel]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change back to (2*Nk(1), 2*Nk(2), 2*Nc, nel)-array, if we used
% frequency symmetry
if useFreqSymmetry
    G = G(:,:,[1:end, end:-1:1],:);
    G(:,:,(Nc+1):end,:) = conj(G(:,:,(Nc+1):end,:));
    if solveCase == 1
        F = F(:,:,[1:end, end:-1:1],:);
        F(:,:,(Nc+1):end,:) = conj(F(:,:,(Nc+1):end,:));
    end
end
% and momentum symmetry
if useSymmetry
    G = G([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    if solveCase == 1
        F = F([1:end, end-1:-1:2],[1:end, end-1:-1:2],:,:);
    end
end
