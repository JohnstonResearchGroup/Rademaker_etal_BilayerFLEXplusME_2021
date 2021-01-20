%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% >>Function: interpSelfE
% Interpolate the self-energy from one Matsubara frequency grid to another
% (meaning, at a different temperature).
%
% >>Usage:
% Two possible modes: if the first variable argument is 1, then we use the
% existing selfenergies in the arrays S, P, WN and mu (input). If the first
% variable argument is 0, then we import the self-energy from a file.
%
% >> Input:
% WNin:     New frequency grid for the self-energy, as defined by
%           WN = pi*(2*(-numwi:numwi-1)+1)/beta.
% varargin{1}: 1 = from input self-energy, 0 = from file;
% varargin{2}: normal self-energy (1) or filename (0)
% varargin{3}: anomalous self-energy
% varargin{4}: Matsubara frequency grid of self-energy
% varargin{5}: chemical potential
%
% >> Output:
% Sut:      Interpolated normal self-energy
% Put:      Interpolated anomalous self-energy
% mu:       Chemical potential (unchanged by this function)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Sut Put mu] = interpSelfE(WNin,varargin)

% From file or from matrices?
flySelfE = varargin{1};
if flySelfE == 1
    S = varargin{2};
    P = varargin{3};
    WN = varargin{4};
    mu = varargin{5};
else
    dataFileName = varargin{2};
    load(dataFileName,'S','P','WN','mu','-mat')
end

% Check the array sizes and input Matsubara grid
[nk1 nk2 nw nel] = size(S);        %nel>=1
nwin = numel(WNin);
if mod(nw,2)~=0 || mod(nwin,2)~=0
    error('Number of Matsubara frequencies must be even!')
end
if nw~=numel(WN)
    error('Wrong grid for self-energy data!')
end

% Reshape the self-energy arrays
S = permute(S,[3 1 2 4]);          %now S becomes [nw nk1 nk2 nel]
P = permute(P,[3 1 2 4]);
% Clear the output self-energy arrays.
Sut(1:nwin,1:nk1,1:nk2,1:nel) = 0;
Put(1:nwin,1:nk1,1:nk2,1:nel) = 0;

% 1D interpolation with spline-extrapolated values
Sut = interp1(WN,S,WNin,'spline','extrap');
Put = interp1(WN,P,WNin,'spline','extrap');

% Replace spline-extrapolated values (for new frequencies outside the 'old'
% grid) by 'nearest' values
% First find boundaries
idxL = find(WNin<WN(1));
idxR = find(WNin>WN(end));
% If there are frequencies 'sticking' out on the lower side
if ~isempty(idxL)
    idxL = idxL(end);
    Sut(1:idxL,:,:,:) = Sut(repmat(idxL+1,1,idxL),:,:,:);
    Put(1:idxL,:,:,:) = Put(repmat(idxL+1,1,idxL),:,:,:);
end
% If there are frequencies 'sticking' out on the higher side
if ~isempty(idxR)
    idxR = idxR(1);
    Sut(idxR:nwin,:,:,:) = Sut(repmat(idxR-1,1,nwin-idxR+1),:,:,:);
    Put(idxR:nwin,:,:,:) = Put(repmat(idxR-1,1,nwin-idxR+1),:,:,:);
end

% Reshape the self-energy arrays in original form
Sut = permute(Sut,[2 3 1 4]);      %now Sut becomes [nk1 nk2 nwin nel]
Put = permute(Put,[2 3 1 4]);