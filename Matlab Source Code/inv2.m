%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% >>Function: inv2
% Invert a twoxtwo matrix that is saved as a 1d array with four elements.
% This is mainly used by solve_dysonbwd.m.
%
% >> Input:
% A:       Input matrix, as a 1d array of length 4.
%
% >> Output:
% A:       Inverted 2x2 matrix, same structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = inv2(A)

% This is how the matrix elements are stored
a11 = A(:,1); a21 = A(:,2); a12 = A(:,3); a22 = A(:,4);

% Inverse Determinant
invd = 1./(a11.*a22 - a21.*a12);

% Output the inverse
A(:,1) = a22.*invd;
A(:,2) = -a21.*invd;
A(:,3) = -a12.*invd;
A(:,4) = a11.*invd;
