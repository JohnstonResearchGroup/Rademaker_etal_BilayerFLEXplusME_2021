%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bilayer FLEX + Phonons
%      S. Johnston & Y. Wang & L. Rademaker & G. Alvarez-Suchini
%      Last update: 3 April 2020
%
% >>Function: mul2
% Product of two 2x2 matrices that are saved as 1d arrays with 4 elements.
% This is mainly used by solve_dysonbwd.m.
%
% >> Input:
% A:       Input matrix 1, as a 1d array of length 4.
% B:       Input matrix 2, as a 1d array of length 4.
%
% >> Output:
% A:       Product AxB, as 1d array with length 4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = mul2(A,B)

% This is how the matrix elements are stored
a11 = A(:,1); a21 = A(:,2); a12 = A(:,3); a22 = A(:,4);
b11 = B(:,1); b21 = B(:,2); b12 = B(:,3); b22 = B(:,4);

% Compute the product
A(:,1) = a11.*b11 + a12.*b21;
A(:,2) = a21.*b11 + a22.*b21;
A(:,3) = a11.*b12 + a12.*b22;
A(:,4) = a21.*b12 + a22.*b22;
