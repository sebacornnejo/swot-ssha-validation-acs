function [L, var, EOFs, PC, Mr] = EOF_(M,n)
% Function to calculate all eigenvectors (EOFs) and eigenvalues
% (lambda) of a data matrix using singular value decomposition (svd).
%
% Usage: 
% [L, var, EOFs, PC, Mr] = EOF_(M,n)
%
% Inputs:
% M : Data matrix of size n x m (time x space). 
% n : Number of EOFs used to reconstruct the data matrix. If not
%     specified, only the first EOF is used.
%
% Outputs:
% L   : Vector containing the eigenvalues of M.
% var : Vector containing the percentage of variance explained by each eigenvalue.
% EOFs: Matrix with eigenvectors (EOFs) in each column. 
% PC  : Matrix with principal components (PCs) in each column.
% Mr  : Reconstructed matrix using the first n EOFs. n=1 by default. 
%
% Darinka Pecarevic E.
% LOOPSO/PUCV
% 2017-08-29
%
% Modified by:
% Mauro Pinto Juica
% 2019-01-22

[r, c] = size(M); % Size of matrix M
if nargin < 2
    n = 1;
end
    
[C, gamma, EOFs] = svd(M, 'econ'); % SVD calculation
L   = diag(gamma).^2 / (r-1);      % Eigenvalues
var = L ./ sum(L) * 100;           % Percentage of explained variance
PC  = M * EOFs;                    % Principal components calculation 
% PC = C * gamma;                  % Can also be calculated this way 
Mr  = PC(:, 1:n) * EOFs(:, 1:n)';  % Reconstructed matrix 
end