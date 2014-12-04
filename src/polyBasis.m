%% Investigating the polynomial basis
% This script has the objective of analysing which terms should be use on
% the polynomial expression of the stream function $\psi$ when we use a 9
% points stencil.

%% Setting up the script
clear, clc, close all

N = 5;       % polynomial degree
numPts = 3;  % number of points in one dimension
%% Constructing the stream function
% Here we construct the polynomial stream function as
%
% $$
% \psi(x,y) = \sum_{i=0}^N \sum_{j=0}^N x^i y^j.
% $$
% 
% However, for the code below (using |N = 5|), we are removing the terms of
% 10th, 9th, 8th and 7th degrees (combining the degrees in $x$ and $y$).
% Moreover, we take out the following 6th degree terms: $x^5y$, $x^4y^2$
% and $x^3y^3$.
% 
% With this combination of terms, we get *full rank (18)* system with *22 
% unknows*. Alternatively, we could remove the following 6th degree terms
% to get a full rank system with 22 unknows:
% 
% * $xy^5$, $x^2y^4$ and $x^3y^3$;
% * $xy^5$, $x^3y^3$ and $x^4y^2$;
% * $x^2y^4$, $x^3y^3$ and $x^5y$;
% 
% It is quite possible to reduce this system to a *full rank* $18\times18$,
% i.e., using only *18 terms*. One example is to remove the following
% terms: $x^5$, $x^4$, $y^4$ and $y^3$.
%
% A important obeservation is that we could choose different terms just
% like we did for the 6th degree terms.

syms x y h
S = kron(x.^(0:N),y.^(0:N));
% Removing redundant terms ------------------------------------------------
S(1) = [];                        % Taking out the constant coefficient
if N == 5
    S(end) = [];                  % Taking out degree 10
    S(end - [5 0]) = [];          % Taking out degree 9
    S(end - [9 4 0]) = [];        % Taking out degree 8
    S(end - [12 7 3 0]) = [];     % Taking out degree 7
    %S(end - [14 9 5 2 0]) = [];  % Taking out degree 6 <== RANK DEFICIENT
    S(end - [5 2 0]) = [];        % Taking out some 6th degree terms
    S(end) = [];                  % Taking a 5th degree term
    S(end - [17 1]) = [];         % Taking out some 4th degree terms
    S(end - 16) = [];             % Taking out a 3rd degree term
end
% -------------------------------------------------------------------------

coeffs = sym('c', [length(S), 1]);
% S = S*coeffs;

Sx = diff(S, x); v = -Sx;
Sy = diff(S, y); u = +Sy;

uSymbolic = matlabFunction(u, 'vars', {[x y], coeffs},'file','uSymbolic');
vSymbolic = matlabFunction(v, 'vars', {[x y], coeffs},'file','vSymbolic');

[gridX, gridY] = meshgrid([-h 0 h],[h 0 -h]);

dSites = [gridX(:) gridY(:)];

U = sym(zeros(numPts^2,length(coeffs)));
for i = 1:numPts^2
    U(i,:) = uSymbolic(dSites(i,:),coeffs);
end
V = sym(zeros(numPts^2,length(coeffs)));
for i = 1:numPts^2
    V(i,:) = vSymbolic(dSites(i,:),coeffs);
end

M = [U; V];
rank(M)
