%% Lebesgue functions of divergence-free RBFs interpolants
% This function calculates cardinal and Lebesgue function (global
% interolation). For more information on how the divergence-free RBF
% interpolant is created, see the file <RBF_DivFreeMatrix.html
% |RBF_DivFreeMatrix.m|> and the paper _"Divergence-Free RBFs on Surfaces"_
% of _Narcowich_, _Ward_ and _Wright_.

%%
function [lebesgueConstU, lebesgueConstV] = ...
                                  lebesgueFunctionsRBF(n, display, ep, rbf)
% Calculates Lebesgue function and its constant. The function uses an
% equispaced grid in x- and y-direction.
%
% n       : Number of points in x- and y-direction
% display : displays the cardinal function (optional, default = false)
% ep      : Shape parameter (optional, default = 2)
% rbf     : Radial basis function (optional,
%                                  default = @(ep,r) exp(-(ep*r).^2) )

%% Setting up the function
if nargin == 1
    display = false;
end
if nargin <= 2
    ep = 2;
end
if nargin <= 3
    rbf = @ (ep,r) exp(-(ep*r).^2);
end

%% Generating the 2d grid of n^2 equispaced nodes
x = linspace(-1,1,n);
xx = linspace(-1,1,7*n);

[X,Y] = meshgrid(x);
[XX,YY] = meshgrid(xx);

dSites = [X(:) Y(:)];   % interpolation points
ePoints = [XX(:) YY(:)]; % evaluation points

%% Calculating cardinal functions
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[A, Dx, Dy, F, G] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);

L = zeros(n,2*n);  % used to construct delta function

lebesgueFunctionU = zeros(size(XX));
lebesgueFunctionV = zeros(size(XX));

for p = 1:2*n^2  % all interpolation points
    L(p) = 1;
    U = L(:,1:n);
    V = L(:,n+1:2*n);
    t =  [U(:) V(:)];  % vector field
    L(p) = 0;          % going back to zero matrix
    
    d = reshape(t', 1, numel(t))';
    coeffs = A\d;
    coeffs = reshape(coeffs, 2, size(t,1))';  % I'm not keeping the coeffs
    
    cardFunction = RBFdivFreeInterp(coeffs, ePoints, dSites, F, G);
    cardFunctionU = reshape(cardFunction(:,1), size(XX));
    cardFunctionV = reshape(cardFunction(:,2), size(XX));

    [i,j] = ind2sub([n, n], rem(p-1,n^2) + 1); 
    if display
        figure(1)
        set(gcf, 'Position', [100,100, 600*2, 600])
        h1 = subplot(1,2,1);
        mesh(XX,YY,cardFunctionU)
        axis([-1 1 -1 1 -.5 1])
        title(['Cardinal Function U for (i,j) = (' num2str(i) ',' ...
               num2str(j) ')'])
    
        h2 = subplot(1,2,2);
        mesh(XX,YY,cardFunctionV)
        axis([-1 1 -1 1 -.5 1])
        title(['Cardinal Function V for (i,j) = (' num2str(i) ',' ...
               num2str(j) ')'])

        set([h1 h2], 'clim', [-.5 1])
        snapnow
    end

    lebesgueFunctionU = lebesgueFunctionU + abs(cardFunctionU);
    lebesgueFunctionV = lebesgueFunctionV + abs(cardFunctionV);
end

%% 
% Above we show the cardinal functions for the case where |n = 3|. Observe
% that the divergence-free interpolant of the zero function is not zero,
% however it interpolates zeros at the interpolation points. This happens
% because the interpolant is divergence-free.

%% Plotting the Lebesgue functions and the Lebesgue constant
lebesgueConstU = max(max(lebesgueFunctionU));
lebesgueConstV = max(max(lebesgueFunctionV));

figure(2)
set(gcf, 'Position', [100,100, 600*2, 600])
subplot(1,2,1)
mesh(XX,YY,lebesgueFunctionU);
title(['Lebesgue Function U, \Lambda_U = ', num2str(lebesgueConstU)])

subplot(1,2,2)
mesh(XX,YY,lebesgueFunctionV);
title(['Lebesgue Function V, \Lambda_V = ', num2str(lebesgueConstV)])
