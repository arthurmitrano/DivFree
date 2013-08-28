%% Limit of divFree RBF interpolant as $\varepsilon$ goes to zero
% This script will show what happens to the divergence-free interpolant as
% the value of $\varepsilon$ goes to zero. We will compare such interpolant
% with the divergence-free interpolant described in <polyBasis.html
% |polyBasis.m|>


%% Setting up the script
clear, clc
n = 3;      % number of points
N = 5;      % degree of the polynomial interpolant
alpha = 1;  % 1: equally-spaced; 0: cheb points
rbf = @(ep,r) exp(-(ep*r).^2);
ep = .01;


%% Generating the 2d grid of n^2 equispaced nodes
xCheb = sort(cos(pi*(0:n-1)/(n-1)));
if (alpha ~= 0)
    x = asin(alpha*xCheb)/asin(alpha);
else
    x = xCheb;
end
h = x(2) - x(1);  % h is only used if alpha = 1 in FD_DivFreeMatrixStream
xx = linspace(-1,1,7*n);

[X,Y] = meshgrid(x);
[XX,YY] = meshgrid(xx);

dSites = [X(:) Y(:)];    % interpolation points
ePoints = [XX(:) YY(:)]; % evaluation points

%% Calculating RBF interpolant
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[A, Dx, Dy, F, G] = RBF_DivFreeMatrix(sym(r), sym(d1), sym(d2), rbf, ep);

U = [0 0 0; 0 1 0; 0 0 0]; V = [0 0 0; 0 0 0; 0 0 0];
t = [U(:) V(:)];
d = reshape(t', 1, numel(t))';
coeffsSym = A\d;
coeffs = double(coeffs);
coeffs = reshape(coeffs, 2, size(t,1))';

interp = RBFdivFreeInterp(coeffs, ePoints, dSites, F, G);
interpU = reshape(interp(:,1), size(XX));
interpV = reshape(interp(:,2), size(XX));

%% Checking if RBF interpolant is close to polynomial
% The idea is to use the divergence-free polynomial and look to the
% residual between those interpolants.

uIdx = find(ones(n));
vIdx = find(ones(n));
interpPts = struct('u', uIdx,'v', vIdx, 'numPts',n); 

[M, polyInterpU, polyInterpV] = FD_DivFreeMatrixStream(h, N, ...
                                                       interpPts, alpha);
coeffsP = M\d;

resU = max(max(abs(polyInterpU(XX,YY,coeffsP) - interpU)));
resV = max(max(abs(polyInterpV(XX,YY,coeffsP) - interpV)));

%% Plotting
set(gcf, 'Position', [100,100, 600*2, 600])
h1 = subplot(1,2,1);
mesh(XX,YY,interpU)
axis([-1 1 -1 1 -.5 1])
title('interpolant for U','FontSize',15)

h2 = subplot(1,2,2);
mesh(XX,YY,interpV)
axis([-1 1 -1 1 -.5 1])
title('interpolant for V','FontSize',15)

set([h1 h2], 'clim', [-.5 1])