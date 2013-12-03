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
% rbf = @(ep,r) exp(-(ep*r).^2);
rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);
ep = sym(1/10);

%% Generating the 2d grid of n^2 equispaced nodes
% xCheb = sort(cos(pi*(0:n-1)/(n-1)));
% if (alpha ~= 0)
%     x = asin(alpha*xCheb)/asin(alpha);
% else
%     x = xCheb;
% end
x = linspace(-1,1,n);
h = x(2) - x(1);  % h is only used if alpha = 1 in FD_DivFreeMatrixStream
xx = linspace(-1,1,11*n); % always use an odd number of points to plot the
                          % interpolation points.

[X,Y] = meshgrid(x);
[XX,YY] = meshgrid(xx);

dSites = sym([X(:) Y(:)]);    % interpolation points
ePoints = sym([XX(:) YY(:)]); % evaluation points

%% Calculating RBF interpolant
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[A, F, G] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);

%%
U = [1 0 0; 0 1 0; 0 0 1]; V = [0 0 0; 0 0 0; 0 0 0];
t = sym([U(:) V(:)]);
d = reshape(t.', 1, numel(t)).';
coeffsSym = A\d;
coeffs = reshape(coeffsSym, 2, size(t,1)).';

interp = RBFdivFreeInterp(coeffs, r, d1, d2, F, G, ep);
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
coeffsP = M\[U(uIdx); V(vIdx)];

resU = max(max(abs(polyInterpU(XX,YY,coeffsP) - double(interpU))))
resV = max(max(abs(polyInterpV(XX,YY,coeffsP) - double(interpV))))

%% Plotting
figure(1)
set(gcf, 'Position', [100,100, 600*2, 600])
h1 = subplot(1,2,1);
mesh(XX,YY,polyInterpU(XX,YY,coeffsP))
axis([-1 1 -1 1 -.5 1])
title('polynomial interpolant for U','FontSize',15)
xlabel('x'), ylabel('y')

h2 = subplot(1,2,2);
mesh(XX,YY,double(interpU))
axis([-1 1 -1 1 -.5 1])
title('RBF interpolant for U','FontSize',15)
xlabel('x'), ylabel('y')

set([h1 h2], 'clim', [-.5 1])

figure(2)
set(gcf, 'Position', [100,100, 600*2, 600])
h3 = subplot(1,2,1);
mesh(XX,YY,polyInterpV(XX,YY,coeffsP))
axis([-1 1 -1 1 -.5 1])
title('polynomial interpolant for V','FontSize',15)
xlabel('x'), ylabel('y')

h4 = subplot(1,2,2);
mesh(XX,YY,double(interpV))
axis([-1 1 -1 1 -.5 1])
title('RBF interpolant for V','FontSize',15)
xlabel('x'), ylabel('y')

set([h3 h4], 'clim', [-.5 1])