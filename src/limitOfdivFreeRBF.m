%% Limit of divFree RBF interpolant as $\varepsilon$ goes to zero
% This script will show what happens to the divergence-free interpolant as
% the value of $\varepsilon$ goes to zero. We will compare such interpolant
% with the divergence-free interpolant described in <polyBasis.html
% |polyBasis.m|>


%% Setting up the script
clear, clc
n = 3;      % number of points
N = 9;      % degree of the polynomial interpolant
alpha = 1;  % 1: equally-spaced; 0: cheb points
% rbf = @(ep,r) exp(-(ep*r).^2);
rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);
ep = sym(1/10000000000);

%% Generating the 2d grid of n^2 equispaced nodes
% xCheb = sort(cos(pi*(0:n-1)/(n-1)));
% if (alpha ~= 0)
%     x = asin(alpha*xCheb)/asin(alpha);
% else
%     x = xCheb;
% end
x = linspace(-1,1,n);
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
U = [0 0 0; 0 1 0; 0 0 0]; V = [0 0 0; 0 0 0; 0 0 0];
t = sym([U(:) V(:)]);
d = reshape(t.', 1, numel(t)).';
coeffsSym = A\d;
coeffs = reshape(coeffsSym, 2, size(t,1)).';

r_ePoints = DistanceMatrix(ePoints, dSites);
d1_ePoints = DifferenceMatrix(ePoints(:,1), dSites(:,1));
d2_ePoints = DifferenceMatrix(ePoints(:,2), dSites(:,2));
interp = RBFdivFreeInterp(coeffs, r_ePoints, d1_ePoints, d2_ePoints, F, ...
                          G, ep);
interpU = reshape(interp(:,1), size(XX));
interpV = reshape(interp(:,2), size(XX));

%% Checking if RBF interpolant is close to polynomial
% The idea is to use the divergence-free polynomial and look to the
% residual between those interpolants.
[M, polyInterpU, polyInterpV] = FD_DivFreeMatrixStream(dSites, N);

coeffsP = M\[U(:); V(:)];

resU = max(max(abs(polyInterpU(XX,YY,coeffsP) - double(interpU))))
resV = max(max(abs(polyInterpV(XX,YY,coeffsP) - double(interpV))))

%% Plotting
figure(1)
set(gcf, 'Position', [100,100, 600*2, 600])
h1 = subplot(1,2,1);
mesh(XX,YY,polyInterpU(XX,YY,coeffsP))
axis([-1 1 -1 1 -.5 1])
title('Polynomial interpolant for $u$', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
ylabel('$y$', 'Interpreter','latex', 'FontSize',18)

h2 = subplot(1,2,2);
mesh(XX,YY,double(interpU))
axis([-1 1 -1 1 -.5 1])
title('RBF interpolant for $u$', 'Interpreter','latex', 'FontSize',20)
xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
ylabel('$y$', 'Interpreter','latex', 'FontSize',18)

set([h1 h2], 'clim', [-.5 1])

figure(2)
set(gcf, 'Position', [100,100, 600*2, 600])
h3 = subplot(1,2,1);
mesh(XX,YY,polyInterpV(XX,YY,coeffsP))
axis([-1 1 -1 1 -.5 1])
title('Polynomial interpolant for $v$', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
ylabel('$y$', 'Interpreter','latex', 'FontSize',18)

h4 = subplot(1,2,2);
mesh(XX,YY,double(interpV))
axis([-1 1 -1 1 -.5 1])
title('RBF interpolant for $v$', 'Interpreter','latex', 'FontSize',20)
xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
ylabel('$y$', 'Interpreter','latex', 'FontSize',18)

set([h3 h4], 'clim', [-.5 1])

%%
% As we can see from the plots the RBF interpolant seems to converge to
% some function that is different than the divergence-free polynomial
% interpolant. However, this doesn't mean that the limit of the
% divergence-free RBFs is not a polynomial, since it might be still a
% polynomial just of a different form than how we are constructing our
% divergence-free polyomial.