%% Comparing Symbolic limit with Contour-Padé limit
% The objective of this script is to compare the Symbolic limit found with
% few interpolation points and the one found using the Contour-Padé limit.

%% Setting up the script
clear, clc, close all

%%

rbf = @(ep,r) exp(-(ep*r).^2);        % gaussian
% rbf = @(ep,r) 1./(1 + (ep*r).^2)^2;   % inverse multiquadric (beta = 2)
% rbf = @(ep,r) 1./(1 + (ep*r).^2);     % inverse quadratic - not + def R^2
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2); % IM - not strictly + def. in R^2
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

pointwise = true;  % false: take the limit ep->0 symbolically at (X,Y)
                   % true : evaluate the interp at (XX,YY) and then take
                   %        the limit.

%% Evaluation and data points
X = [0 -1 -1  1 1 1 0 -1  0];
Y = [0 -1  1 -1 1 0 1  0 -1];
U = [1  0  0  0 0 0 0  0  0];
V = [0  0  0  0 0 0 0  0  0];

t = [U(:) V(:)];
d = zeros(2*numel(U),1);
d(1:2:end) = t(:,1);
d(2:2:end) = t(:,2);

dSites = [X(:) Y(:)];    % data points

% xx = linspace(-1,1,40);
% [XX, YY] = meshgrid(xx);
[XX, YY] = chebpts2(32); % using Chebyshev points to evaluate polynomial
ePoints = [XX(:) YY(:)]; % evaluation points

%% Contour-Padé
s = linspace(0,2*pi,1000);  % to plot the circle
rho = 1/(2*sqrt(2)) - .03;  % radius of the circle
e0 = 0;
ep = rho*exp(1i*s) + e0;

r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[~, F, G, Aep] = RBF_DivFreeMatrix(r, d1, d2, rbf, 2); %Just for Aep

r_ePoints = DistanceMatrix(ePoints, dSites);
d1_ePoints = DifferenceMatrix(ePoints(:,1), dSites(:,1));
d2_ePoints = DifferenceMatrix(ePoints(:,2), dSites(:,2));
interp = RBFdivFreeInterpComplex(Aep, d, r_ePoints, d1_ePoints, ...
                                 d2_ePoints, F, G, ep);
interpU = zeros([size(ep) size(XX)]); interpV = interpU;

for i = 1:size(ep,1)
    for j = 1:size(ep,2)
        interpU(i,j,:,:) = reshape(interp.u{i,j}, size(XX));
        interpV(i,j,:,:) = reshape(interp.v{i,j}, size(XX));
    end
end

if (size(ep,1) == 1) || (size(ep,2) == 1)
    interpUatEps0 = 1/(2*pi) * squeeze(trapz(s, interpU));
    interpVatEps0 = 1/(2*pi) * squeeze(trapz(s, interpV));

    figure(1)
    set(gcf, 'Position',[100 100 2*500 1*500])
    
    h1 = subplot(1,2,1);
    mesh(XX,YY,real(interpUatEps0))
    title(['interpolant of U at \epsilon = ' num2str(e0), ', \rho = ' ...
            num2str(rho,'%2.2e')], 'FontSize', 15)
    xlabel('x', 'FontSize',13), ylabel('y', 'FontSize',13)
    zlim([-.5 1])
    
    h2 = subplot(1,2,2);
    mesh(XX,YY,real(interpVatEps0))
    title(['interpolant of V at \epsilon = ' num2str(e0), ', \rho = ' ...
            num2str(rho,'%2.2e')], 'FontSize', 15)
    xlabel('x', 'FontSize',13), ylabel('y', 'FontSize',13)
    zlim([-.5 1])
    
    set([h1 h2], 'clim', [-.5 1])
end

pU = chebfun2(interpUatEps0);
pV = chebfun2(interpVatEps0);
pUcoeffs = chebpoly2(pU);
pVcoeffs = chebpoly2(pV);
ppU = chebpoly2function(pUcoeffs);
ppV = chebpoly2function(pVcoeffs);

syms X Y
disp('ppU')
pretty(expand(sym(ppU(X,Y))))
disp('ppV')
pretty(expand(sym(ppV(X,Y))))
%% Symbolic
syms ep r x1 x2 y1 y2 X Y
assume(X,'real')
assume(Y,'real')
x = [x1; x2];
y = [y1; y2];

x_e = [X;Y];   % evaluation points
x_i = dSites'; % interpolation points

% Constructing divergence-free kernel
F = matlabFunction(diff(rbf(ep,r),r)/r);
G = matlabFunction(diff(F(ep,r),r)/r);
R = sym([0 -1; 1 0]);
K = @(ep,x,y) -F(ep,norm(x-y))*eye(2) ...
              -G(ep,norm(x-y))*R*(x-y)*(x-y).'*R.';
          
A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

disp('taking limit')
if ~pointwise
    L = simplify(limit(simplify(t),ep,0));
    disp('Limit when ep -> 0')
    pretty(L)
    L1 = matlabFunction(L(1),'vars',{X,Y});
    L2 = matlabFunction(L(2),'vars',{X,Y});
else
    epSmall = 1e-6;
    L1 = double(subs( subs(t(1),{X,Y},{XX,YY}), ep,epSmall));
    L2 = double(subs( subs(t(2),{X,Y},{XX,YY}), ep,epSmall));
    disp('Creating polynomial via chebfun2')
    pL1 = chebfun2(L1);
    pL2 = chebfun2(L2);
    pL1coeffs = chebpoly2(pL1);
    pL2coeffs = chebpoly2(pL2);
    L1 = chebpoly2function(pL1coeffs);
    L2 = chebpoly2function(pL2coeffs);

%     disp(['L1 ''limit'' (epSmall = ' num2str(epSmall) ')'])
%     pretty(expand(L1(X,Y)))
%     disp(['L2 ''limit'' (epSmall = ' num2str(epSmall) ')'])
%     pretty(expand(L2(X,Y)))
end

%% Plotting
figure(2)
set(gcf, 'Position', [100 100 2*500 1*500])

h1 = subplot(1,2,1);
mesh(XX,YY,L1(XX,YY))
title('interpolant of U (symbolic)', 'FontSize', 15)
xlabel('x', 'FontSize',13), ylabel('y', 'FontSize',13)
zlim([-.5 1])

h2 = subplot(1,2,2);
mesh(XX,YY,L2(XX,YY))
title('interpolant of V (symbolic)', 'FontSize', 15)
xlabel('x', 'FontSize',13), ylabel('y', 'FontSize',13)
zlim([-.5 1])

set([h1 h2], 'clim', [-.5 1])

%% Mesuaring the difference between both methods
ErrU = max(max(abs(interpUatEps0 - L1(XX,YY))))
ErrV = max(max(abs(interpVatEps0 - L2(XX,YY))))