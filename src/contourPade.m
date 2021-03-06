%% Contour-Pad� for divergence-free RBFs
% In this script first we look to the cardinal functions interpolants on
% the complex plane ($\varepsilon \in \mathbb{C}$). Then we calculate the
% Cachy's integral around a circle to get the value of the interpolant when
% $\varepsilon \rightarrow 0$.

%% Setting up the script
clear, clc

n = 3;      % sqrt root of number of data sites
rbf = @(ep,r) exp(-(ep*r).^2);
% rbf = @(ep,r) 1./(1 + (ep*r).^2)^2;  % beta = 2
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

s = linspace(0,2*pi,1000);   % to plot the circle
rho = 1/(2*sqrt(2)) - 0.03;  % radius of the circle
e0 = 0;
complexPlane = false;
allPoints = false;

if complexPlane
    e = 0.001:0.05/4:1;
    [epX, epY] = meshgrid([-fliplr(e) e]);
    ep = epX + 1i*epY;
else
    ep = rho*exp(1i*s) + e0;
end

%% Generating the 2d grid, finner grid and U and V at data sites
if allPoints
    x = linspace(-1,1,n);
    [X,Y] = meshgrid(x);
    center = ceil(n/2);
    U = zeros(n); U(center, center) = 1;
    V = zeros(n);
else
    X = [0 0 1  0 -1];
    Y = [0 1 0 -1  0];
    U = [1 0 0  0  0];
    V = [0 0 0  0  0];
end
dSites = [X(:) Y(:)];    % data points

if complexPlane
    xx = .75;
else
    xx = linspace(-1,1,11*n);
end
[XX,YY] = meshgrid(xx);
ePoints = [XX(:) YY(:)]; % evaluation points

%% Calculating RBF interpolant
% We are finding the intepolant of the cardinal function where $U \neq 0$
% at center point, all the rest is zero.
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[~, F, G, Aep] = RBF_DivFreeMatrix(r, d1, d2, rbf, 2, true); %Just for Aep


t = [U(:) V(:)];

d = zeros(2*size(t,1),1);
d(1:2:end) = t(:,1);
d(2:2:end) = t(:,2);

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

%% Plotting the interpolant
if (size(ePoints,1) == 1)
%     interpPlotU = ([rot90(interpU,2), flipud(interpU); ...
%                         fliplr(interpU), interpU]);
%     interpPlotV = ([rot90(interpV,2), flipud(interpV); ...
%                         fliplr(interpV), interpV]);
    interpPlotU = interpU; interpPlotV = interpV;
    figure(1)
    set(gcf, 'Position',[50 50 2*700 1*700])
    
    subplot(1,2,1)
    contour(epX,epY,abs(interpPlotU), (0:100))
    set(gca, 'FontSize', 14)  % Increasing ticks fontsize
    zlim([0 10])
    title(['$U^\varepsilon(' num2str(ePoints(1,1)) ', ' ...
           num2str(ePoints(1,2)) ')$ with $\sqrt{N}$ = ' num2str(n)], ...
           'Interpreter','latex', 'FontSize', 20)
    xlabel('$\Re$', 'Interpreter','latex', 'FontSize',18)
    ylabel('$\Im$', 'Interpreter','latex', 'FontSize',18)
    hold on
    plot(rho*cos(s)+e0,rho*sin(s),'k--','LineWidth',2)
    plot([0 0], [-1 1]/(2*sqrt(2)),'m.','MarkerSize',14)
    axis square
    hold off
    
    subplot(1,2,2)
    contour(epX,epY,abs(interpPlotV), (0:100))
    set(gca, 'FontSize', 14)  % Increasing ticks fontsize
    title(['$V^\varepsilon(' num2str(ePoints(1,1)) ', ' ...
           num2str(ePoints(1,2)) ')$ with $\sqrt{N}$ = ' num2str(n)], ...
           'Interpreter','latex', 'FontSize', 20)
    xlabel('$\Re$', 'Interpreter','latex', 'FontSize',18)
    ylabel('$\Im$', 'Interpreter','latex', 'FontSize',18)
    hold on
    plot(rho*cos(s)+e0,rho*sin(s),'k--','LineWidth',2)
    plot([0 0], [-1 1]/(2*sqrt(2)),'m.','MarkerSize',14)
    axis square
    hold off
end
%% Calculating the interpolant when the shape parameter is zero
% For this we will use the Cauchy's intergral formula

if (size(ep,1) == 1) || (size(ep,2) == 1)
    interpUatEps0 = 1/(2*pi) * squeeze(trapz(s, interpU));
    interpVatEps0 = 1/(2*pi) * squeeze(trapz(s, interpV));

    figure(2)
    set(gcf, 'Position',[50 50 2*500 1*300])
    h1 = subplot(1,2,1);
    mesh(XX,YY,real(interpUatEps0))
    title(['Interpolant $u$ at $\varepsilon = ' num2str(e0), '$'], ...
           'Interpreter','latex', 'FontSize',20)
    xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
    ylabel('$y$', 'Interpreter','latex', 'FontSize',18)
    zlim([-.5 1])
    
    h2 = subplot(1,2,2);
    mesh(XX,YY,real(interpVatEps0))
    title(['Interpolant $v$ at $\varepsilon = ' num2str(e0), '$'], ...
           'Interpreter','latex', 'FontSize',20)
    xlabel('$x$', 'Interpreter','latex', 'FontSize',18)
    ylabel('$y$', 'Interpreter','latex', 'FontSize',18)
    zlim([-.5 1])
    
    set([h1 h2], 'clim', [-.5 1])
end
