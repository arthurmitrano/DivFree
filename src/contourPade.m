%% Contour-Padé for divergence-free RBFs
% In this script first we look to the cardinal functions interpolants on
% the complex plane ($\varepsilon \in \mathbb{C}$). Then we calculate the
% Cachy's integral around a circle to get the value of the interpolant when
% $\varepsilon \rightarrow 0$.

%% Setting up the script
clear, clc
n = 3;      % sqrt root of number of data sites
alpha = 1;  % Kosloff & Tal-Ezer map parameter (1- equispaced; 0- Cheb)
% rbf = @(ep,r) exp(-(ep*r).^2);
rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

s = linspace(0,2*pi,100);  % to plot the circle
rho = 1/(2*sqrt(2)) -0.01;                 % radius of the circle
e0 = 1/10;
complexPlane = false;

if complexPlane
    e = 0.01:0.1/4:1;
%     [epX, epY] = meshgrid(e);
    [epX, epY] = meshgrid([-fliplr(e) e]);
    ep = epX + 1i*epY;
else
    ep = rho*exp(1i*s) + e0;
end

%% Generating the 2d grid of n^2 equispaced nodes
% xCheb = sort(cos(pi*(0:n-1)/(n-1)));
% if (alpha ~= 0)
%     x = asin(alpha*xCheb)/asin(alpha);
% else
%     x = xCheb;
% end
x = linspace(-1,1,n);
if complexPlane
    xx = 0.75;
else
    xx = linspace(-1,1,11*n);
end

[X,Y] = meshgrid(x);
[XX,YY] = meshgrid(xx);

dSites = [X(:) Y(:)];    % interpolation points
ePoints = [XX(:) YY(:)]; % evaluation points

%% Calculating RBF interpolant
% We are finding the intepolant of the cardinal function where $U \neq 0$
% at center point, all the rest is zero.
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[A, Dx, Dy, F, G, Aep] = RBF_DivFreeMatrix(r, d1, d2, rbf, 2);%Just for Aep

U = [0 0 0; 0 1 0; 0 0 0];
V = zeros(3,3);

t = [U(:) V(:)];
d = reshape([U(:) V(:)]',1,numel(t))';

interp = RBFdivFreeInterpComplex(Aep, d, ePoints, dSites, F, G, ep);
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
    
    subplot(1,2,1)
    contour(epX,epY,abs(interpPlotU), (0:100))
    zlim([0 10])
    title(['U(' num2str(ePoints(1,1)) ', ' num2str(ePoints(1,2)) ')'])
    xlabel('Re'), ylabel('Im')
    hold on
    plot(rho*cos(s)+e0,rho*sin(s),'k--','LineWidth',2)
    plot([0 0], [-1 1]/sqrt(2),'r*','MarkerSize',10)
    hold off
    axis square
    
    subplot(1,2,2)
    contour(epX,epY,abs(interpPlotV), (0:100))
    title(['V(' num2str(ePoints(1,1)) ', ' num2str(ePoints(1,2)) ')'])
    xlabel('Re'), ylabel('Im')
    hold on
    plot(rho*cos(s)+e0,rho*sin(s),'k--','LineWidth',2)
    plot([0 0], [-1 1]/sqrt(2),'r*','MarkerSize',10)
    hold off
    axis square
end
%% Calculating the interpolant when the shape parameter is zero
% For this we will use the Cauchy's intergral formula


if (size(ep,1) == 1) || (size(ep,2) == 1)
    interpUatEps0 = 1/(2*pi) * squeeze(trapz(s, interpU));
    interpVatEps0 = 1/(2*pi) * squeeze(trapz(s, interpV));

    figure(20)
    
    h1 = subplot(1,2,1);
    mesh(XX,YY,real(interpUatEps0))
    title(['interpolant of U at \epsilon = ' num2str(e0), ', \rho = ' ...
            num2str(rho,'%2.2e')])
    xlabel('x'), ylabel('y')
    zlim([-.5 1])
    
    h2 = subplot(1,2,2);
    mesh(XX,YY,real(interpVatEps0))
    title(['interpolant of V at \epsilon = ' num2str(e0), ', \rho = ' ...
            num2str(rho,'%2.2e')])
    xlabel('x'), ylabel('y')
    zlim([-.5 1])
    
    set([h1 h2], 'clim', [-.5 1])
end
