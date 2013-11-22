%% Contour-Padé for 1-dimension

%% Setting up the script

N = 10; % number of data sites
% rbf = @(ep,r) exp(-(ep*r).^2);
rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

s = linspace(0,2*pi,1000);  % to plot the circle
rho = 0.55;                 % radius of the circle

complexPlane = false;

if (complexPlane)
    e = 0.01:.05:5;
    [epX, epY] = meshgrid(e);
    ep = epX + 1i*epY;
    [epX, epY] = meshgrid([-fliplr(e) e]);
else
    ep = rho*exp(1i*s);
end

% ep = 1;
%% Generating the domain

x = linspace(-1,1,N);
if complexPlane
    xx = 0;
else
    xx = linspace(-1,1,2*N);
end
dSites = x(:);
ePoints = xx(:);

%% Constructing the interpolant
r = DistanceMatrix(dSites,dSites);
rr = DistanceMatrix(ePoints,dSites);
% b = zeros(size(x')); b(floor(N/2)) = 1;
b = sin(pi*x');
% b = 1./(1+(1/16)*x.^2)';
interp = zeros([size(ep) size(xx)]);
for i = 1:size(ep,1)
    disp(['i = ' num2str(i) ' of ' num2str(size(ep,1))])
    for j = 1:size(ep,2)
        A = rbf(ep(i,j),r);
        coeffs = A\b;
        interp(i,j,:,:) = rbf(ep(i,j),rr) * coeffs;
    end
end

%% Plotting the interpolant in the complex plane

if (size(ePoints,1) == 1)
    interpPlot = abs([rot90(interp,2), flipud(interp); ...
                       fliplr(interp), interp]);
    figure(1)
    contour(epX,epY,interpPlot,1:100),shg
    title(['U(' num2str(ePoints(1)) ')'])
    xlabel('Re'), ylabel('Im')
    hold on
    plot(rho*cos(s),rho*sin(s),'k--','LineWidth',2)
    plot(zeros(size(r(2:end,1))),1./rr(1:end,1),'r*','MarkerSize',10)
    plot(zeros(size(r(2:end,1))),-1./rr(1:end,1),'r*','MarkerSize',10)
    hold off
    axis square
end
%% Calculating the interpolant when the shape parameter is zero
% For this we will use the Cauchy's intergral formula

if (size(ep,1) == 1) || (size(ep,2) == 1)
    interpEps0 = 1/(2*pi) * squeeze(trapz(s, interp));

    figure(2)

    plot(xx,real(interpEps0))
    title('interpolant at \epsilon = 0')
    xlabel('x'), ylabel('y')
end