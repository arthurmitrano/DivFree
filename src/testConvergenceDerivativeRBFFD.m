%% Divergence-free RBF-FD derivatives test
% In this script we test the convergence rate of the derivatives using
% the divergence-free RBF-FD method.

%% Setting up the script
clc, clear

k1 = 3; k2 = 7;  % Control the amout of vortices on the testFunction

rbf = @(e,r) exp(-(e*r).^2);
ep = 3;  % Shape parameter

numPts = 3;   % numPts^2 points in the stencil
nn = 11:10:400; % Works for n odd only

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

RX = rand(3,3);
RY = rand(3,3);
ZZ = ones(3);
ZZ(2,2) = 0;
%%
% theta = pi/7;
theta = 0;
for n = nn
    n, tic
    
    % Generating the center grid ------------------------------------------
    x = linspace(-1,1,n);
    y = linspace(-1,1,n);
    x = x(-floor(numPts/2) + ceil(n/2) : floor(numPts/2) + ceil(n/2));
    y = y(-floor(numPts/2) + ceil(n/2) : floor(numPts/2) + ceil(n/2));
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    [X,Y] = meshgrid(x);
    
    % Perturbing the points
    s = 10;
    X = X - (dx/(2*s) + dx/s*RX).*ZZ;
    Y = Y - (dy/(2*s) + dy/s*RY).*ZZ;
    Xnew = X(:)*cos(theta) + Y(:)*sin(theta);
    Ynew = -X(:)*sin(theta) + Y(:)*cos(theta);
    X(:) = Xnew;
    Y(:) = Ynew;

    dSites = [X(:) Y(:)];
    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    r = DistanceMatrix(dSites, dSites);
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    [u, v, ux, vx, uy, vy, origin] = testFunction(X, Y, k1, k2);
%     f = @(x,y) cos(x^3 + y^3) + sin(x^2 + y^2) + x + y;
%     f = @(x,y) sin(2*(x-.1))*cos(3*(y+.2));
%     [u, v, ux, vx, uy, vy, origin] = testFunction2(X,Y,f);
    uxAtO_Exact = ux;
    uyAtO_Exact = uy;
    vxAtO_Exact = vx;
    vyAtO_Exact = vy;
    t = [u(:) v(:)];
    % ---------------------------------------------------------------------
    
    % Calculating interpolation matrix and its derivatives matrices -------
    [A, F, G, ~, Dx, Dy] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);
    d = reshape(t', 1, numel(t))';
    c = A\d;
    c = reshape(c, 2, size(t,1))';
    % ---------------------------------------------------------------------
    
    % Getting the derivatives at the origin--------------------------------
    Dx_t = reshape(Dx*d, 2, size(t,1))';
    Dy_t = reshape(Dy*d, 2, size(t,1))';
    
    ux = reshape(Dx_t(:,1),size(X));
    uy = reshape(Dy_t(:,1),size(X));
    vx = reshape(Dx_t(:,2),size(Y));
    vy = reshape(Dy_t(:,2),size(Y));
    
    gridCenterX = ceil(size(X,2)/2); gridCenterY = ceil(size(Y,1)/2);
    uxAtO = ux(gridCenterX,gridCenterY);
    uyAtO = uy(gridCenterX,gridCenterY);
    vxAtO = vx(gridCenterX,gridCenterY);
    vyAtO = vy(gridCenterX,gridCenterY);
    % ---------------------------------------------------------------------
    
    % Measuring the error -------------------------------------------------
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    % ---------------------------------------------------------------------
    
    toc
end

%% Plotting
figure(1)
loglog(nn,uxErr,'ro-', nn,uyErr,'bo-', nn,nn.^-2,'b--', nn,nn.^-4,'r--')
axis tight
id = legend('u_x','u_y', 'Location','Best');
set(id, 'FontSize', 12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on u derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)

figure(2)
loglog(nn,vxErr,'ro-', nn,vyErr,'bo-', nn,nn.^-2,'r--', nn,nn.^-4,'b--')
axis tight
id = legend('v_x','v_y', 'Location','Best');
set(id, 'FontSize', 12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on v derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)