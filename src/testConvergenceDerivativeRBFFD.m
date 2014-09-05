%% Divergence-free RBF-FD derivatives test
% In this script we test the convergence rate of the derivatives using
% the divergence-free RBF-FD method.

%% Setting up the script
clc, clear

k1 = 3; k2 = 7;  % Control the amout of vortices on the testFunction
p = [0 0];       % Point to measure the error
nn = 11:10:200;

rbf = @(e,r) exp(-(e*r).^2);
ep = 2;  % Shape parameter

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

%%
for n = nn
    n, tic
    
    % Generating the center grid ------------------------------------------
    dSites = [p; -1+2*rand(n^2-1,2)];
    dSites = nearstNeighbors(dSites, p, 9);
    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    r = DistanceMatrix(dSites, dSites);
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    [u, v, ux, vx, uy, vy] = testFunction(dSites, p, k1, k2);
%     f = @(x,y) cos(x^3 + y^3) + sin(x^2 + y^2) + x + y;
%     f = @(x,y) sin(2*(x-.1))*cos(3*(y+.2));
%     [u, v, ux, vx, uy, vy] = testFunction2(dSites, p, f);
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
    
    % Getting the derivatives at p ----------------------------------------
    Dx_t = reshape(Dx*d, 2, size(t,1))';
    Dy_t = reshape(Dy*d, 2, size(t,1))';
    
    uxAtO = Dx_t(1,1);
    uyAtO = Dy_t(1,1);
    vxAtO = Dx_t(1,2);
    vyAtO = Dy_t(1,2);
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