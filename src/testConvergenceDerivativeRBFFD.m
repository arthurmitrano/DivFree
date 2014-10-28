%% Divergence-free RBF-FD derivatives test
% In this script we test the convergence rate of the derivatives using
% the divergence-free RBF-FD method.

%% Setting up the script
clc, clear

k1 = 7; k2 = 7;       % Control the amout of vortices on the testFunction
p = [0 0];            % Point to measure the error
nn = 11:10:500;
h = zeros(size(nn));  % Distance of the furthest point from p
rbf = @(e,r) exp(-(e*r).^2);
ep = 2;  % Shape parameter

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

%%
i = 1;
for n = nn
    n, tic
    
    % Generating the center grid ------------------------------------------
    % Random points -------------------
    % dSites = [p; -1+2*rand(n^2-1,2)];
    % ---------------------------------

    % Rectangular grid -----------------
    [X, Y] = meshgrid(linspace(-1,1,n));
    dSites = [X(:) Y(:)];
    % ----------------------------------

    dSites = nearstNeighbors(dSites, p, 9);

    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    r = DistanceMatrix(dSites, dSites);
    h(i) = max(r(1,:));
    i = i + 1;
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    [u, v, ux, vx, uy, vy] = testFunction(dSites, p, k1, k2);
    % f = @(x,y) cos(x^3 + y^3) + sin(x^2 + y^2) + x + y;
    % f = @(x,y) sin(2*(x-.1))*cos(3*(y+.2));
    % [u, v, ux, vx, uy, vy] = testFunction2(dSites, p, f);
    uxAtO_Exact = ux;
    uyAtO_Exact = uy;
    vxAtO_Exact = vx;
    vyAtO_Exact = vy;
    t = [u(:) v(:)];
    % ---------------------------------------------------------------------
    
    % Calculating interpolation matrix and its derivatives matrices -------
    [A, F, G, ~, Ax, Ay] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);
    d = reshape(t', 1, numel(t))';
    c = A\d;
    % c = reshape(c, 2, size(t,1))';
    % ---------------------------------------------------------------------
    
    % Getting the derivatives at p ----------------------------------------
    Dx_t = reshape(Ax*c, 2, size(t,1))';
    Dy_t = reshape(Ay*c, 2, size(t,1))';
    
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
loglog(h,uxErr,'ro', h,uyErr,'bo', h,h.^2,'b--', h,h.^4,'r--')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('$$u_x$$','$$u_y$$', 'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
text(h(5), h(7)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18)
text(h(5), h(6)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18)
title('Error on $$u$$ derivatives', 'Interpreter','latex', 'FontSize',20)
xlabel('$$h$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)


figure(2)
loglog(h,vxErr,'ro', h,vyErr,'bo', h,h.^2,'r--', h,h.^4,'b--')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('$$v_x$$','$$v_y$$', 'Location','Best');
set(id, 'Interpreter','latex', 'FontSize', 18)
text(h(5), h(7)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18)
text(h(5), h(6)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18)
title('Error on $$v$$ derivatives', 'Interpreter','latex', 'FontSize',20)
xlabel('$h$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)