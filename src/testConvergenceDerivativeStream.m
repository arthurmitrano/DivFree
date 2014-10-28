%% Testing the convergence of derivatives using divFree FD method
clc
close all

k1 = 7; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 3;                % Degree of the bivariate polynomial stream function
p = [0 0];            % Point to measure the error
nn = 11:10:500;
h = zeros(size(nn));  % Distance of the furthest point from p
uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

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
    h(i) = max(DistanceMatrix(dSites,p));
    i = i + 1;
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    [u, v, ux, vx, uy, vy] = testFunction(dSites, p, k1, k2);
    uxAtO_Exact = ux;
    uyAtO_Exact = uy;
    vxAtO_Exact = vx;
    vyAtO_Exact = vy;
    % ---------------------------------------------------------------------
    
    % Not using extra interpolation points anymore!
    
    % Not using interpPts anymore.
    
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrixStream(dSites, N);
    % ---------------------------------------------------------------------

    coeffs = M\[u(:); v(:)];  % Get polynomial coefficients
    
    % Numerical derivatives (inefficient code, will replace when we decide 
    % on the format of the interpolant):
    numUnknows = length(coeffs);
    uxAtO = ux(p(1),p(2),coeffs);
    uyAtO = uy(p(1),p(2),coeffs);
    vxAtO = vx(p(1),p(2),coeffs);
    vyAtO = vy(p(1),p(2),coeffs);
    uAtO = t1(p(1),p(2),coeffs);
    vAtO = t2(p(1),p(2),coeffs);
    % ---------------------------------------------------------------------
    
    % Measuring the error for divFree method -
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    % ----------------------------------------
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