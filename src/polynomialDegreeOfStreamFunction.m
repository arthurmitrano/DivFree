%% Effect of polynomial degree on accuracy of derivatives
% This script measuares what happen to the error when we increase the
% degree of the polynomial stream function.

%% Setting up the scrip
clc, clear

k1 = 7; k2 = 7;     % control the amout of vortices on the testFunction
NN = 2:1:10;        % Degree of the bivariate polynomial stream function
p = [0 0];          % Point to measure the error
randomGrid = true;  % Use a random grid
n = 191;            % Number of points

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

%% Generating the center grid
if randomGrid
    % Random points -------------------
    dSites = [p; -1+2*rand(n^2-1,2)];
    % ---------------------------------
else
    % Rectangular grid -----------------
    [X, Y] = meshgrid(linspace(-1,1,n));
    dSites = [X(:) Y(:)];
    % ----------------------------------
end
dSites = nearstNeighbors(dSites, p, 9);
r = DistanceMatrix(dSites, dSites);
R = r.*(diag(inf*ones(size(r,1),1)) + ones(size(r)));
h = max(min(R));

%% Getting testFunction values
[u, v, ux, vx, uy, vy] = testFunction(dSites, p, k1, k2);
% f = @(x,y) cos(x^3 + y^3) + sin(x^2 + y^2) + x + y;
% f = @(x,y) sin(2*(x-.1))*cos(3*(y+.2));
% [u, v, ux, vx, uy, vy] = testFunction2(dSites, p, f);
uxAtO_Exact = ux;
uyAtO_Exact = uy;
vxAtO_Exact = vx;
vyAtO_Exact = vy;

%% Main loop
for N = NN
    N, tic
    % Calculating interpolant of degree N ---------------------------------
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrixStream(dSites, N);
    coeffs = M\[u(:); v(:)];  % Get polynomial coefficients
    % ---------------------------------------------------------------------
    
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
    
    % Measuring the error for divFree method ------------------------------
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    % ---------------------------------------------------------------------
    toc
end

%% Plotting
figure(1)
semilogy(NN,uxErr,'ro-', NN,uyErr,'bs-')
xlim([NN(1) NN(end)])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('$u_x$', '$u_y$', 'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)

title(['Error on $$u$$ derivatives: $h = ', num2str(h), '$'], ...
       'Interpreter','latex', 'FontSize',20)
xlabel('$n$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(2)
semilogy(NN,vxErr,'ro-', NN,vyErr,'bs-')
xlim([NN(1) NN(end)])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('$v_x$', '$v_y$', 'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)

title(['Error on $$v$$ derivatives: $h = ', num2str(h), '$'], ...
       'Interpreter','latex', 'FontSize',20)
xlabel('$n$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)