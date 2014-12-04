%% testConvergenceRBF
% This script test the convergence of the regular RBF interpolation method
% and the divergence-free RBF interpolation method.

%% Setting up the script
clc, clear

k1 = 7; k2 = 7;       % Control the amout of vortices on the testFunction
p = [0 0];            % Point to measure the error on the derivatives
ep = 2;               % Shape parameter
randomGrid = false;   % Use a random grid
nn = 7:2:41;          % Points on the unit square
rbf = @(e,r) exp(-(e*r).^2);

h = zeros(size(nn));  % Fill distance vector

regular = struct('u',{[]}, 'v',{[]}, 'ux',{[]}, 'vx',{[]}, 'uy',{[]}, ...
                 'vy',{[]},'cond',{[]});
divFree = struct('u',{[]}, 'v',{[]}, 'ux',{[]}, 'vx',{[]}, 'uy',{[]}, ...
                 'vy',{[]},'cond',{[]});

%%
i = 1;
for n = nn
    n, tic

    % Generating the grid -------------------------------------------------
    if randomGrid
        % Random points -----------------
        dSites = [p; -1+2*rand(n^2-1,2)];
        ePoints = -1+2*rand(2*n^2,2);
        % -------------------------------
    else
        % Rectangular grid ---------------------
        [X, Y] = meshgrid(linspace(-1,1,n));
        dSites = [X(:) Y(:)];
        [XX, YY] = meshgrid(linspace(-1,1,2*n));
        ePoints = [XX(:) YY(:)];
        % --------------------------------------
    end
    dSites = nearstNeighbors(dSites, p, n^2);

    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    r = DistanceMatrix(dSites, dSites);

    % Finner grid ------------------------------------
    dd1 = DifferenceMatrix(ePoints(:,1), dSites(:,1));
    dd2 = DifferenceMatrix(ePoints(:,2), dSites(:,2));
    rr = DistanceMatrix(ePoints, dSites);
    % ------------------------------------------------
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

    % Calculating interpolation matrices, derivatives and interpolants ----
    % Divergence-free method --------------------------------------------
    [A_divFree, F, G, ~, Ax, Ay] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);
    d = reshape(t', 1, numel(t))';
    c_divFree = A_divFree\d;

    Dx_t = reshape(Ax*c_divFree, 2, size(t,1))';
    Dy_t = reshape(Ay*c_divFree, 2, size(t,1))';

    uxAtO_divFree = Dx_t(1,1);
    uyAtO_divFree = Dy_t(1,1);
    vxAtO_divFree = Dx_t(1,2);
    vyAtO_divFree = Dy_t(1,2);

    c_divFree = reshape(c_divFree, 2, size(t,1))';
    % -------------------------------------------------------------------

    % Tradicional method --------------------------------
    [A_regular, Ax, Ay] = RBF_Matrix(r, d1, d2, rbf, ep);
    c_regular = A_regular\[t(:,1) t(:,2)];

    Dx_t = Ax*c_regular;
    Dy_t = Ay*c_regular;

    uxAtO_regular = Dx_t(1,1);
    uyAtO_regular = Dy_t(1,1);
    vxAtO_regular = Dx_t(1,2);
    vyAtO_regular = Dy_t(1,2);
    % ---------------------------------------------------
    % ---------------------------------------------------------------------

    % Generating the interpolant in a finner grid -------------------------
    ttDivFree = RBFdivFreeInterp(c_divFree, rr, dd1, dd2, F, G, ep);

    ttRegular = rbf(ep,rr)*c_regular;

    [u, v] = testFunction(ePoints, [0 0], k1, k2);
    ttExact = [u(:) v(:)];
    % ---------------------------------------------------------------------

    % Calculating the condition number of interpolation matrices ----------
    divFree.cond = [divFree.cond; cond(A_divFree)];
    regular.cond = [regular.cond; cond(A_regular)];
    % ---------------------------------------------------------------------

    % Measuring the error -------------------------------------------------
    divFree.u = [divFree.u; max(abs(ttExact(:,1) - ttDivFree(:,1)))];
    divFree.v = [divFree.v; max(abs(ttExact(:,2) - ttDivFree(:,2)))];
    divFree.ux = [divFree.ux; abs(uxAtO_divFree - uxAtO_Exact)];
    divFree.uy = [divFree.uy; abs(uyAtO_divFree - uyAtO_Exact)];
    divFree.vx = [divFree.vx; abs(vxAtO_divFree - vxAtO_Exact)];
    divFree.vy = [divFree.vy; abs(vyAtO_divFree - vyAtO_Exact)];

    regular.u = [regular.u; max(abs(ttExact(:,1) - ttRegular(:,1)))];
    regular.v = [regular.v; max(abs(ttExact(:,2) - ttRegular(:,2)))];
    regular.ux = [regular.ux; abs(uxAtO_regular - uxAtO_Exact)];
    regular.uy = [regular.uy; abs(uyAtO_regular - uyAtO_Exact)];
    regular.vx = [regular.vx; abs(vxAtO_regular - vxAtO_Exact)];
    regular.vy = [regular.vy; abs(vyAtO_regular - vyAtO_Exact)];
    % ---------------------------------------------------------------------
    toc
end

%% Plotting error decay
figure(1)
semilogy(nn,divFree.u,'ro', nn,regular.u,'bs')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('Divergence-free method','Traditional method', ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
title('Error on $$u$$ component', 'Interpreter','latex', 'FontSize',20)
xlabel('$$\sqrt{N}$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(2)
semilogy(nn,divFree.v,'ro', nn,regular.v,'bs')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('Divergence-free method','Traditional method', ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
title('Error on $$v$$ component', 'Interpreter','latex', 'FontSize',20)
xlabel('$$\sqrt{N}$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

%% Plotting interpolant and vector field at data sites
figure(3)
quiver(ePoints(:,1),ePoints(:,2), ttDivFree(:,1), ttDivFree(:,2), 'b')
hold on
quiver(dSites(:,1),dSites(:,2), t(:,1), t(:,2), 'r', 'LineWidth', 2.5)
axis([-1 1 -1 1])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
title('Interpolant (blue) and vector field {\bf (red)}', ...
      'Interpreter','latex', 'FontSize',20)
xlabel('$$x$$', 'Interpreter','latex', 'FontSize',18)
ylabel('$$y$$', 'Interpreter','latex', 'FontSize',18)
hold off

%% Printing condition numbers
fprintf('\n\n')
disp('N     DF       R')
fprintf('%2i %4.2e %4.2e\n', [nn' divFree.cond regular.cond]')

%% Print accuracy of derivatives of $u$ (1st component)
fprintf('\n\n')
disp('N   ux-DF     ux-R    uy-DF     uy-R')
fprintf('%2i %4.2e %4.2e %4.2e %4.2e\n', ...
        [nn' divFree.ux regular.ux divFree.uy regular.uy]')