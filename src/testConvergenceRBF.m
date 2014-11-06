%% testConvergenceRBF
% This script test the convergence of the regular RBF interpolation method
% and the divergence-free RBF interpolation method.

%% Setting up the script
clc, clear

k1 = 7; k2 = 7;       % Control the amout of vortices on the testFunction
ep = 2;               % Shape parameter
randomGrid = false;   % Use a random grid
nn = 7:2:41;        % Points on the unit square
rbf = @(e,r) exp(-(e*r).^2);

h = zeros(size(nn));  % Fill distance vector

ErrRegular = struct('u',{[]}, 'v',{[]});
ErrDivFree = struct('u',{[]}, 'v',{[]});

%%
i = 1;
for n = nn
    n, tic

    % Generating the center grid ------------------------------------------
    if randomGrid
        % Random points -------------
        dSites = -1+2*rand(n^2,2);
        ePoints = -1+2*rand(2*n^2,2);
        % ---------------------------
    else
        % Rectangular grid ---------------------
        [X, Y] = meshgrid(linspace(-1,1,n));
        dSites = [X(:) Y(:)];
        [XX, YY] = meshgrid(linspace(-1,1,2*n));
        ePoints = [XX(:) YY(:)];
        % --------------------------------------
    end

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
    [u, v] = testFunction(dSites, [0 0], k1, k2);
    % f = @(x,y) cos(x^3 + y^3) + sin(x^2 + y^2) + x + y;
    % f = @(x,y) sin(2*(x-.1))*cos(3*(y+.2));
    % [u, v] = testFunction2(dSites, ~, f);
    t = [u(:) v(:)];
    % ---------------------------------------------------------------------

    % Calculating interpolation matrices and interpolants -----------------
    % Divergence-free method ---------------------------------------
    [A_divFree, F, G] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);
    d = reshape(t', 1, numel(t))';
    c_divFree = A_divFree\d;
    c_divFree = reshape(c_divFree, 2, size(t,1))';
    % --------------------------------------------------------------

    % Tradicional method -----------------
    A_regular = rbf(ep,r);
    c_regular = A_regular\[t(:,1) t(:,2)];
    % ------------------------------------
    % ---------------------------------------------------------------------

    % Generating the interpolant in a finner grid -------------------------
    ttDivFree = RBFdivFreeInterp(c_divFree, rr, dd1, dd2, F, G, ep);

    ttRegular = rbf(ep,rr)*c_regular;

    [u, v] = testFunction(ePoints, [0 0], k1, k2);
    ttExact = [u(:) v(:)];
    % ---------------------------------------------------------------------

    % Measuring the error -------------------------------------------------
    ErrDivFree.u = [ErrDivFree.u; max(abs(ttExact(:,1) - ttDivFree(:,1)))];
    ErrDivFree.v = [ErrDivFree.v; max(abs(ttExact(:,2) - ttDivFree(:,2)))];
    ErrRegular.u = [ErrRegular.u; max(abs(ttExact(:,1) - ttRegular(:,1)))];
    ErrRegular.v = [ErrRegular.v; max(abs(ttExact(:,2) - ttRegular(:,2)))];
    % ---------------------------------------------------------------------
    toc
end

%% Plotting
figure(1)
semilogy(nn,ErrDivFree.u,'ro', nn,ErrRegular.u,'bo')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('Divergence-free method','Traditional method', ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
title('Error on $$u$$ component', 'Interpreter','latex', 'FontSize',20)
xlabel('$$\sqrt{N}$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(2)
semilogy(nn,ErrDivFree.v,'ro', nn,ErrRegular.v,'bo')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend('Divergence-free method','Traditional method', ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
title('Error on $$v$$ component', 'Interpreter','latex', 'FontSize',20)
xlabel('$$\sqrt{N}$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)