%% Testing the convergence of derivatives using divFree FD method
clc, clear

k1 = 7; k2 = 7;     % control the amout of vortices on the testFunction
N = 5;              % Degree of the bivariate polynomial stream function
p = [0 0];          % Point to measure the error
randomGrid = true;  % Use a random grid
nn = 11:10:500;

h = zeros(size(nn));  % Fill distance
uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

i = 1;
for n = nn
    n, tic
    
    % Generating the center grid ------------------------------------------
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
    h(i) = min(max(R));
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

%% Sorting the error vectors for plotting
[h, I] = sort(h);
uxErr = uxErr(I);
uyErr = uyErr(I);
vxErr = vxErr(I);
vyErr = vyErr(I);

%% Rate of decay
p_ux = polyfit(log10(h'),log10(uxErr),1);
p_uy = polyfit(log10(h'),log10(uyErr),1);
p_vx = polyfit(log10(h'),log10(vxErr),1);
p_vy = polyfit(log10(h'),log10(vyErr),1);

fprintf('Rate of decay for ux: h^(%f)\n', p_ux(1))
fprintf('Rate of decay for uy: h^(%f)\n', p_uy(1))
fprintf('Rate of decay for vx: h^(%f)\n', p_vx(1))
fprintf('Rate of decay for vy: h^(%f)\n', p_vy(1))

%% Plotting
k = floor(length(nn)*(.9));  % index to place h^2 and h^4 on the plot

figure(1)
loglog(h,uxErr,'ro', h,uyErr,'bo', h,h.^2,'b--', h,h.^4,'r--')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend(['$$u_x$$: $$h^{', num2str(p_ux(1),3), '}$$'], ...
            ['$$u_y$$: $$h^{', num2str(p_uy(1),3), '}$$'], ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
text(h(k), h(k)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','blue', 'BackgroundColor','white', 'color','blue')
text(h(k), h(k)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','red', 'BackgroundColor','white', 'color','red')
title('Error on $$u$$ derivatives: FD', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$$h$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(2)
loglog(h,vxErr,'ro', h,vyErr,'bo', h,h.^2,'r--', h,h.^4,'b--')
axis tight
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend(['$$v_x$$: $$h^{', num2str(p_vx(1),3), '}$$'], ...
            ['$$v_y$$: $$h^{', num2str(p_vy(1),3), '}$$'], ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize', 18)
text(h(k), h(k)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','red', 'BackgroundColor','white', 'color','red')
text(h(k), h(k)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','blue', 'BackgroundColor','white', 'color','blue')
title('Error on $$v$$ derivatives: FD', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$h$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)