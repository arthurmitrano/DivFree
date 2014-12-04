%% 3d divergence-free RBF-FD test
% In this script we test the 3d divergence-free RBF-FD method. More
% precisely, we test the convergence rates of the derivatives of the
% divergence-free interpolant.

%% Setting up the script
clc, clear

p = [0 0 0];         % Point to measure the error
ep = 8;              % Shape parameter
randomGrid = false;  % Use a random grid
nn = 11:4:150;       % Points on the unit cube
rbf = @(e,r) exp(-(e*r).^2);

h = zeros(size(nn));  % Fill distance vector

% Derivatives approximation -----------------
DuAtO = struct('x',{[]}, 'y',{[]}, 'z',{[]});
DvAtO = struct('x',{[]}, 'y',{[]}, 'z',{[]});
DwAtO = struct('x',{[]}, 'y',{[]}, 'z',{[]});
% -------------------------------------------

% Derivatives errors ------------------------
DuErr = struct('x',{[]}, 'y',{[]}, 'z',{[]});
DvErr = struct('x',{[]}, 'y',{[]}, 'z',{[]});
DwErr = struct('x',{[]}, 'y',{[]}, 'z',{[]});
% -------------------------------------------

%%
i = 1;
for n = nn
    n, tic
    % Generating the center grid ------------------------------------------
    if randomGrid
        % Random points -------------------
        dSites = [p; -1+2*rand(n^3-1,3)];
        % ---------------------------------
    else
        % Rectangular grid -----------------
        [X, Y, Z] = meshgrid(linspace(-1,1,n));
        dSites = [X(:) Y(:) Z(:)];
        % ----------------------------------
    end
    dSites = nearstNeighbors(dSites, p, 27);

    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    d3 = DifferenceMatrix(dSites(:,3), dSites(:,3));
    r = DistanceMatrix(dSites, dSites);
    R = r.*(diag(inf*ones(size(r,1),1)) + ones(size(r)));
    h(i) = max(min(R));
    i = i + 1;
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    F1 = @(x,y,z) sin((x - 1).^2 + (y - 1).^2 + (z - 1).^2);
    F2 = @(x,y,z) cos((x - .5).^2 + (y - .5).^2 + (z - .5).^2);
    F3 = @(x,y,z) tanh((x + .5).^2 + (y - 1).^2 + (z + 1).^2);
    
    [u, v, w, Du, Dv, Dw] = testFunction3d(dSites, p, F1,F2,F3);
    t = [u(:) v(:) w(:)];
    % ---------------------------------------------------------------------
    
    % Calculating interpolation matrix and its derivatives matrices -------
    [A, Ax, Ay, Az] = rbfDivFreeMatrix3d(r, d1, d2, d3, ep);
    d = [u(:); v(:); w(:)];
    c = A\d;
    % ---------------------------------------------------------------------
    
    % Calculating x-derivatives ---------
    d_x = Ax*c;
    t_x = reshape(d_x, size(t,1), 3);
    DuAtO.x = t_x(1,1);
    DvAtO.x = t_x(1,2);
    DwAtO.x = t_x(1,3);
    % -----------------------------------
    
    % Calculating y-derivatives ---------
    d_y = Ay*c;
    t_y = reshape(d_y, size(t,1), 3);
    DuAtO.y = t_y(1,1);
    DvAtO.y = t_y(1,2);
    DwAtO.y = t_y(1,3);
    % -----------------------------------
    
    % Calculating z-derivatives ---------
    d_z = Az*c;
    t_z = reshape(d_z, size(t,1), 3);
    DuAtO.z = t_z(1,1);
    DvAtO.z = t_z(1,2);
    DwAtO.z = t_z(1,3);
    % -----------------------------------
    
    % Measuring the errors ------------------
    DuErr.x = [DuErr.x; abs(DuAtO.x - Du.x)];
    DuErr.y = [DuErr.y; abs(DuAtO.y - Du.y)];
    DuErr.z = [DuErr.z; abs(DuAtO.z - Du.z)];
    
    DvErr.x = [DvErr.x; abs(DvAtO.x - Dv.x)];
    DvErr.y = [DvErr.y; abs(DvAtO.y - Dv.y)];
    DvErr.z = [DvErr.z; abs(DvAtO.z - Dv.z)];
    
    DwErr.x = [DwErr.x; abs(DwAtO.x - Dw.x)];
    DwErr.y = [DwErr.y; abs(DwAtO.y - Dw.y)];
    DwErr.z = [DwErr.z; abs(DwAtO.z - Dw.z)];
    % --------------------------------------
    toc
end

%% Sorting the error vectors for plotting
[h, I] = sort(h);
DuErr.x = DuErr.x(I);
DuErr.y = DuErr.y(I);
DuErr.z = DuErr.z(I);

DvErr.x = DvErr.x(I);
DvErr.y = DvErr.y(I);
DvErr.z = DvErr.z(I);

DwErr.x = DwErr.x(I);
DwErr.y = DwErr.y(I);
DwErr.z = DwErr.z(I);

%% Rate of decay
p_ux = polyfit(log10(h'),log10(DuErr.x),1);
p_uy = polyfit(log10(h'),log10(DuErr.y),1);
p_uz = polyfit(log10(h'),log10(DuErr.z),1);
p_vx = polyfit(log10(h'),log10(DvErr.x),1);
p_vy = polyfit(log10(h'),log10(DvErr.y),1);
p_vz = polyfit(log10(h'),log10(DvErr.z),1);
p_wx = polyfit(log10(h'),log10(DwErr.x),1);
p_wy = polyfit(log10(h'),log10(DwErr.y),1);
p_wz = polyfit(log10(h'),log10(DwErr.z),1);


fprintf('Rate of decay for ux: h^(%f)\n', p_ux(1))
fprintf('Rate of decay for uy: h^(%f)\n', p_uy(1))
fprintf('Rate of decay for uz: h^(%f)\n', p_uz(1))
fprintf('Rate of decay for vx: h^(%f)\n', p_vx(1))
fprintf('Rate of decay for vy: h^(%f)\n', p_vy(1))
fprintf('Rate of decay for vz: h^(%f)\n', p_vz(1))
fprintf('Rate of decay for wx: h^(%f)\n', p_wx(1))
fprintf('Rate of decay for wy: h^(%f)\n', p_wy(1))
fprintf('Rate of decay for wz: h^(%f)\n', p_wz(1))

%% Plotting
k = floor(length(nn)*(.9));  % index to place h^2 and h^4 on the plot

figure(1)
loglog(h,DuErr.x,'ro', h,DuErr.y,'bs', h,DuErr.z,'k^', ...
       h,h.^2,'c--', h,h.^4,'m--')
xlim([h(1) h(end)])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend(['$$u_x$$: $$h^{', num2str(p_ux(1),3), '}$$'], ...
            ['$$u_y$$: $$h^{', num2str(p_uy(1),3), '}$$'], ...
            ['$$u_z$$: $$h^{', num2str(p_uz(1),3), '}$$'], ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
text(h(k), h(k)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','cyan', 'BackgroundColor','white', 'color','cyan')
text(h(k), h(k)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','magenta', 'BackgroundColor','white', 'color','magenta')
title('Error on $$u$$ derivatives: RBF-FD', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$$h$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(2)
loglog(h,DvErr.x,'ro', h,DvErr.y,'bs', h,DvErr.z,'k^', ...
       h,h.^2,'c--', h,h.^4,'m--')
xlim([h(1) h(end)])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend(['$$v_x$$: $$h^{', num2str(p_vx(1),3), '}$$'], ...
            ['$$v_y$$: $$h^{', num2str(p_vy(1),3), '}$$'], ...
            ['$$v_z$$: $$h^{', num2str(p_vz(1),3), '}$$'], ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
text(h(k), h(k)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','cyan', 'BackgroundColor','white', 'color','cyan')
text(h(k), h(k)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','magenta', 'BackgroundColor','white', 'color','magenta')
title('Error on $$v$$ derivatives: RBF-FD', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$$h$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)

figure(3)
loglog(h,DwErr.x,'ro', h,DwErr.y,'bs', h,DwErr.z,'k^', ...
       h,h.^2,'c--', h,h.^4,'m--')
xlim([h(1) h(end)])
set(gca, 'FontSize', 14)  % Increasing ticks fontsize
id = legend(['$$w_x$$: $$h^{', num2str(p_wx(1),3), '}$$'], ...
            ['$$w_y$$: $$h^{', num2str(p_wy(1),3), '}$$'], ...
            ['$$w_z$$: $$h^{', num2str(p_wz(1),3), '}$$'], ...
            'Location','Best');
set(id, 'Interpreter','latex', 'FontSize',18)
text(h(k), h(k)^2, '$$h^{2}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','cyan', 'BackgroundColor','white', 'color','cyan')
text(h(k), h(k)^4, '$$h^{4}$$', 'Interpreter','latex', 'FontSize',18, ...
     'EdgeColor','magenta', 'BackgroundColor','white', 'color','magenta')
title('Error on $$w$$ derivatives: RBF-FD', 'Interpreter','latex', ...
      'FontSize',20)
xlabel('$$h$$', 'Interpreter','latex', 'FontSize',18)
ylabel('Error', 'Interpreter','latex', 'FontSize',18)