%% 3d divergence-free RBF-FD test
% In this script we test the 3d divergence-free RBF-FD method. More
% precisely, we test the convergence rates of the derivatives of the
% divergence-free interpolant.

%% Setting up the script
clc, clear

e = 2;  % Shape parameter

numPts = 3;   % numPts^2 points in the stencil
nn = 11:10:200; % Works for n odd only

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

for n = nn
    n, tic
    % Generating the center grid ------------------------------------------
    x = linspace(-1,1,n);
    h = x(2) - x(1);
    x = x(-floor(numPts/2) + ceil(n/2) : floor(numPts/2) + ceil(n/2));
    [X,Y,Z] = meshgrid(x);
    
    % Perturbing the points
%     s = 2;
%     X = X - (h/(2*s) + h/s*RX).*ZZ;
%     Y = Y - (h/(2*s) + h/s*RY).*ZZ;
%     Xnew = X(:)*cos(theta) + Y(:)*sin(theta);
%     Ynew = -X(:)*sin(theta) + Y(:)*cos(theta);
%     X(:) = Xnew;
%     Y(:) = Ynew;

    dSites = [X(:) Y(:) Z(:)];
    d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
    d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
    d3 = DifferenceMatrix(dSites(:,3), dSites(:,3));
    r = DistanceMatrix(dSites, dSites);
    % ---------------------------------------------------------------------
    
    % Getting testFunction values -----------------------------------------
    F1 = @(x,y,z) sin((x-1).^2 + (y-1).^2 + (z-1).^2);
    F2 = @(x,y,z) cos((x-.5).^2 + (y-.5).^2 + (z-.5).^2);
    F3 = @(x,y,z) tanh((x+.5).^2 + (y-1).^2 + (z+1).^2);
    
    [u, v, w, Du, Dv, Dw, O] = testFunction3d(X,Y,Z, F1,F2,F3);
    t = [u(:) v(:) w(:)];
    % ---------------------------------------------------------------------
    
    % Calculating interpolation matrix and its derivatives matrices -------
    [A, Ax, Ay, Az] = rbfDivFreeMatrix3d(r, d1, d2, d3, e);
    d = [u(:); v(:); w(:)];
    c = A\d;
    % ---------------------------------------------------------------------
    
    % Calculating x-derivatives ---------
    d_x = Ax*c;
    t_x = reshape(d_x, size(t,1), 3);
    DuAtO.x = reshape(t_x(:,1), size(X));
    DvAtO.x = reshape(t_x(:,2), size(X));
    DwAtO.x = reshape(t_x(:,3), size(X));
    
    DuAtO.x = DuAtO.x(O.i,O.j,O.k);
    DvAtO.x = DvAtO.x(O.i,O.j,O.k);
    DwAtO.x = DwAtO.x(O.i,O.j,O.k);
    % -----------------------------------
    
    % Calculating y-derivatives ---------
    d_y = Ay*c;
    t_y = reshape(d_y, size(t,1), 3);
    DuAtO.y = reshape(t_y(:,1), size(X));
    DvAtO.y = reshape(t_y(:,2), size(X));
    DwAtO.y = reshape(t_y(:,3), size(X));
    
    DuAtO.y = DuAtO.y(O.i,O.j,O.k);
    DvAtO.y = DvAtO.y(O.i,O.j,O.k);
    DwAtO.y = DwAtO.y(O.i,O.j,O.k);
    % -----------------------------------
    
    % Calculating z-derivatives ---------
    d_z = Az*c;
    t_z = reshape(d_z, size(t,1), 3);
    DuAtO.z = reshape(t_z(:,1), size(X));
    DvAtO.z = reshape(t_z(:,2), size(X));
    DwAtO.z = reshape(t_z(:,3), size(X));
    
    DuAtO.z = DuAtO.z(O.i,O.j,O.k);
    DvAtO.z = DvAtO.z(O.i,O.j,O.k);
    DwAtO.z = DwAtO.z(O.i,O.j,O.k);
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

%% Plotting
figure(1)
loglog(nn,DuErr.x,'ro-', nn,DuErr.y,'bo-', nn,DuErr.z,'k.-', ...
       nn,nn.^-2,'c--', nn,nn.^-4,'m--')
axis tight
id = legend('u_x','u_y','u_z', 'Location','Best');
set(id, 'FontSize',12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on u derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)

figure(2)
loglog(nn,DvErr.x,'ro-', nn,DvErr.y,'bo-', nn,DvErr.z,'k.-', ...
       nn,nn.^-2,'c--', nn,nn.^-4,'m--')
axis tight
id = legend('v_x','v_y','v_z', 'Location','Best');
set(id, 'FontSize',12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on v derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)

figure(3)
loglog(nn,DwErr.x,'ro-', nn,DwErr.y,'bo-', nn,DwErr.z,'k.-', ...
       nn,nn.^-2,'c--', nn,nn.^-4,'m--')
axis tight
id = legend('w_x','w_y','w_z', 'Location','Best');
set(id, 'FontSize',12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on w derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)