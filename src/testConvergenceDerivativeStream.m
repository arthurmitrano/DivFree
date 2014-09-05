%% Testing the convergence of derivatives using divFree FD method
clc
close all

k1 = 3; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 5;           % Degree of the bivariate polynomial
p = [0 0];       % Point to measure the error

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors

nn = 11:10:200; % The code only works for n odd, will get non-integer index
                % otherwise
for n = nn
    n, tic
    
    % Generating the center grid ------- --
    dSites = [p; -1+2*rand(n^2-1,2)];
    dSites = nearstNeighbors(dSites, p, 9);
    % --------------------------------------------
    
    % Getting testFunction values ------------------------------
    [u, v, ux, vx, uy, vy] = testFunction(dSites, p, k1, k2);
    uxAtO_Exact = ux;
    uyAtO_Exact = uy;
    vxAtO_Exact = vx;
    vyAtO_Exact = vy;
    % ----------------------------------------------------------
    
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

figure(1)
loglog(nn,uxErr,'r.-', nn,uyErr,'b.-', nn,nn.^-2,'b--', nn,nn.^-4,'r--')
id = legend('u_x','u_y', 'Location','Best');
set(id, 'FontSize', 12)
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on u derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)

figure(2)
loglog(nn,vxErr,'r.-', nn,vyErr,'b.-', nn,nn.^-2,'r--', nn,nn.^-4,'b--')
id = legend('v_x','v_y', 'Location','Best');
text(nn(5), nn(5)^-2, 'N^{-2}', 'FontSize', 12, 'FontWeight', 'bold')
text(nn(5), nn(5)^-4, 'N^{-4}', 'FontSize', 12, 'FontWeight', 'bold')
title('Error on v derivatives', 'FontSize', 14)
xlabel('N', 'FontSize', 12), ylabel('Error', 'FontSize', 12)