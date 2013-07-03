%% Testing the convergence of derivatives using divFree FD method
clc
close all

k1 = 3; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 5;           % Degree of the bivariate polynomial
numPts = 3;      % The main stecil will have numPts^2 points
k = 1;           % Order of the derivative

uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors
uxErrFD = []; uyErrFD = []; vxErrFD = []; vyErrFD = []; % Trad diff errors

nn = 11:20:200; % The code only works for n odd, will get non-integer index
                % otherwise
for n = nn
    tic
    
    % Generating 2d-grid n^2 pts equally-spaced --
    xx = linspace(-1,1,n);  yy = linspace(-1,1,n);
    [X, Y] = meshgrid(xx, yy);
    h = abs(xx(2) - xx(1));
    % --------------------------------------------
    
    % Getting testFunction values ------------------------------
    [u, v, ux, vx, uy, vy, origin] = testFunction(X, Y, k1, k2,k);
    uxAtO_Exact = ux;
    uyAtO_Exact = uy;
    vxAtO_Exact = vx;
    vyAtO_Exact = vy;
    % ----------------------------------------------------------
    
    % Not using extra interpolation points anymore!
    
    % Selecting points to interpolate ------------------------------
    uIdx = find(load('uInterpPts.txt') == 1);
    vIdx = find(load('vInterpPts.txt') == 1); % using linear indexes
    interpPts = struct('u', uIdx,'v', vIdx, 'numPts',numPts); 
    % --------------------------------------------------------------
    
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrixStream(h,N,interpPts,k);
    
%     M = M(h);       % Fix M for the stencil %% MIGHT be useful later on.
    
    % Selecting interpolation points of the stencil -----------------------
    [Iu, Ju] = ind2sub([numPts numPts], uIdx); % Converting linear indexes
    [Iv, Jv] = ind2sub([numPts numPts], vIdx); % to traditional ones
    
    Iu = (Iu-2) + origin.i;
    Ju = (Ju-2) + origin.j;
    Iv = (Iv-2) + origin.i;
    Jv = (Jv-2) + origin.j;
    % NOTE: this only works for a stencil of 3x3, for now. Need to use
    % numPts to generalize.
    U = u(sub2ind([n n], Iu, Ju)); % Converting back to linear 
    V = v(sub2ind([n n], Iv, Jv)); % indexes for selection issues 
    % ---------------------------------------------------------------------
    
    % No more extra interpolation points code!
    
    coeffs = M\[U; V;];  % Get polynomial coefficients
    
    % Numerical derivatives (inefficient code, will replace when we decide 
    % on the format of the interpolant):
    numUnknows = length(coeffs);
    uxAtO = ux(0,0,coeffs);
    uyAtO = uy(0,0,coeffs);
    vxAtO = vx(0,0,coeffs);
    vyAtO = vy(0,0,coeffs);
    uAtO = t1(0,0,coeffs);
    vAtO = t2(0,0,coeffs);
    % ---------------------------------------------------------------------
    
    % Measuring the error for divFree method -
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    % ----------------------------------------
    toc
    n
end

figure(1)
loglog(nn,uxErr,'r.-', nn,uyErr,'b.-', nn,nn.^-2,'b--', nn,nn.^-4,'r--')
legend('uxErr','uyErr','N^{-2}','N^{-4}', 'Location','Best')
title('Error on u derivatives')
xlabel('N'), ylabel('Error')

figure(2)
loglog(nn,vxErr,'r.-', nn,vyErr,'b.-', nn,nn.^-2,'r--', nn,nn.^-4,'b--')
legend('vxErr','vyErr','N^{-2}','N^{-4}', 'Location','Best')
title('Error on v derivatives')
xlabel('N'), ylabel('Error')