%% Testing the convergence of derivatives using divFree FD method

k1 = 3; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 2;           % Degree of the bivariate polynomial
m = 2;           % Degree of the univariate polynomial
numPts = 3;      % The main stecil will have numPts^2 points


uxErr = []; uyErr = []; vxErr = []; vyErr = [];  % Derivative errors

nn = 3:10:200; % The code only works for N odd, will get non-integer index
              % otherwise
for n = nn
    tic
    xx = linspace(-1,1,n);
    yy = linspace(-1,1,n);
    [X, Y] = meshgrid(xx, yy);
    h = abs(xx(2) - xx(1));
    
    [u, v, ux, vx, uy, vy] = testFunction(X, Y, k1, k2);
    uxAtO_Exact = ux( (length(xx) + 1)/2 , (length(yy) + 1)/2 );
    uyAtO_Exact = uy( (length(xx) + 1)/2 , (length(yy) + 1)/2 );
    vxAtO_Exact = vx( (length(xx) + 1)/2 , (length(yy) + 1)/2 );
    vyAtO_Exact = vy( (length(xx) + 1)/2 , (length(yy) + 1)/2 );
    % NOTE: We are using length to get the analytical derivative at zero, 
    % for now that what we want to focus on. Later on, we will generalize
    % this code for different stencil position.
    
    %[dxFD, dyFD] = FD_DivFreeDiff(u, v, h, numPts, degree, extraTerms);
    uInterpPts = []; % Extra interpolation points for component U
    vInterpPts = []; % Extra interpolation points for component V
    
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrix(N, m, numPts, uInterpPts, ...
                                           vInterpPts);
    
    M = M(h);       % Fix M for the stencil 
    Minv = pinv(M); % Calculate the inverse on the least-square sense
    
    idx = (-(numPts-1)/2:(numPts-1)/2); % stencil on x (or y) direction
    I = (length(xx) + 1)/2 + idx; % This will change on the future. Using 
    J = (length(yy) + 1)/2 + idx; % length to get a stencil center at 0
    
    U = u(I,J); V = v(I,J);  % Selects the points of the stencil
    coeff = Minv*[U(:); V(:); zeros([numPts^2,1])];  % Get poly coeffs
    
    % Numerical derivatives (ineficient code, will replace when we decide 
    % on the format of the interpolant):
    numUnknows = length(coeff);
    uxAtO = ux(0,0)*coeff(1:numUnknows/2);
    uyAtO = uy(0,0)*coeff(1:numUnknows/2);
    vxAtO = vx(0,0)*coeff(numUnknows/2 + 1:end);
    vyAtO = vy(0,0)*coeff(numUnknows/2 + 1:end);
    uAtO = t1(0,0)*coeff(1:numUnknows/2);
    vAtO = t2(0,0)*coeff(numUnknows/2 + 1:end);
    % Mesuring the error:
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    toc
    n
end

figure(1)
loglog(nn,uxErr,'r.-', nn,uyErr,'b.-', nn,nn.^-2,'b--', nn,nn.^-4,'r--')
legend('uxErr','uyErr','N^{-2}','N^{-4}','Location','Best')
title('Error on u derivatives')
xlabel('N'), ylabel('Error')
