%% Testing the convergence of derivatives using divFree FD method

k1 = 3; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 3;           % Degree of the bivariate polynomial
m = 3;           % Degree of the univariate polynomial
numPts = 3;      % The main stecil will have numPts^2 points


uxErr = []; uyErr = []; vxErr = []; vyErr = []; % Derivative errors
uxErrFD = []; uyErrFD = []; vxErrFD = []; vyErrFD = []; % Trad diff errors

nn = 11:20:200; % The code only works for n odd, will get non-integer index
                % otherwise
for n = nn
    tic
    xx = linspace(-1,1,n);
    yy = linspace(-1,1,n);
    [X, Y] = meshgrid(xx, yy);
    h = abs(xx(2) - xx(1));
    
    [u, v, ux, vx, uy, vy] = testFunction(X, Y, k1, k2);
    
    gridCenterX = (length(xx) + 1)/2; gridCenterY = (length(yy) + 1)/2;
    uxAtO_Exact = ux(gridCenterX, gridCenterY);
    uyAtO_Exact = uy(gridCenterX, gridCenterY);
    vxAtO_Exact = vx(gridCenterX, gridCenterY);
    vyAtO_Exact = vy(gridCenterX, gridCenterY);
    % NOTE: We are using length to get the analytical derivative at zero, 
    % for now that is what we want to focus on. Later, we will generalize
    % this code for different stencil position.
    
    % Extra interpolation points for component U
%     yPts = [gridCenterY-2 gridCenterY+2]; % Add 2 points to the stencil
%     uInterpPts = [xx(gridCenterX)*[1; 1] yy(yPts)']; 
    uInterpPts = [];
    % Extra interpolation points for component V
%     xPts = [gridCenterX-2 gridCenterX+2]; % Add 2 points to the stencil
%     vInterpPts = [xx(xPts)' yy(gridCenterY)*[1; 1]]; 
    vInterpPts = [];
    
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrixStream(h, N, m, numPts, ...
                                                   uInterpPts, vInterpPts);
    
%     M = M(h);       % Fix M for the stencil 
    Minv = pinv(M); % Calculate the inverse on the least-square sense
    
    idx = (-(numPts-1)/2:(numPts-1)/2); % stencil on x (or y) direction
    I = gridCenterX + idx; % This will change on the future. Using 
    J = gridCenterY + idx; % length to get a stencil center at 0
    U = u(I,J); V = v(I,J);  % Selects the points of the stencil
    % Adding extra interpolation conditions
    U = U(:); V = V(:);
%     U = [U; u(yPts,gridCenterX)]; V = [V; v(gridCenterY,xPts)' ];
    
    coeff = M\[U; V;];  % Get poly coeffs
    
    % Numerical derivatives (ineficient code, will replace when we decide 
    % on the format of the interpolant):
    numUnknows = length(coeff);
    uxAtO = ux(0,0)*coeff;
    uyAtO = uy(0,0)*coeff;
    vxAtO = vx(0,0)*coeff;
    vyAtO = vy(0,0)*coeff;
    uAtO = t1(0,0)*coeff;
    vAtO = t2(0,0)*coeff;
    
    % Measuring the error for divFree method:
    uxErr = [uxErr; abs(uxAtO - uxAtO_Exact)];
    uyErr = [uyErr; abs(uyAtO - uyAtO_Exact)];
    vxErr = [vxErr; abs(vxAtO - vxAtO_Exact)];
    vyErr = [vyErr; abs(vyAtO - vyAtO_Exact)];
    
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