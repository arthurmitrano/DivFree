%% Testing the convergence of derivatives using divFree FD method

k1 = 3; k2 = 7;  % control the amout of vortices on the testFunction

% Setting up the FD_DivFreeMatrix function:
N = 3;           % Degree of the bivariate polynomial
m = 3;           % Degree of the univariate polynomial
numPts = 3;      % The main stecil will have numPts^2 points


uxErr = []; uyErr = []; vxErr = []; vyErr = [];  % Derivative errors

nn = 7:10:100; % The code only works for n odd, will get non-integer index
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
    
    yPts = gridCenterY-2:gridCenterY+2;
    xPts = gridCenterX-2:gridCenterX+2;
    uInterpPts = [xx(gridCenterX)*ones(size((-2:2)')) yy(yPts)']; % Extra interpolation points for component U
    %uInterpPts = [];
    vInterpPts = [xx(gridCenterX)*ones(size((-2:2)')) yy(yPts)']; % Extra interpolation points for component V
    %vInterpPts = [];
    
    [M, t1, t2, ux, uy, vx, vy] = FD_DivFreeMatrix(N, m, numPts, ...
                                                   uInterpPts, vInterpPts);
    
    M = M(h);       % Fix M for the stencil 
    Minv = pinv(M); % Calculate the inverse on the least-square sense
    
    idx = (-(numPts-1)/2:(numPts-1)/2); % stencil on x (or y) direction
    I = gridCenterX + idx; % This will change on the future. Using 
    J = gridCenterY + idx; % length to get a stencil center at 0
    
    U = u(I,J); V = v(I,J);  % Selects the points of the stencil
    U = U(:); V = V(:);
    U = [U; u(yPts,gridCenterX)]; V = [V; v(yPts,gridCenterX) ];
    coeff = Minv*[U; V; zeros([numPts^2,1])];  % Get poly coeffs
    
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

figure(2)
loglog(nn,vxErr,'r.-', nn,vyErr,'b.-', nn,nn.^-2,'r--', nn,nn.^-4,'b--')
legend('vxErr','vyErr','N^{-2}','N^{-4}','Location','Best')
title('Error on v derivatives')
xlabel('N'), ylabel('Error')
