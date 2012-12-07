function [M, u, v, ux, uy, vx, vy] = ...
         FD_DivFreeMatrixStream(h, N, interpPts, k)
% Generate the interpolation matrix to get the coeffs of the interpolating
% polynomial on a stencil with equally spaced points in both directions 
% with numPts^2 points.
% INPUT:
% h          : dx and dy of the stencil
% N          : degree of the bivariate polynomial
% interpPts  : struct with the interpolation points indexes
% k          : order of the derivative to be taken (default = 1)
% OUTPUT: all annonymous functions
% M              : matrix(h) to be inverted to find the desired polynomial
% ux, uy, vx, vy : depends on (x,y), useful to calculate the numerical
%                  derivatives at a point. Vector containing the monomials
% NOTE: the output derivatives might be removed in the future, when we
% decide on the final format of the interpolant.

if (nargin < 4)
    k = 1;
end
numPts = interpPts.numPts;  % The sqrt of the total number of pts

syms x y

% Defining stream function and divFree interpolant
S = kron(x.^(0:N),y.^(0:N));
S([1]) = []; % Taking out the constant coefficient
Sx = diff(S, x); v = -Sx;
Sy = diff(S, y); u = +Sy;   % div(rot( (0,0,S) )) = 0

ux = diff(u, x, k); % ux and vy are used to impose div-free condition, and, 
vy = diff(v, y, k); % just like uy and vx, are returned as numerical 
uy = diff(u, y, k); % derivatives.
vx = diff(v, x, k);

u = matlabFunction(u,'vars',[x y]);
v = matlabFunction(v,'vars',[x y]);

ux = matlabFunction(ux,'vars',[x y]);
vy = matlabFunction(vy,'vars',[x y]);
uy = matlabFunction(uy,'vars',[x y]);
vx = matlabFunction(vx,'vars',[x y]);

% Generates 3x3 stecil and selecting interpolation points
[X, Y] = meshgrid( h*(-(numPts-1)/2:1:(numPts-1)/2) ); 
Xu = X(interpPts.u); Yu = Y(interpPts.u); 
Xv = X(interpPts.v); Yv = Y(interpPts.v);
% -------------------------------------------------------

% The interpolation matrix will be constructed line by line ---------------
Du = zeros(length(Xu(:)),length(S));
Dv = zeros(length(Xv(:)),length(S)); % There are (N+1)^2 - 1 coeffs

for i = 1:length(Xu(:))
    Du(i,:) = u(Xu(i),Yu(i)); % Add interpolation condition for u
end
for i = 1:length(Xv(:))
    Dv(i,:) = v(Xv(i),Yv(i)); % Add interpolation condition for v
end
% NOTE: Using 2 for loops to be able to use different amount of points for
% each component of the vector field. Probably it will be just 1 loop,
% since symetry tends to give better accuracy.

M = [Du; ...     % This matrix has this particular format to impose the
     Dv];        % interpolation conditions. Divergence free condition 
                 % comes from the stream function now. See below for 
                 % more information.
% M = matlabFunction(M);
% -------------------------------------------------------------------------


% The idea here is that M*coeffs = f, where
% coeffs = [a(:) b(:)], 
%       with a(:) coeffs of p1 (polynomial for 1st coordnate function) 
%        and b(:) coeffs of p2 (polynomial for 2nd coordnate function)
% f = [f1(X(:),Y(:)) f2(X(:),Y(:)) zeros(length(X(:)))]',
% where f1 and f2 are the coordnate functions that we are interpolating.
% The zeros are added to f to impose the divergence free condition.
end