function [M, u, v, ux, uy, vx, vy] = FD_DivFreeMatrix(N, m, numPts, ...
                                                uInterpPts, vInterpPts)
% Generate the interpolation matrix to get the coeffs of the interpolating
% polynomial on a stencil with equally spaced points in both directions 
% with numPts^2 points.
% Input:
% N          : degree of the bivariate polynomial
% m          : degree of the univariate polynomial
% numPts     : sqrt of the total number of points of the stecil (default=3)
% uInterpPts : additional interpolation points, totalPts x 2 matrix
% vInterpPts : additional interpolation points, totalPts x 2 matrix
% Output: all annonymous functions
% M              : matrix(h) to be inverted to find the desired polynomial
% ux, uy, vx, vy : depends on (x,y), useful to calculate the numerical
%                  derivatives at a point. Vector containing the monomials
% NOTE: the output derivatives might be removed in the future, when we
% decide on the final format of the interpolant.

if (nargin < 3)
    numPts = 3;
end
if (nargin < 4)
    uInterpPts = [];
    vInterpPts = [];
end

syms x y h

% Interpolating polynomial of the form sum_i(sum_j(a_{ij} * x^i * y^j))
u = [kron(x.^(0:N),y.^(1:N)) y.^(0:m)];
v = [kron(x.^(1:N),y.^(0:N)) x.^(0:m)];
ux = diff(u, x); % ux and vy are used to impose div-free condition, they 
vy = diff(v, y); % are also returned as output and, just like
uy = diff(u, y); % uy and vx, are returned as output for calculating the 
vx = diff(v, x); % numerical derivatives.

u = matlabFunction(u);
v = matlabFunction(v);
ux = matlabFunction(ux);
vy = matlabFunction(vy);
uy = matlabFunction(uy);
vx = matlabFunction(vx);

[X, Y] = meshgrid( h*(-(numPts-1)/2:1:(numPts-1)/2) ); % Generates stecil

% The interpolation matrix will be constructed line by line
Du = [];
Dv = [];
Dux = [];
Dvy = [];
for i = 1:length(X(:))
    Du = [Du; u(X(i),Y(i))];   % Add interpolation condition for @ (X,Y)(i)
    Dv = [Dv; v(X(i),Y(i))];   % on both components of the vector field
    Dux = [Dux; ux(X(i),Y(i))]; % This last two lines add the divFree 
    Dvy = [Dvy; vy(X(i),Y(i))]; % condition.
end
% Adding interpolation condition for function u, using uInterpPts
for i = 1:length(uInterpPts)
    Du = [Du; u(uInterpPts(i,1),uInterpPts(i,2))];
end
% Adding interpolation condition for function v, using vInterpPts
for i = 1:length(vInterpPts)
    Dv = [Dv; v(vInterpPts(i,1),vInterpPts(i,2))];
end

Ou = zeros(size(Du));
Ov = zeros(size(Dv));
M = [Du  Ou; ...     % This matrix has this particular format to impose the
     Ov  Dv; ...     % interpolation conditions and the divergence free
     Dux Dvy];       % condition. See below for more information.
M = matlabFunction(M);

% The idea here is that M*coeffs = f, where
% coeffs = [a(:) b(:)], 
%       with a(:) coeffs of p1 (polynomial for 1st coordnate function) 
%        and b(:) coeffs of p2 (polynomial for 2nd coordnate function)
% f = [f1(X(:),Y(:)) f2(X(:),Y(:)) zeros(length(X(:)))]',
% where f1 and f2 are the coordnate functions that we are interpolating.
% The zeros are added to f to impose the divergence free condition.
end