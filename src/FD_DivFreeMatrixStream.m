function [M, u, v, ux, uy, vx, vy] = FD_DivFreeMatrixStream(h, N, m, numPts, ...
                                                uInterpPts, vInterpPts)
% Generate the interpolation matrix to get the coeffs of the interpolating
% polynomial on a stencil with equally spaced points in both directions 
% with numPts^2 points.
% Input:
% N          : degree of the bivariate polynomial
% m          : degree of the univariate polynomial
% numPts     : sqrt of the total number of points of the stecil (default=3)
% uInterpPts : additional interpolation points, additionalPts x 2 matrix
% vInterpPts : additional interpolation points, additionalPts x 2 matrix
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

syms x y

% Defining stream function
S = kron(x.^(0:N),y.^(0:N));
Sx = diff(S, x); Sy = diff(S, y);
u = Sy; v = -Sx;  % div(rot( (0,0,S) )) = 0
%  u(u == 0) = []; v(v == 0) = []; % Taking out the zero coeffs

ux = diff(u, x); % ux and vy are used to impose div-free condition, they 
vy = diff(v, y); % are also returned as output and, just like
uy = diff(u, y); % uy and vx, are returned as output for calculating the 
vx = diff(v, x); % numerical derivatives.

u = matlabFunction(u,'vars',[x y]);
v = matlabFunction(v,'vars',[x y]);

ux = matlabFunction(ux,'vars',[x y]);
vy = matlabFunction(vy,'vars',[x y]);
uy = matlabFunction(uy,'vars',[x y]);
vx = matlabFunction(vx,'vars',[x y]);

% Generates 3x3 stecil
[X, Y] = meshgrid( h*(-(numPts-1)/2:1:(numPts-1)/2) ); 

% The interpolation matrix will be constructed line by line
Du = zeros(length(X(:)),(N+1)*(N+1));
Dv = zeros(length(X(:)),(N+1)*(N+1)); % There are (N+1)*N coeffs

for i = 1:length(X(:))
    Du(i,:) = u(X(i),Y(i)); % Add interpolation condition for @ (X,Y)(i)
    Dv(i,:) = v(X(i),Y(i)); % on both components of the vector field
end

% % Adding interpolation condition for function u, using uInterpPts
% for i = 1:length(uInterpPts)
%     Du = [Du; u(uInterpPts(i,1),uInterpPts(i,2))];
% end
% % Adding interpolation condition for function v, using vInterpPts
% for i = 1:length(vInterpPts)
%     Dv = [Dv; v(vInterpPts(i,1),vInterpPts(i,2))];
% end

M = [Du; ...     % This matrix has this particular format to impose the
     Dv];        % interpolation conditions. Divergence free condition 
                     % comes from the stream function now. See below for 
                     % more information.
% M = matlabFunction(M);

% The idea here is that M*coeffs = f, where
% coeffs = [a(:) b(:)], 
%       with a(:) coeffs of p1 (polynomial for 1st coordnate function) 
%        and b(:) coeffs of p2 (polynomial for 2nd coordnate function)
% f = [f1(X(:),Y(:)) f2(X(:),Y(:)) zeros(length(X(:)))]',
% where f1 and f2 are the coordnate functions that we are interpolating.
% The zeros are added to f to impose the divergence free condition.
end