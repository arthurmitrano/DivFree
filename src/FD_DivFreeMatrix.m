function M = FD_DivFreeMatrix(N, numPts, extraTerms)
% Generate the interpolation matrix to get the coeffs of the interpolating
% polynomial on a stencil with equally spaced points in both directions 
% with numPts^2 points.
% N          : Degree of the polynomial
% numPts     : sqrt of the total number of points of the stecil (default=3)
% extraTerms : extra polynomial terms (might remove this later)

if (nargin == 1)
    numPts = 3;
end
if (nargin == 3)
    N = N + length(extraTerms);
else
    extraTerms = [];
end

syms x y h

% Interpolating polynomial of the form sum_i(sum_j(a_{ij} * x^i * y^j))
A = [kron(x.^(0:N),y.^(0:N)) extraTerms];
Ax = diff(A, x);
Ay = diff(A, y);  % Ax and Ay are used to impose div-free condition

A = matlabFunction(A);
Ax = matlabFunction(Ax);
Ay = matlabFunction(Ay);

[X, Y] = meshgrid( h*(-(numPts-1)/2:1:(numPts-1)/2) ); % Generates stecil

% The interpolation matrix will be constructed line by line
D = [];
Dx = [];
Dy = [];
for i = 1:length(X(:))
    D = [D; A(X(i),Y(i))];    % Add interpolation condition for @ (X,Y)(i)
    Dx = [Dx; Ax(X(i),Y(i))]; % This last two lines add the divFree 
    Dy = [Dy; Ay(X(i),Y(i))]; % condition.
end
O = zeros(size(D));
M = [D  O; ...     % This matrix has this particular format to impose the
     O  D; ...     % interpolation conditions and the divergence free
     Dx Dy];       % condition. See below for more information.
M = matlabFunction(M);

% The idea here is that M*coeffs = f, where
% coeffs = [a(:) b(:)], 
%       with a(:) coeffs of p1 (polynomial for 1st coordnate function) 
%        and b(:) coeffs of p2 (polynomial for 2nd coordnate function)
% f = [f1(X(:),Y(:)) f2(X(:),Y(:)) zeros(length(X(:)))]',
% where f1 and f2 are the coordnate functions that we are interpolating.
% The zeros are added to f to impose the divergence free condition.
end