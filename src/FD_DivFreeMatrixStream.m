%% FD_DivFreeMatrixStream
% This function creates a matrix $M$ that contains the necessary
% information to create a divergence-free polynomial interpolant of the
% form:
%
% $$u(x,y) = \frac{\partial \psi}{\partial y}, \quad v(x,y) =
% -\frac{\partial \psi}{\partial x},$$
%
% $$\psi(x,y) = \sum_{i=0}^N \sum_{j=0}^N c_{ij}x^iy^j,$$
%
% where $N$ is the degree of the polynomial approximation for the stream
% function.
%
% Using a polynomial of this form, the vector field generated by $(u,v)$
% will be divergence-free. However, numerical experiments showed (so far)
% that the polynomials do not interpolate the test function exactly. A
% least-square solution is find.

%%
function [M, u, v, ux, uy, vx, vy] = ...
         FD_DivFreeMatrixStream(h, N, interpPts, k)
% Generate the interpolation matrix to get the coeffs of the interpolating
% polynomial on a stencil with equally spaced points in both directions 
% with numPts^2 points.
%
% INPUT:
% h          : dx and dy of the stencil
% N          : degree of the bivariate polynomial
% interpPts  : struct with the interpolation points indexes
% k          : order of the derivative to be taken (default = 1)
%
% OUTPUT:
% M              : matrix(h) to be inverted to find the desired polynomial
% u, v           : polynomials that depends on its coefficients.
% ux, uy, vx, vy : depends on (x,y), useful to calculate the numerical
%                  derivatives at a point.
%
% NOTE: the output derivatives might be removed in the future, when we
% decide on the final format of the interpolant.

%% Setting up the function
if (nargin < 4)
    k = 1;
end
numPts = interpPts.numPts;  % The sqrt of the total number of pts

%% Defining the stream function and divergence-free interpolant
% Note that we are using the symbolic toolbox, which make the code slower.
% However, it allows some flexibility for testing. Later we can rewrite the
% code to be computationally efficient.
syms x y

S = kron(x.^(0:N),y.^(0:N));
% Removing redundant terms ------------------------------------------------
S(1) = [];                        % Taking out the constant coefficient
if N == 5
    S(end) = [];                  % Taking out degree 10
    S(end - [5 0]) = [];          % Taking out degree 9
    S(end - [9 4 0]) = [];        % Taking out degree 8
    S(end - [12 7 3 0]) = [];     % Taking out degree 7
    %S(end - [14 9 5 2 0]) = [];  % Taking out degree 6 <== RANK DEFICIENT
    S(end - [5 2 0]) = [];        % Taking out some 6th degree terms
    S(end) = [];                  % Taking a 5th degree term
    S(end - [17 1]) = [];         % Taking out some 4th degree terms
    S(end - 16) = [];             % Taking out a 3rd degree term
end
% -------------------------------------------------------------------------
Sx = diff(S, x); v = -Sx;
Sy = diff(S, y); u = +Sy;   % div(rot( (0,0,S) )) = 0
coeffs = sym('c', [length(u), 1]);

ux = diff(u, x, k); % ux and vy are used to impose div-free condition, and, 
vy = diff(v, y, k); % just like uy and vx, are returned as numerical 
uy = diff(u, y, k); % derivatives.
vx = diff(v, x, k);

ux = matlabFunction(ux*coeffs,'vars',{x y coeffs});
vy = matlabFunction(vy*coeffs,'vars',{x y coeffs});
uy = matlabFunction(uy*coeffs,'vars',{x y coeffs});
vx = matlabFunction(vx*coeffs,'vars',{x y coeffs});

u = matlabFunction(u,'vars',[x y]);
v = matlabFunction(v,'vars',[x y]);

%% Creating numPts x numPts stecil and selecting interpolation points
% The function can only deals with grids that are equally-spaced in both
% directions. Later we might implement a more general function if
% necessary.
%
% The interpolations points are denoted by $(X_u,Y_u)$ for the component
% $u$ and $(X_v,Y_v)$ for the component $v$.
[X, Y] = meshgrid( h*(-floor(numPts/2):1:floor(numPts/2)) );
Xu = X(interpPts.u); Yu = Y(interpPts.u); 
Xv = X(interpPts.v); Yv = Y(interpPts.v);

%% The interpolation matrix is constructed line by line
% At this point, $u$ and $v$ are symbolic row vectors with entries being
% monomials terms of the polynomials that "interpolates" $(f_u,f_v)$. The
% matrix M is constructed such that:
%
%
% $$
% M * \tt{coeffs} = \left[
% \begin{array}{c}
%   U \\
%   V
% \end{array} \right] * \tt{coeffs} = \left[
% \begin{array}{c}
%   f_u(X_u,Y_u) \\
%   f_v(X_v,Y_v)
% \end{array}\right],
% $$
%
% where the matrices $U = [u(X_u,Y_u)]$ and $V = [v(X_v,Y_v)]$ and,
% $(X_u,Y_u)$ and $(X_v,Y_v)$ are a list of interpolation points.
U = zeros(length(Xu(:)),length(S));
V = zeros(length(Xv(:)),length(S)); % There are (N+1)^2 - 1 coeffs

for i = 1:length(Xu(:))
    U(i,:) = u(Xu(i),Yu(i)); % Add interpolation condition for u
end
for i = 1:length(Xv(:))
    V(i,:) = v(Xv(i),Yv(i)); % Add interpolation condition for v
end

% NOTE: Using 2 for loops to be able to use different amount of points for
% each component of the vector field. Probably it will be just 1 loop,
% since symmetry tends to give better results.

M = [U; ...     % This matrix has this particular format to impose the
     V];        % interpolation conditions. Divergence-free condition
                % comes from the stream function now. See below for
                % more information.

% Transforming u and v from symbolic vectors to annonymous functions that
% depends on their coefficients and independent variables.
u = matlabFunction(u(x,y)*coeffs,'vars',{x y coeffs});
v = matlabFunction(v(x,y)*coeffs,'vars',{x y coeffs});
end

%% Comments
% One thing worth mentioning is the fact that, after some numerical
% experiments, we observed that the interpolation condition is not always
% satisfied. The linear system $M*\tt{coeffs} = \mathbf{f}$ is solved
% through least-squares, so the solution that we get is the best
% divergence-free polynomial that minimizes the error on the interpolation
% condition.
%
% Another notable fact is that if we choose the degree $N$ of the
% polynomial to be higher or equal to the number of points $n$ we get
% rank-deficient problem. Matlab displays an warning since there are
% multiple least-squares solution and *returns the one with as few non-zero
% entries as possible*. Moreover, since our polynomial has more degrees of
% freedom to ajust to the vector field to be interpolated, we get *double
% the accuracy* on the derivatives related to the divergence-free condition
% ($u_x$ and $v_y$). See <testConvergenceDerivativeStream.m
% |testConvergenceDerivativeStream.m|>