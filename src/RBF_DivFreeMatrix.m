%% RBF_DivFreeMatrix
% Construct the interpolation matrix *$A$* for the RBF *divFree*
% interpolant and return divergence-free RBF differentiation matrices
% assistants ($A_x$ and $A_y$), where $D_x = A_x*A^{-1}$ and $D_y =
% A_y*A^{-1}$.
%
%  INPUT:
%  r       : distance matrix
%  d1      : difference matrix on the first coordinate (x?)
%  d2      : difference matrix on the second coordinate (y?)
%  rbf     : annonymous rbf function
%  e       : rbf shape parameter
%  calcAep : calculate symbolic matrix Aep
%
%  OUTPUT:
%  A      : interpolation matrix. Invert to find interpolant coeffs
%  Ax, Ay : divergence-free RBF differentiation matrices
%  F, G   : functions used to create the kernel and evaluate the
%           interpolant
%  Aep    : interpolation matrix depending on the shape parameter

%%
function [A, F, G, Aep, Ax, Ay] = RBF_DivFreeMatrix(r, d1, d2, rbf, e, ...
                                                    calcAep)

%% Setting up the function
n = size(r,1); % number of interpolation points
[F, G, dF, dG] = FandG(rbf);
if nargin < 6
    calcAep = false;
end

if calcAep
    syms ep
else
    ep = e;
end

%%
% Click <FandG.html here> for more info on the |FandG| function.

%% Creating divergence-free RBF interpolation matrix
A11 = -F(ep,r) - G(ep,r) .* d2.^2;
A22 = -F(ep,r) - G(ep,r) .* d1.^2;
A12 = +G(ep,r) .* d1 .* d2;
% A21 = A12;

A(1:2:2*n, 1:2:2*n) = A11;
A(1:2:2*n, 2:2:2*n) = A12;
A(2:2:2*n, 1:2:2*n) = A12;
A(2:2:2*n, 2:2:2*n) = A22;

if calcAep
    Aep = matlabFunction(A);
    A = Aep(e);
else
    Aep = [];
end

%% Creating support matrices Ax and Ay to calculate Dx and Dy
if nargout > 4
    r0 = r + eye(size(r)); % To avoid division by zero.

    temp1 = (-dF(e,r) - dG(e,r) .* d2.^2)./r0;
    temp2 = (-dF(e,r) - dG(e,r) .* d1.^2)./r0;
    temp3 = dG(e,r) .* d1.*d2 ./r0;

    A11 = d1 .* temp1;
    A22 = d1 .* ( temp2  - 2*G(e,r) );
    A12 = temp3 .* d1 + G(e,r) .* d2;
    % A21 = A12;

    Ax = zeros(2*n,2*n);
    Ax(1:2:2*n, 1:2:2*n) = A11;
    Ax(1:2:2*n, 2:2:2*n) = A12;
    Ax(2:2:2*n, 1:2:2*n) = A12;
    Ax(2:2:2*n, 2:2:2*n) = A22;

end

if nargout > 5
    A11 = d2 .* ( temp1 - 2*G(e,r) );
    A22 = d2 .* temp2;
    A12 = temp3 .* d2 + G(e,r) .* d1;
    % A21 = A12;

    Ay = zeros(2*n,2*n);
    Ay(1:2:2*n, 1:2:2*n) = A11;
    Ay(1:2:2*n, 2:2:2*n) = A12;
    Ay(2:2:2*n, 1:2:2*n) = A12;
    Ay(2:2:2*n, 2:2:2*n) = A22;

end

% NOTE: The matrices A11, A22, A12, and the temp variables are not
% important. They are used just for calculations.

end
