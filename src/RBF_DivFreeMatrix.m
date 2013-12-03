%% RBF_DivFreeMatrix
% This function construct a matrix $A$ that will be used to get the
% divergence-free interpolant.

%%
function [A, F, G, Aep, Dx, Dy] = RBF_DivFreeMatrix(r, d1, d2, rbf, e)
% Construct the interpolation matrix for the RBF divFree interpolant and
% return the divergence-free RBF differentiation matrices (Dx and Dy)
%
% INPUT:
% r    : distance matrix
% d1   : difference matrix on the first coordinate (x?)
% d2   : difference matrix on the second coordinate (y?)
% rbf  : annonymous rbf function
% e    : rbf shape parameter
%
% OUTPUT:
% A      : interpolation matrix. Invert to find interpolant coeffs
% Dx, Dy : divergence-free RBF differentiation matrices
% F, G   : functions used to create the kernel and evaluate the interpolant
% Aep    : interpolation matrix depending on the shape parameter

%% Setting up the function
n = size(r,1); % number of interpolation points
[F, G, dF, dG] = FandG(rbf);
if nargout >= 4
    syms ep
else
    ep = e;
end

%%
% Click <FandG.html here> for more info on the |FandG| function.
%% Creating divergence-free RBF interpolation matrix
A11 = -F(ep,r) - G(ep,r) .* d1.^2;
A22 = -F(ep,r) - G(ep,r) .* d2.^2;
A12 = -G(ep,r) .* d1 .* d2;
% A21 = A12;

% H1 = reshape([A11(:) A12(:)].', 2*n, n).';
% H2 = reshape([A12(:) A22(:)].', 2*n, n).';
% A = reshape([H1 H2].', 2*n, 2*n);
A(1:2:2*n, 1:2:2*n) = A11;
A(1:2:2*n, 2:2:2*n) = A12;
A(2:2:2*n, 1:2:2*n) = A12;
A(2:2:2*n, 2:2:2*n) = A22;

if nargout >= 4
    Aep = matlabFunction(A);
    A = Aep(e);
end


%% Creating differentiation matrices Dx and Dy
if nargout > 4
r0 = r + eye(size(r)); % To avoid division by zero.

temp1 = (-dF(e,r) - dG(e,r) .* d2.^2)./r0;
temp2 = (-dF(e,r) - dG(e,r) .* d1.^2)./r0;
temp3 = dG(e,r) .* d1.*d2 ./r0;

A11 = d1 .* temp1;
A22 = d1 .* ( temp2  - 2*G(e,r) );
A12 = temp3 .* d1 + G(e,r) .* d2;
% A21 = A12;
% H1 = reshape([A11(:) A12(:)].', 2*n, n).';
% H2 = reshape([A12(:) A22(:)].', 2*n, n).';
% Ax = reshape([H1 H2].', 2*n, 2*n);
Ax(1:2:2*n, 1:2:2*n) = A11;
Ax(1:2:2*n, 2:2:2*n) = A12;
Ax(2:2:2*n, 1:2:2*n) = A12;
Ax(2:2:2*n, 2:2:2*n) = A22;

Dx = Ax/A;  % Differentiate in respect to x
end

if nargout > 5
A11 = d2 .* ( temp1 - 2*G(e,r) );
A22 = d2 .* temp2;
A12 = temp3 .* d2 + G(e,r) .* d1;
% A21 = A12;
% H1 = reshape([A11(:) A12(:)].', 2*n, n).';
% H2 = reshape([A12(:) A22(:)].', 2*n, n).';
% Ay = reshape([H1 H2].', 2*n, 2*n);
Ay = zeros(2*n,2*n);
Ay(1:2:2*n, 1:2:2*n) = A11;
Ay(1:2:2*n, 2:2:2*n) = A12;
Ay(2:2:2*n, 1:2:2*n) = A12;
Ay(2:2:2*n, 2:2:2*n) = A22;

Dy = Ay/A;  % Differentiate in respect to y
end

% NOTE: The matrices A11, A22, A12, H1, H2 and the temp variables are not
% important. They are used just for calculations.

end