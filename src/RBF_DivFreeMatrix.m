%% RBF_DivFreeMatrix
% This function construct a matrix $A$ that will be used to get the
% divergence-free interpolant.

%%
function [A, Dx, Dy, F, G] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep)
% Construct the interpolation matrix for the RBF divFree interpolant and
% return the divergence-free RBF differentiation matrices (Dx and Dy)
%
% INPUT:
% r    : distance matrix
% d1   : difference matrix on the first coordinate (x?)
% d2   : difference matrix on the second coordinate (y?)
% rbf  : annonymous rbf function
% ep   : shape parameter of the rbf function
%
% OUTPUT:
% A      : interpolation matrix to be inverted to find interpolant coeffs
% Dx, Dy : divergence-free RBF differentiation matrices
% F, G   : functions used to create the kernel and evaluate the interpolant

%% Setting up the function
n = size(r,1); % number of interpolation points
[F, G, dF, dG] = FandG(ep, rbf);

%%
% Click <FandG.html here> for more info on the |FandG| function.
%% Creating divergence-free RBF interpolation matrix
A11 = -F(r) - G(r) .* d2.^2; 
A22 = -F(r) - G(r) .* d1.^2; 
A12 = +G(r) .* d1 .* d2;
% A21 = A12;

H1 = reshape([A11(:) A12(:)]', 2*n, n)';
H2 = reshape([A12(:) A22(:)]', 2*n, n)';
A = reshape([H1 H2]', 2*n, 2*n);

%% Creating differentiation matrices Dx and Dy
r0 = r + eye(size(r)); % To avoid division by zero.

temp1 = (-dF(r) - dG(r) .* d2.^2)./r0; 
temp2 = (-dF(r) - dG(r) .* d1.^2)./r0; 
temp3 = dG(r) .* d1.*d2 ./r0;

A11 = d1 .* temp1;
A22 = d1 .* ( temp2  - 2*G(r) );
A12 = temp3 .* d1 + G(r) .* d2;
% A21 = A12;
H1 = reshape([A11(:) A12(:)]', 2*n, n)';
H2 = reshape([A12(:) A22(:)]', 2*n, n)';
Ax = reshape([H1 H2]', 2*n, 2*n);

Dx = Ax/A;  % Differentiate in respect to x

A11 = d2 .* ( temp1 - 2*G(r) );
A22 = d2 .* temp2;
A12 = temp3 .* d2 + G(r) .* d1;
% A21 = A12;
H1 = reshape([A11(:) A12(:)]', 2*n, n)';
H2 = reshape([A12(:) A22(:)]', 2*n, n)';
Ay = reshape([H1 H2]', 2*n, 2*n);

Dy = Ay/A;  % Differentiate in respect to y

% NOTE: The matrices A11, A22, A12, H1, H2 and the temp variables are not
% important. They are used just for calculations.

end