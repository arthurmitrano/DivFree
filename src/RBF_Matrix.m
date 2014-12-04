%% RBF_DivFreeMatrix
% Construct the interpolation matrix *$A$* for the *regular* RBF
% interpolant and return regular RBF differentiation matrices assistants
% ($A_x$ and $A_y$), where $D_x = A_x*A^{-1}$ and $D_y = A_y*A^{-1}$.
%
%  INPUT:
%  r       : distance matrix
%  d1      : difference matrix on the first coordinate (x?)
%  d2      : difference matrix on the second coordinate (y?)
%  rbf     : annonymous rbf function
%  e       : rbf shape parameter
%
%  OUTPUT:
%  A      : interpolation matrix. Invert to find interpolant coeffs
%  Ax, Ay : divergence-free RBF differentiation matrices

%%
function [A, Ax, Ay] = RBF_Matrix(r, d1, d2, rbf, e)
%% Calculating the derivative of the rbf function

syms ep R x
dxrbf = matlabFunction(diff(rbf(ep,R),R) * x/R, 'vars',[ep R x]);

%% Creating RBF interpolation matrix
A = rbf(e,r);

%% Creating support matrices Ax and Ay to calculate Dx and Dy
if nargout > 1
    Ax = dxrbf(e,r,d1);
end

if nargout > 2
    Ay = dxrbf(e,r,d2);
end

end