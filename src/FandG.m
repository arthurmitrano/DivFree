%% FandG
% This function calculates the auxiliary functions deduced on the paper
% "Divergence-Free RBFs on Surfaces" from Narcowich, Ward and Wright. In
% terms of equations:
%
% $$
% F(r) = \frac{\phi^{\prime}(r)}{r} \quad {\rm and} \quad
% G(r) = \frac{F^{\prime}(r)}{r}.
% $$
%%
function [F, G, dF, dG] = FandG(rbf)
% Find the expressions for the annonymous functions F, G and it's
% derivatives.
%
% INPUT:
% rbf : radial basis function (annonymous)
%
% OUTPUT:
% F, G   : functions discribed in "Divergence-Free RBFs on Surfaces"
% dF, dG : derivatives of F and G, respectively
%
% NOTE: One can avoid the use of this function by calculating this
% expression once for a particular basic RBF, such as gaussians or
% multiquadrics.

syms ep r

F = simplify( 1/r * diff(rbf(ep,r), r) );
dF = simplify( diff(F,r) );
G = dF/r;
dG = simplify( diff(G,r) );

F = matlabFunction(F);
G = matlabFunction(G);
dF = matlabFunction(dF);
dG = matlabFunction(dG);

end