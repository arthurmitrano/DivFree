%% Symbolic flat limit of the divergence-free RBF interpolant
% The objective of this script is to look to the symbolic expressions of
% the divergence-free RBF and see what happens when $\varepsilon
% \rightarrow 0$. We are looking for:
%
% $$ \lim_{\varepsilon \rightarrow 0} t^\varepsilon(\mathbf{x}) =
% \sum_{k=1}^N \lim_{\varepsilon \rightarrow 0}
% \Psi^\varepsilon(\mathbf{x},\mathbf{x}_k)\mathbf{s}_k.$$
%
% Where
%
% $$\Psi^\varepsilon(\mathbf{x},\mathbf{y}) = -F^\varepsilon(r)I -
% G^\varepsilon(r)(\mathbf{x} - \mathbf{y})(\mathbf{x} - \mathbf{y})^T,$$
%
% with
%
% $$ F^\varepsilon(r) = \frac{\varphi^\prime(r,\varepsilon)}{r} \quad and
% \quad G^\varepsilon(r) = \frac{{F^\varepsilon}^\prime(r)}{r}.$$
%
% Here $r$ is the distance between two centers and $\varphi$ the basic
% function (gaussian, multiquadrics, inverse multiquadrics, etc).
%
% On the end, what really matters is:
% 
% $$ \lim_{\varepsilon \rightarrow 0} F^\varepsilon(r) \quad and \quad
% \lim_{\varepsilon \rightarrow 0} G^\varepsilon(r).$$
%
% The problem with this idea is that we are not taking in consideration the
% coefficients vectors $\mathbf{s}_k$, they actually depend on
% $\varepsilon$. (11/22/2013)

%% Setup for the script
clear, clc
syms ep r

rbf = @(ep,r) exp(-(ep*r).^2);
% rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

%% Calculating the limits
F = diff(rbf(ep,r),r)/r;
G = diff(F,r)/r;
pretty(simplify(F))
limit(F,ep,0)
pretty(simplify(G))
limit(G,ep,0)