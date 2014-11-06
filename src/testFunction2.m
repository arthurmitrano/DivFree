%% testFunction2
% This function generates a divergence-free vector field given some
% function $F(x,y)$.
%
%  INPUT:
%  P : list of points where the test function will be evaluated.
%  p : point where we calculate the derivatives of the test function.
%  F : anonymous function (default: @(x,y) x^2 + y^2).
%
%  OUTPUT:
%  u, v   : divergence-free vector field at the points P.
%  u?, v? : derivatives of the divergence-free vector field at point p.

%%
function [u, v, ux, vx, uy, vy] = testFunction2(P,p,F)

%% Setting up the function
if nargin < 3
    F = @(x,y) x.^2 + y.^2;
end

%% Defining a divergence-free vector field
syms x y
u = diff(F(x,y),y);
v = diff(-F(x,y),x);

%% Calculating derivatives
if nargout > 2
    ux = matlabFunction(diff(u,x), 'vars',[x,y]);
    uy = matlabFunction(diff(u,y), 'vars',[x,y]);
    vx = matlabFunction(diff(v,x), 'vars',[x,y]);
    vy = matlabFunction(diff(v,y), 'vars',[x,y]);

    % Evaluating derivatives at p
    ux = ux(p(1),p(2));
    vy = vy(p(1),p(2));
    vx = vx(p(1),p(2));
    uy = uy(p(1),p(2));
end

%% Calculating the divergece-free vector field at grid points P
u = matlabFunction(u, 'vars',[x,y]);
v = matlabFunction(v, 'vars',[x,y]);
u = u(P(:,1),P(:,2));
v = v(P(:,1),P(:,2));