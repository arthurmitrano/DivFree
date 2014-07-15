%% testFunction2
% This function generates a divergence-free vector field given some
% function $F(x,y)$.
%
%  INPUT:
%  X,Y : matrices of coordinates where the divergence-free vector field
%        will be evaluated.
%  F   : anonymous function (default: @(x,y) x^2 + y^2).
%
%  OUTPUT:
%  u, v   : divergence vector field at the points (X,Y).
%  u?, v? : derivatives at the center of the grid.
%  origin : struct containing the indexes of the center of the grid.
%%
function [u, v, ux, vx, uy, vy, origin] = testFunction2(X,Y,F)

%% Setting up the function
if nargin < 3
    F = @(x,y) x.^2 + y.^2;
end

syms x y

% Defining a divergence-free vector field
u = diff(F(x,y),y);
v = diff(-F(x,y),x);

% Calculating derivatives ------------------------------------------------
ux = matlabFunction(diff(u,x), 'vars',[x,y]); ux = ux(X,Y).*ones(size(X));
uy = matlabFunction(diff(u,y), 'vars',[x,y]); uy = uy(X,Y).*ones(size(X));
vx = matlabFunction(diff(v,x), 'vars',[x,y]); vx = vx(X,Y).*ones(size(X));
vy = matlabFunction(diff(v,y), 'vars',[x,y]); vy = vy(X,Y).*ones(size(X));
% ------------------------------------------------------------------------

% Evaluating derivatives at the origin --------------------------
gridCenterX = ceil(size(X,2)/2); gridCenterY = ceil(size(Y,1)/2);
origin = struct('i',gridCenterX,'j',gridCenterY);
ux = ux(origin.i,origin.j);
vy = vy(origin.i,origin.j);
vx = vx(origin.i,origin.j);
uy = uy(origin.i,origin.j);
% ---------------------------------------------------------------

u = matlabFunction(u, 'vars',[x,y]);
v = matlabFunction(v, 'vars',[x,y]);
u = u(X,Y); v = v(X,Y); % Calculating the functions at grid points
