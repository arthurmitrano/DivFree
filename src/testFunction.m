%% testFunction
% Generates a test function that has divergence zero, evaluate it at grid
% points and calculate its derivatives at some point.
%
%  INPUT:
%  P      : list of points where the test function will be evaluated.
%  p      : point where we calculate the derivatives of the test function.
%  k1, k2 : controls the amount of vortex in the test function.
%
%  OUTPUT:
%  u, v   : divergence-free vector field at the points P.
%  u?, v? : derivatives of the divergence-free vector field at point p.

%%
function [u, v, ux, vx, uy, vy] = testFunction(P,p,k1,k2)

%% Defining test function
sx = 0.1; sy = 0.2; % Shift in X and Y
syms x y
u = @(x,y) +1/k1 * sin(k1*(x - sx)) .* cos(k2*(y - sy));
v = @(x,y) -1/k2 * cos(k1*(x - sx)) .* sin(k2*(y - sy));

%% Calculating derivatives
ux = matlabFunction(diff(u(x,y), x));
vx = matlabFunction(diff(v(x,y), x));
uy = matlabFunction(diff(u(x,y), y));
vy = matlabFunction(diff(v(x,y), y));

%% Evaluating derivatives at p
ux = ux(p(1),p(2));
uy = uy(p(1),p(2));
vx = vx(p(1),p(2));
vy = vy(p(1),p(2));

%% Calculating the divergece-free vector field at grid points P
u = u(P(:,1),P(:,2));
v = v(P(:,1),P(:,2)); 
end