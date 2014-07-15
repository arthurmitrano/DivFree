function [u, v, ux, vx, uy, vy, origin] = testFunction(X,Y,k1,k2,k)
% Evaluate the test function and it's derivatives
% INPUT:
% k : order of the derivative (default = 1)
% NOTE: Values of the derivatives at the center of the grid.

if (nargin < 5)
    k = 1;
end

sx = 0.1; sy = 0.2; % Shift in X and Y
% X = X+sx; Y = Y+sy; % Shifted variables

syms x y
u = @(x,y) +1/k1 * sin(k1*(x - sx)) .* cos(k2*(y - sy));
v = @(x,y) -1/k2 * cos(k1*(x - sx)) .* sin(k2*(y - sy));

% Calculating derivatives ----------------------------
ux = matlabFunction(diff(u(x,y), x, k)); ux = ux(X,Y);
vx = matlabFunction(diff(v(x,y), x, k)); vx = vx(X,Y);
uy = matlabFunction(diff(u(x,y), y, k)); uy = uy(X,Y);
vy = matlabFunction(diff(v(x,y), y, k)); vy = vy(X,Y);
% ----------------------------------------------------

% Evaluating derivatives at the origin --------------------------
gridCenterX = ceil(size(X,2)/2); gridCenterY = ceil(size(Y,1)/2);
origin = struct('i',gridCenterX,'j',gridCenterY);
ux = ux(origin.i,origin.j);
vy = vy(origin.i,origin.j);
vx = vx(origin.i,origin.j);
uy = uy(origin.i,origin.j);
% ux = ux(0,0);
% vy = vy(0,0);
% vx = vx(0,0);
% uy = uy(0,0);
% ---------------------------------------------------------------

u = u(X,Y); v = v(X,Y); % Calculating the functions at grid points
end