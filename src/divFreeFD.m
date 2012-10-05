%% Divergence-free differentiation using polynomials

%% Defining a divergence-free vector field
N = 30; % square root of the total number of points
x = linspace(0, 1, N); y = x;  
[X, Y] = meshgrid(x, y);
k1 = 7; k2 = 7;
[tX, tY, DxtX, DxtY, DytX, DytY] = testFunction(X, Y, k1, k2);
figure(1), quiver(X,Y,tX,tY), axis tight

%% Calculating the derivatives using
