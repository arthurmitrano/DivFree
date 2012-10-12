function [tX, tY, DxtX, DxtY, DytX, DytY] = testFunction(X, Y, k1, k2)
% Evaluate the test function and it's derivatives

sx = 0.1; sy = 0.1; % Shift in X and Y
X = X+sx; Y = Y+sy; % Shifted variables

tX = +1/k1 * sin(k1*X) .* cos(k2*Y);
tY = -1/k2 * cos(k1*X) .* sin(k2*Y);
DxtX = +cos(k1*X) .* cos(k2*Y);
DytY = -cos(k1*X) .* cos(k2*Y);
DxtY = +k1/k2 * sin(k1*X) .* sin(k2*Y);
DytX = -k2/k1 * sin(k1*X) .* sin(k2*Y);
end