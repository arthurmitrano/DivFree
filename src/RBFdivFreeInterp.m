%% RBFdivFreeInterporlant
% This function evaluates the RBF divergence-free interpolant on the points
% |ePoints| with interpolation points |dSites| and coefficients |coeffs|.

%%
function t = RBFdivFreeInterp(coeffs, ePoints, dSites, F, G, ep)
% Evaluate the RBF divergence-free interpolant in the evaluation points
% coeffs  : coefficients of the divFree RBF interpolant
% ePoints : evaluation points
% dSites  : interpolation points
% F, G    : functions used to create the kernel and evaluate the 
%           interpolant (anonymous functions)
% ep      : shape parameter of the rbf

r = DistanceMatrix(ePoints, dSites);
d1 = DifferenceMatrix(ePoints(:,1), dSites(:,1));
d2 = DifferenceMatrix(ePoints(:,2), dSites(:,2));

s1 = coeffs(:,1);
s2 = coeffs(:,2);

temp1 = G(ep,r) .* d1;
temp2 = G(ep,r) .* d2;
temp3 = temp2 .* d1;

t(:,1) = -(F(ep,r) + temp1.*d1) * s1 - temp3 * s2;
t(:,2) = -(F(ep,r) + temp2.*d2) * s2 - temp3 * s1;

% OLD OBSERVATION: the result is numerically different when N is large.
% Maybe it is because of the ill condition of the divFree interpolant
% matrix. I DON'T RECALL WHY I MADE THIS OBSERVATION, IT WAS ON THE OLD
% CODE. MIGHT BE HELPFUL FOR DEBUGGING.

end