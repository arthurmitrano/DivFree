%% Lebesgue constant growth
% This script has the objective of showing how the Lebesgue constant grows
% as the number of points in the domain increases. We work with the
% polynomial case and the RBF case for different values of $\varepsilon$
% and for the basic function $\varphi(r) = exp(-(\varepsilon r)^2)$.

%% Setting up the script
clear, clc, close all

rbf = @(e,r) exp(-(e*r).^2);
shapeParameters = [2 5 10];
totalPoints = 3:2:11;

%% Calculating the constants for the RBF interpolant
display('RBF divergence-free interpolant')
constRBF = zeros(length(totalPoints),length(shapeParameters));
figure(1)
legends = cell(length(shapeParameters),1);
j = 1;
for ep = shapeParameters
    i = 1;
    for n = totalPoints
        [constU, constV] = lebesgueFunctionsRBF(n, false, ep, rbf);
        constRBF(i,j) = max([constU, constV]);
        i = i + 1;
    end
    legends(j) = {['\epsilon = ', num2str(ep)]};
    j = j + 1;
end

%% Plotting the Lebesgue constants for the RBF interpolant
semilogy(totalPoints,constRBF,'.-', 'MarkerSize',12)
title('Lebesgue constant - RBF divergence-free interpolant')
xlabel('N')
ylabel('Lebesgue constant')
legend(legends, 'Location','NorthWest')

%%
% We can see here that the shape parameter chages the rate of growth of the
% Lebesgue constant. You can change the parameters and let it run for more
% points, however, when the interpolation matrix starts to become
% ill-conditioned the script runs very slowly.

%% Calculating the constants for the polynomial interpolant
display('Polynomial divergence-free interpolant')
constPoly = zeros(length(totalPoints),1);
figure(2)
i = 1;
for n = totalPoints
    [constU, constV] = lebesgueFunctions(n,n+4);
    constPoly(i) = max([constU, constV]);
    i = i + 1;
end

%% Plotting the Lebesgue constant for the polynomial interpolant
semilogy(totalPoints,constPoly,'.-', 'MarkerSize',12)
title('Lebesgue constant - Polynomial divergence-free interpolant')
xlabel('N')
ylabel('Lebesgue constant')

%%
% For this case we can see that the Lebesgue constant grows faster then the
% RBF divergence-free polynomial. Note that we are not choosing a
% polynomial with minimal degree, we are letting Matlab find the least
% square solution of a rank deficient problem. Matlab returns a solution
% with as few non-zero entries as possible.