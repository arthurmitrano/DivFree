%% Lebesgue constant growth
% This script has the objective of showing how the Lebesgue constant grows
% as the number of points in the domain increases. We work with the
% polynomial case and the RBF case for different values of $\varepsilon$
% and for the basic function $\varphi(r) = exp(-(\varepsilon r)^2)$.

%% Setting up the script
clear, clc, close all

rbf = @(e,r) exp(-(e*r).^2);
shapeParameters = [2 5 10];
totalPoints = 3:2:13;

%% Calculating the constants for the RBF interpolant
display('RBF divergence-free interpolant')
constRBF = zeros(length(totalPoints),length(shapeParameters));
legends = cell(length(shapeParameters) + 1,1);  % For plotting purposes
j = 1;
for ep = shapeParameters
    i = 1;
    for n = totalPoints
        constRBF(i,j) = lebesgueFunctionsRBF(n, false, ep, rbf);
        i = i + 1;
    end
    legends(j) = {['\epsilon = ', num2str(ep)]};
    j = j + 1;
end

%% Plotting the Lebesgue constants for the RBF interpolant
figure(1)
semilogy(totalPoints,constRBF,'.-', 'MarkerSize',12)
hold on

%%
% We can see here that the shape parameter chages the rate of growth of the
% Lebesgue constant. You can change the parameters and let it run for more
% points, however, when the interpolation matrix starts to become
% ill-conditioned the script runs very slowly.

%% Calculating the constants for the polynomial interpolant
display('Polynomial divergence-free interpolant')
constPoly = zeros(length(totalPoints),1);

i = 1;
for n = totalPoints
    constPoly(i) = lebesgueFunctions(n,n+4);
    i = i + 1;
end

%% Plotting the Lebesgue constant for the polynomial interpolant
semilogy(totalPoints,constPoly,'.--k', 'MarkerSize',12)
legends(end) = {'polynomial'};
title('RBF vs Polynomial divergence-free interpolants')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('Lebesgue constant')
legend(legends, 'Location','NorthWest')
hold off

%%
% For this case we can see that the Lebesgue constant grows faster then the
% RBF divergence-free polynomial. Note that we are not choosing a
% polynomial with minimal degree, we are letting Matlab find the least
% square solution of a rank deficient problem. Matlab returns a solution
% with as few non-zero entries as possible.

%% Optimizing node this distribution for RBF interpolant
% Now we use the Kosloff & Tal-Ezer mapping and the function |fminbnd| to
% find a distribution of points which is the tensor product of the nodes:
%
% $$
% x_j^{kte(\alpha)} := \frac{\arcsin(\alpha x_j^{cheb})}{\arcsin(\alpha)},
% \quad j = 1,\ldots,n,
% $$
%
% where $x_j^{cheb} = cos(\frac{\pi j}{n - 1})$.
display('Optimizing node distribution for RBF interpolant')
minConstRBF = zeros(length(totalPoints),length(shapeParameters));
alphasRBF = zeros(length(totalPoints),length(shapeParameters));

j = 1;
for ep = shapeParameters
    i = 1;
    for n = totalPoints
        [a, L] = fminbnd(@(alpha) lebesgueFunctionsRBF(n, false, ep, ...
                                  rbf, alpha), 0.1, 1);
        alphasRBF(i,j) = a;
        minConstRBF(i,j) = L;
        
        i = i + 1;
    end
    j = j + 1;
end

%% Table of Kosloff & Tal-Ezer parameters and Lebesgue constants (RBF)
rows = num2str(totalPoints); columns = num2str(shapeParameters);
display('minConstRBF')
printmat(minConstRBF,'N\e', rows, columns)
display('alphasRBF')
printmat(alphasRBF,'N\e', rows, columns)

%% Optimizing node this distribution for polynomial interpolant
display('Optimizing node distribution for RBF interpolant')
minConstPoly = zeros(length(totalPoints),1);
alphasPoly = zeros(length(totalPoints),1);

i = 1;
for n = totalPoints
    [a, L] = fminbnd(@(alpha) lebesgueFunctions(n,n+4,n,false,alpha), 0,1);
    alphasPoly(i) = a;
    minConstPoly(i) = L;
    i = i + 1;
end

%% Table of Kosloff & Tal-Ezer parameters and Lebesgue constants (Poly)
rows = num2str(totalPoints); columns = 'LebesgueConst alpha';
printmat([minConstPoly, alphasPoly],'PolyCase', rows, columns)