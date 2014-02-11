%% Lebesgue functions of divergence-free RBFs interpolants
% This function calculates cardinal and Lebesgue function (global
% interolation). For more information on how the divergence-free RBF
% interpolant is created, see the file <FreeMatrix.html
% |RBF_DivFreeMatrix.m|> and the paper _"Divergence-Free RBFs on Surfaces"_
% of _Narcowich_, _Ward_ and _Wright_.
%
% Calculates Lebesgue function and its constant. The function uses an
% equispaced grid in x- and y-direction.
%
% n       : Number of points in x- and y-direction
% display : displays the cardinal function (optional, default = false)
% ep      : Shape parameter (optional, default = 2)
% rbf     : Radial basis function (optional,
%                                  default = @(ep,r) exp(-(ep*r).^2) )
% alpha   : Parameter for the Kosloff & Tal-Ezer mapping (default = 1);

%%
function lebesgueConst = lebesgueFunctionsRBF(n, display, ep, rbf, alpha)
% Calculates Lebesgue function and its constant. The function uses an
% equispaced grid in x- and y-direction.
%
% n       : Number of points in x- and y-direction
% display : displays the cardinal function (optional, default = false)
% ep      : Shape parameter (optional, default = 2)
% rbf     : Radial basis function (optional,
%                                  default = @(ep,r) exp(-(ep*r).^2) )
% alpha   : Parameter for the Kosloff & Tal-Ezer mapping (default = 1);

%% Setting up the function
if nargin == 1
    display = false;
end
if nargin <= 2
    ep = 2;
end
if nargin <= 3
    rbf = @ (ep,r) exp(-(ep*r).^2);
end
if nargin <= 4
    alpha = 1;
end

%% Generating the 2d grid of n^2 equispaced nodes
xCheb = sort(cos(pi*(0:n-1)/(n-1)));
if (alpha ~= 0)
    x = asin(alpha*xCheb)/asin(alpha);
else
    x = xCheb;
end
xx = linspace(-1,1,7*n);

[X,Y] = meshgrid(x);
[XX,YY] = meshgrid(xx);

dSites = [X(:) Y(:)];   % interpolation points
ePoints = [XX(:) YY(:)]; % evaluation points

%% Calculating cardinal functions
r = DistanceMatrix(dSites, dSites);
d1 = DifferenceMatrix(dSites(:,1), dSites(:,1));
d2 = DifferenceMatrix(dSites(:,2), dSites(:,2));
[A, F, G] = RBF_DivFreeMatrix(r, d1, d2, rbf, ep);

r_ePoints = DistanceMatrix(ePoints, dSites);
d1_ePoints = DifferenceMatrix(ePoints(:,1), dSites(:,1));
d2_ePoints = DifferenceMatrix(ePoints(:,2), dSites(:,2));

L = zeros(n,2*n);  % used to construct delta function

lebesgueFunctionU = zeros(size(XX));
lebesgueFunctionV = zeros(size(XX));

for p = 1:2*n^2  % all interpolation points
    L(p) = 1;
    U = L(:,1:n);
    V = L(:,n+1:2*n);
    t =  [U(:) V(:)];  % vector field
    L(p) = 0;          % going back to zero matrix

    d = reshape(t', 1, numel(t))';
    coeffs = A\d;
    coeffs = reshape(coeffs, 2, size(t,1))';  % I'm not keeping the coeffs
    
    cardFunction = RBFdivFreeInterp(coeffs, r_ePoints, d1_ePoints, ...
                                    d2_ePoints, F, G, ep);
    cardFunctionU = reshape(cardFunction(:,1), size(XX));
    cardFunctionV = reshape(cardFunction(:,2), size(XX));

    [i,j] = ind2sub([n, n], rem(p-1,n^2) + 1); 
    if display
        figure(1)
%         set(gcf, 'Position', [100,100, 600*2, 600])
%         h1 = subplot(1,2,1);
%         mesh(XX,YY,cardFunctionU)
        normcard = sqrt(cardFunctionU.^2+cardFunctionV.^2);
        quiver(XX,YY,cardFunctionU,cardFunctionV,max(max(normcard)))
        %axis tight
        axis([-1 1 -1 1])
        hold on, plot(X,Y,'.k','markersize',20)
        plot(X(i,j),Y(i,j),'or','markersize',20), hold off
%         axis([-1 1 -1 1 -.5 1])
        title(['Cardinal Function U for (i,j) = (' num2str(i) ',' ...
               num2str(j) ')'],'FontSize',15)
    
%         h2 = subplot(1,2,2);
%         mesh(XX,YY,cardFunctionV)
%         axis([-1 1 -1 1 -.5 1])
%         title(['Cardinal Function V for (i,j) = (' num2str(i) ',' ...
%                num2str(j) ')'])
%
%         set([h1 h2], 'clim', [-.5 1])
        snapnow
        pause
    end

    lebesgueFunctionU = lebesgueFunctionU + abs(cardFunctionU);
    lebesgueFunctionV = lebesgueFunctionV + abs(cardFunctionV);
end

%% 
% Above we show the cardinal functions for the case where |n = 3|. Observe
% that the divergence-free interpolant of the zero function is not zero,
% however it interpolates zeros at the interpolation points. This happens
% because the interpolant is divergence-free.

%% Plotting the Lebesgue functions and the Lebesgue constant
lebesgueConstU = max(max(lebesgueFunctionU));
lebesgueConstV = max(max(lebesgueFunctionV));
lebesgueConst = max([lebesgueConstU, lebesgueConstV]);

if display
    figure(2)
    set(gcf, 'Position', [100,100, 600*2, 600])
    subplot(1,2,1)
    mesh(XX,YY,lebesgueFunctionU);
    title(['Lebesgue Function U, \Lambda_U = ', num2str(lebesgueConstU)])

    subplot(1,2,2)
    mesh(XX,YY,lebesgueFunctionV);
    title(['Lebesgue Function V, \Lambda_V = ', num2str(lebesgueConstV)])
end
