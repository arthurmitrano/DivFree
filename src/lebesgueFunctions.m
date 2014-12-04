%% Lebesgue functions of divergence-free polynomial interpolants
% The purpose of this script is to visualize the cardinal and Lebesgue
% functions of the *polynomial* divergence-free method based on the *stream
% function* idea. The user can change the degree $N$ of the polynomial
% interpolant $(u_{interp}, v_{interp})$. For more information on how this
% polynomial is constructed see <FD_DivFreeMatrixStream.html
% FD_divFreeMatrixStream.html>.
%
%  Calculates Lebesgue function and it's constant.
%
%  n       : Number of grid points in one dimension
%  N       : Degree of the bivariate interpolation polynomial
%  display : Display the cardinal function plots (optinal)
%  alpha   : Kosloff & Tal-Ezer parameter (default = 1)
%
%  NOTE: n must be odd numbers.
%%
function lebesgueConst = lebesgueFunctions(n,N,display,alpha)
%  Calculates Lebesgue function and it's constant.
%
%  n       : Number of grid points in one dimension
%  N       : Degree of the bivariate interpolation polynomial
%  display : Display the cardinal function plots (optinal)
%  alpha   : Kosloff & Tal-Ezer parameter (default = 1)
%
%  NOTE: n must be odd numbers.

%% Setting up the function
if (nargin <= 2)
    display = false;
end

if (nargin < 4)
    alpha = 1;
end

%% Generating 2d-grid of n^2 pts equally-spaced or mapped via K&T
xCheb = sort(cos(pi*(0:n-1)/(n-1)));
if (alpha ~= 0)
    xx = asin(alpha*xCheb)/asin(alpha);
else
    xx = xCheb;
end
h = abs(xx(2) - xx(1));
[X,Y] = meshgrid(xx);
dSites = [X(:), Y(:)];

hh = h/7;  % Spacement of the refined grid
xxRefined = -1:hh:1;
[XX, YY] = meshgrid(xxRefined); % Equal refinement on both directions

%% Calculating the cardinal functions
[M, uInterp, vInterp] = FD_DivFreeMatrixStream(dSites, N);

L = zeros(n,2*n);  % used to construct delta function

% coeffs = zeros(size(M,2), length(uIdx) + length(vIdx));
    % lines  : number of columns of the divFree matrix M
    % columns: total amount of interpolation points for u and v

lebesgueFunctionU = zeros(size(XX));
lebesgueFunctionV = zeros(size(XX));

rhs = zeros(2*n^2);
for p = 1:2*n^2 % all interpolation points
    L(p) = 1;
    U = L(:,1:n);
    V = L(:,n+1:2*n);
    rhs(:,p) = [U(:); V(:)];
    L(p) = 0; % going back to zero matrix
end

coeffs = M\rhs;

for p = 1:2*n^2 % 1:2*numPts^2 minus some pts.
    cardFunctionU = uInterp(XX,YY,coeffs(:,p));
    cardFunctionV = vInterp(XX,YY,coeffs(:,p));

    lebesgueFunctionU = lebesgueFunctionU + abs(cardFunctionU);
    lebesgueFunctionV = lebesgueFunctionV + abs(cardFunctionV);

    [i,j] = ind2sub([n, n], rem(p-1,n^2) + 1);
    % p = 1, ..., numPts^2, 1, ..., numPts^2
    if display
        figure(1)
        set(gcf, 'Position', [100,100, 600*2, 600])
        h1 = subplot(1,2,1);
        mesh(XX,YY,cardFunctionU);
        axis([xxRefined(1) xxRefined(end) xxRefined(1) xxRefined(end) ...
              -.3 1])
        title(['Cardinal Function U for (i,j) = (' num2str(i) ',' ...
               num2str(j) ')'])
    
        h2 = subplot(1,2,2);
        mesh(XX,YY,cardFunctionV);
        axis([xxRefined(1) xxRefined(end) xxRefined(1) xxRefined(end) ...
              -.3 1])
        title(['Cardinal Function V for (i,j) = (' num2str(i) ',' ...
               num2str(j) ')'])

        set([h1 h2], 'clim', [-.3 1])
        snapnow
        pause(0.01)
    end
end

%%
% Observe on the plots (case |n = 3| and |N = 3|) that the cardinal
% functions do not satisfy the interpolation condition exactly. This
% happens because we are looking for the least-squares solution when we use
% the backslash operator |"\"|.
%
% Here we are using the interpolants based on the stream function, then, by
% construction, they form a divergence-free vector field. However, as
% mentioned on <FD_DivFreeMatrixStream.html |FD_DivFreeMatrixStream.m|> the
% interpolation condition might not hold. The interpolant that we find is a
% least-square solution, so it gives a pair of divergence-free polynomials
% of certain degree that best interpolates the vector field.
%
% The rank deficient warning happens because we have repeated information
% for the columns of the interpolation matrix $M$, i.e., the base we are
% using to represent our function is redundant. Matlab returns the
% *least-squares solution with as few non-zero entries as possible*.

%% Plotting the Lebesgue functions and the Lebesgue constant
lebesgueConstU = max(max(lebesgueFunctionU));
lebesgueConstV = max(max(lebesgueFunctionV));
lebesgueConst = max([lebesgueConstU, lebesgueConstV]);

if display
    figure(2)
    set(gcf, 'Position', [100,100, 600*2, 600])
    subplot(1,2,1)
    mesh(XX,YY,lebesgueFunctionU);
    title(['Lebesgue Function U for (i,j) = (' num2str(i) ',' num2str(j)...
        '), \Lambda_U = ' num2str(lebesgueConstU)])
    
    subplot(1,2,2)
    mesh(XX,YY,lebesgueFunctionV);
    title(['Lebesgue Function V for (i,j) = (' num2str(i) ',' num2str(j)...
        '), \Lambda_V = ' num2str(lebesgueConstV)])
end

end