%% Lebesgue functions of divergence-free polynomial interpolants
% The purpose of this script is to visualize the cardinal and Lebesgue
% functions of the *local and global polynomial* divergence-free method
% based on the *stream function* idea. The user can change the
% interpolation points on the files <../uInterpPts.txt |uInterpPts.txt|>
% and <../vInterpPts.txt |vInterpPts.txt|> (for the local case only), as
% well as the degree $N$ of the polynomial
% interpolant $(u_{interp},v_{interp})$. For more information on how this
% polynomial is constructed see <FD_DivFreeMatrixStream.html
% FD_divFreeMatrixStream.html>.

%%
function [lebesgueConstU, lebesgueConstV] = ...
                                      lebesgueFunctions(n,N,numPts,display)
% Calculates Lebesgue function and it's constant. If numPts is not
% provided, numPts = n. If numPts present then display = true else false.
%
% n       : Number of grid points in one dimension
% N       : Degree of the bivariate interpolation polynomial
% numPts  : The main stecil will have numPts^2 points (optional)
% display : Display the cardinal function plots (optinal)
%
% NOTE: n and numPts must be odd numbers.

%% Setting up the function:
if (nargin <= 2)
    numPts = n;             % Using all the grid points to interpolate
    local = false;
    display = false;
else
    local = (n ~= numPts);  % The main stecil will have numPts^2 points
end

if (nargin == 3)
    display = true;
end
%% Generating 2d-grid of n^2 pts equally-spaced
xx = linspace(-1,1,n);
h = abs(xx(2) - xx(1));

hh = h/7;  % Spacement of the refined grid
if local
    xxRefined = -floor(numPts/2)*h:hh:floor(numPts/2)*h;
else
    xxRefined = -1:hh:1;
end
[XX, YY] = meshgrid(xxRefined); % Equal refinement on both directions

%% Selecting the interpolation points
% If local, we center the stencil at the origin. To choose specific points
% of interpolations use the files <../uInterpPts.txt |uInterpPts.txt|>
% and <../vInterpPts.txt |vInterpPts.txt|>. Make sure that the matrix is on
% the correct dimension for a particular numPts _(numPts x numPts)_.
if local
    uIdx = find(load('uInterpPts.txt'));
    vIdx = find(load('vInterpPts.txt')); % using linear indexes
else % global, use all points of the grid
    uIdx = find(ones(numPts));
    vIdx = find(ones(numPts));
end
interpPts = struct('u', uIdx,'v', vIdx, 'numPts',numPts);

%% Calculating the cardinal functions
[M, uInterp, vInterp] = FD_DivFreeMatrixStream(h,N,interpPts);

L = zeros(numPts,2*numPts);  % used to construct delta function
coeffs = zeros(size(M,2), length(uIdx) + length(vIdx));
    % lines  : number of columns of the divFree matrix M
    % columns: total amount of interpolation points for u and v

lebesgueFunctionU = zeros(size(XX));
lebesgueFunctionV = zeros(size(XX));
for p = [uIdx; vIdx+numPts^2]' % 1:2*numPts^2 minus some pts.
    L(p) = 1;
    U = L(:,1:numPts);
    V = L(:,numPts+1:2*numPts);
    L(p) = 0; % going back to zero matrix
    coeffs(:,p) = M\[U(uIdx); V(vIdx)];  % WHY DO I KEEP THOSE COEFFS
    [i,j] = ind2sub([numPts, numPts], rem(p-1,numPts^2)+1); 
                                   % p = 1, ..., numPts^2, 1, ..., numPts^2
    cardFunctionU = uInterp(XX,YY,coeffs(:,p));
    cardFunctionV = vInterp(XX,YY,coeffs(:,p));
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
        pause
    end

    lebesgueFunctionU = lebesgueFunctionU + abs(cardFunctionU);
    lebesgueFunctionV = lebesgueFunctionV + abs(cardFunctionV);
end

%%
% Observe on the plots (case |n = 51|, |N = 3| and |numPts = 3|) that the
% cardinal functions do not satisfy the interpolation condition exactly.
% This happens because we are looking for the least-squares solution when
% we use the backslash operator |"\"|.
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

if display
    figure(2)
    set(gcf, 'Position', [100,100, 600*2, 600])
    subplot(1,2,1)
    mesh(XX,YY,lebesgueFunctionU);
    title(['Lebesgue Function U for (i,j) = (' num2str(i) ',' num2str(j) ...
        '), \Lambda_U = ' num2str(lebesgueConstU)])
    
    subplot(1,2,2)
    mesh(XX,YY,lebesgueFunctionV);
    title(['Lebesgue Function V for (i,j) = (' num2str(i) ',' num2str(j) ...
        '), \Lambda_V = ' num2str(lebesgueConstV)])
end
end