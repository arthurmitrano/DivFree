%% Lebesgue functions of divergence-free polynomial interpolants
clc, clear, close all
%% Setting up the FD_DivFreeMatrixStream function:
N = 3;           % Degree of the bivariate polynomial
numPts = 3;      % The main stecil will have numPts^2 points
k = 1;           % Order of the derivative

%% Generating 2d-grid n^2 pts equally-spaced
n = 50; % number of points in one dimension
xx = linspace(-1,1,n);  yy = linspace(-1,1,n);
[X, Y] = meshgrid(xx, yy);
h = abs(xx(2) - xx(1));
hh = h/20;
xxRefined = -1.0*h:hh:1.0*h;  yyRefined = xxRefined;
[XX, YY] = meshgrid(xxRefined, yyRefined);

%% Selecting the interpolation points
uIdx = find(load('uInterpPts.txt') == 1);
vIdx = find(load('vInterpPts.txt') == 1); % using linear indexes
interpPts = struct('u', uIdx,'v', vIdx, 'numPts',numPts);

%% Calculating the cardinal functions
[M, t1, t2] = FD_DivFreeMatrixStream(h,N,interpPts,k);

syms x y
L = zeros(numPts,2*numPts);
coeffs = zeros(size(M,2), length(uIdx) + length(vIdx));
lebesgueFunctionU = zeros(size(XX)); lebesgueFunctionV = zeros(size(XX));
for p = [uIdx; vIdx+numPts^2]' % 1:2*numPts^2 minus some pts.  
    L(p) = 1;
    U = L(:,1:numPts);  
    V = L(:,numPts+1:2*numPts);
    L(p) = 0; % going back to zero matrix
        
    coeffs(:,p) = M\[U(uIdx); V(vIdx)];
    
    cardFunctionU = matlabFunction(t1(x,y)*coeffs(:,p));
    cardFunctionV = matlabFunction(t2(x,y)*coeffs(:,p));
    
    [i,j] = ind2sub([numPts, numPts], rem(p-1,numPts^2)+1); 
                                   % p = 1, ..., numPts^2, 1, ..., numPts^2
    
    figure(1)
    set(gcf, 'Position', [100,100, 600*2, 600])
    subplot(1,2,1)
    mesh(XX,YY,cardFunctionU(XX,YY));
    title(['Cardinal Function U for (i,j) = (' num2str(i) ',' num2str(j)...
           ')'])
    
    subplot(1,2,2)
    mesh(XX,YY,cardFunctionV(XX,YY));
    title(['Cardinal Function V for (i,j) = (' num2str(i) ',' num2str(j)...
           ')'])
    
    snapnow
    pause(.1)
    
    lebesgueFunctionU = lebesgueFunctionU + abs(cardFunctionU(XX,YY));
    lebesgueFunctionV = lebesgueFunctionV + abs(cardFunctionV(XX,YY));
    
end

%% Plotting the Lebesgue functions and the Lebesgue constant

figure(2)
set(gcf, 'Position', [100,100, 600*2, 600])
subplot(1,2,1)
mesh(XX,YY,lebesgueFunctionU);
title(['Lebesgue Function U for (i,j) = (' num2str(i) ',' num2str(j) ...
       '), \Lambda_U = ' num2str(max(max(lebesgueFunctionU)))])

subplot(1,2,2)
mesh(XX,YY,lebesgueFunctionV);
title(['Cardinal Function V for (i,j) = (' num2str(i) ',' num2str(j) ...
       '), \Lambda_V = ' num2str(max(max(lebesgueFunctionV)))])