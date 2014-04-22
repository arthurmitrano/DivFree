%% f = chebpoly2function
% f = chebpoly2function(A) returns an anonymous function using the
% coefficients given by matrix A where any number less than 10^{-12} is
% considered 0.
%
% f = chebpoly2function(A,tol) do the same task but using a specific
% tolerance value tol.

%%
function f = chebpoly2function(A,tol)
% f = chebpoly2function(A) returns an anonymous function using the
% coefficients given by matrix A where any number less than 10^{-12} is
% considered 0. f = sum_i ( sum_j Y(i,j) T_i(y) T_j(x) ), where
% Y=rot90(A,2)
%
% f = chebpoly2function(A,tol) do the same task but using a specific
% tolerance value tol.

if nargin < 2
    tol = 1e-12;
end

trunc = (abs(A)>tol);
Y = rot90(real(A.*trunc),2);

J = find(sum(abs(Y)) ~= 0);
I = find(sum(abs(Y),2)' ~= 0);

% function to create the Chebyshev polynomial T_i anonymous
if isempty(I) && isempty(J)
    N = 0;
elseif ~isempty(I) && ~isempty(J)
    N = max([I(end) J(end)]);
elseif ~isempty(I)
    N = I(end);
else
    N = J(end);
end
     
T = cell(N,1);
for k = 0:N-1
    T{k+1} = matlabFunction(chebpoly(k), 'vars', 'x');
end

syms x y
f = 0;
for i = I
    for j = J
        f = f + Y(i,j).*T{i}(y).*T{j}(x);
    end
end
f = matlabFunction(f, 'vars', [x y]);