%% chebpoly : Returns symbolic expression of Chebyshev polynomial
function T = chebpoly(n)
syms x
if n == 0
    T = sym(1);
elseif n == 1
    T = x;
elseif mod(n,2) == 0
    T = 2*chebpoly(n/2)^2 - 1;
else
    T = 2*chebpoly((n-1)/2)*chebpoly((n+1)/2) - x;
end