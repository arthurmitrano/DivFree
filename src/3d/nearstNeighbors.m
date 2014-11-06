%% nearstNeighbors
% Find the $n$ nearst neighbors of a given point.
%
%  INPUT:
%  P : matrix with all points of the domain.
%  p : point that we want to know the nearst neighbors.
%  n : amount of nearst neighbors that we want to find of point p.
%
%  OUTPUT:
%  neighbors : n nearst neighbors of point p.

%%
function neighbors = nearstNeighbors(P, p, n)

%% Setting up the function
d = zeros(size(P,1),1);

%% Calculating the nearst neighbors
for k = 1:size(P,1)
    d(k) = norm(P(k,:) - p(:)',2);
end
[~, i] = sort(d);
neighbors = P(i(1:n),:);
