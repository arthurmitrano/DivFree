%% Symbolic flat limit of the divergence-free RBF interpolant
% The objective of this script is to look to the symbolic expressions of
% the divergence-free RBF and see what happens when $\varepsilon
% \rightarrow 0$. We are looking for:
%
% $$ \lim_{\varepsilon \rightarrow 0} t^\varepsilon(\mathbf{x}) =
% \sum_{k=1}^N \lim_{\varepsilon \rightarrow 0}
% \Psi^\varepsilon(\mathbf{x},\mathbf{x}_k)\mathbf{s}^{\varepsilon}_k.$$
%
% Where
%
% $$\Psi^\varepsilon(\mathbf{x},\mathbf{y}) = -F^\varepsilon(r)I -
% G^\varepsilon(r)(\mathbf{x} - \mathbf{y})(\mathbf{x} - \mathbf{y})^T,$$
%
% with
%
% $$ F^\varepsilon(r) = \frac{\varphi^\prime(r,\varepsilon)}{r} \quad and
% \quad G^\varepsilon(r) = \frac{{F^\varepsilon}^\prime(r)}{r}.$$
%
% Here $r$ is the distance between $\mathbf{x}$ and $\mathbf{y}$ and
% $\varphi$ the basic function (gaussian, multiquadrics, inverse
% multiquadrics, etc).
%
% On the end, what really matters is:
% 
% $$ \lim_{\varepsilon \rightarrow 0}
% F^\varepsilon(r)\mathbf{s}^{\varepsilon}_k \quad and \quad
% \lim_{\varepsilon \rightarrow 0} G^\varepsilon(r)(\mathbf{x} -
% \mathbf{y})(\mathbf{x} - \mathbf{y})^T\mathbf{s}^{\varepsilon}_k.$$

%% Setup for the script
clear, clc
syms ep r x1 x2 y1 y2

x = [x1; x2];
y = [y1; y2];

x_e = [1/3;1/5];  % evaluation point

% rbf = @(ep,r) exp(-(ep*r).^2);
rbf = @(ep,r) 1./(1 + (ep*r).^2);
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

F = matlabFunction(diff(rbf(ep,r),r)/r);
G = matlabFunction(diff(F,r)/r);
K = @(ep,x,y) -F(ep,norm(x-y))*eye(2) ...
              -G(ep,norm(x-y))*(x-y)*(x-y).';  % divergence-free kernel

%% 1 interpolation point U(0,0) = 1

disp('1 interpolation point')

U = [1];
V = [0];
d = [U; V];

A = K(ep,[0;0],[0;0]);
s = A\d;
t = K(ep,x_e,[0;0])*s;
limit(t,ep,0)

%% 2 interpolation points U(0,0) = 1 & U(1,0) = 0

disp('2 interpolation points')

U = [1 0];
V = [0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0]';  % interpolation points

A = [K(ep,x_i(:,1),x_i(:,1)) K(ep,x_i(:,1),x_i(:,2)); ...
     K(ep,x_i(:,2),x_i(:,1)) K(ep,x_i(:,2),x_i(:,2))];
s = A\d;
t = K(ep,x_e,x_i(:,1))*s(1:2) + K(ep,x_e,x_i(:,2))*s(3:4);
limit(t,ep,0)

%% 3 interpolation points U(0,0) = 1, U(1,0) = 0 & U(-1,0) = 0

disp('3 interpolation points')

U = [1 0 0];
V = [0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0; -1 0]';  % interpolation points

A = [K(ep,x_i(:,1),x_i(:,1)) K(ep,x_i(:,1),x_i(:,2)) K(ep,x_i(:,1),x_i(:,3)); ...
     K(ep,x_i(:,2),x_i(:,1)) K(ep,x_i(:,2),x_i(:,2)) K(ep,x_i(:,2),x_i(:,3)); ...
     K(ep,x_i(:,3),x_i(:,1)) K(ep,x_i(:,3),x_i(:,2)) K(ep,x_i(:,3),x_i(:,3))];
s = A\d;
t = K(ep,x_e,x_i(:,1))*s(1:2) + K(ep,x_e,x_i(:,2))*s(3:4) + K(ep,x_e,x_i(:,3))*s(5:6);
limit(t,ep,0)

%% 4 interpolation points U(0,0) = 1, U(1,0) = 0, U(-1,0) = 0 & U(1,1) = 0

disp('4 interpolation points')

U = [1 0 0 0];
V = [0 0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0; -1 0; 1 1]';  % interpolation points

A = [K(ep,x_i(:,1),x_i(:,1)) K(ep,x_i(:,1),x_i(:,2)) K(ep,x_i(:,1),x_i(:,3)) K(ep,x_i(:,1),x_i(:,4)); ...
     K(ep,x_i(:,2),x_i(:,1)) K(ep,x_i(:,2),x_i(:,2)) K(ep,x_i(:,2),x_i(:,3)) K(ep,x_i(:,2),x_i(:,4)); ...
     K(ep,x_i(:,3),x_i(:,1)) K(ep,x_i(:,3),x_i(:,2)) K(ep,x_i(:,3),x_i(:,3)) K(ep,x_i(:,3),x_i(:,4)); ...
     K(ep,x_i(:,4),x_i(:,1)) K(ep,x_i(:,4),x_i(:,2)) K(ep,x_i(:,4),x_i(:,3)) K(ep,x_i(:,4),x_i(:,4))];
s = A\d;
t = K(ep,x_e,x_i(:,1))*s(1:2) + K(ep,x_e,x_i(:,2))*s(3:4) + ...
    K(ep,x_e,x_i(:,3))*s(5:6) + K(ep,x_e,x_i(:,4))*s(7:8);
limit(t,ep,0)


%% 5 interpolation points U=1 @ (0,0), U = 0 at the other points

disp('5 interpolation points')

U = [1 0 0 0 0];
V = [0 0 0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 0 1; 0 -1; -1 0; 1 0]';  % interpolation points

for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end
limit(t,ep,0)


%% 9 interpolation points U=1 @ (1,1), U = 0 at the other points
tic
disp('9 interpolation points')

U = [0 0 0 1 0 0 0 0 0];
V = [0 0 0 0 0 0 0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1]';  % interpolation points

for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d; % Takes 700 seconds for inverse quadrics
t = [0;0];
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end
toc
tic
limit(t,ep,0)
toc