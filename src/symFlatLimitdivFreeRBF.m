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
% G^\varepsilon(r)\left[\matrix{(x_2 - y_2)^2&-(x_1 - y_1)(x_2 - y_2)\cr
% -(x_1 - y_1)(x_2 - y_2)&(x_1 - y_1)^2\cr}\right],$$
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
syms ep r x1 x2 y1 y2 X Y
assume(X,'real')
assume(Y,'real')
x = [x1; x2];
y = [y1; y2];

x_e = [X;Y];  % evaluation point

takeLimit = true;

rbf = @(ep,r) exp(-(ep*r).^2);
% rbf = @(ep,r) 1./(1 + (ep*r).^2)^2;
% rbf = @(ep,r) 1./sqrt(1 + (ep*r).^2);
% rbf = @(ep,r) sqrt(1 + (ep*r).^2);

% Constructing divergence-free kernel
F = matlabFunction(diff(rbf(ep,r),r)/r);
G = matlabFunction(diff(F(ep,r),r)/r);
R = sym([0 -1; 1 0]);
K = @(ep,x,y) -F(ep,norm(x-y))*eye(2) ...
              -G(ep,norm(x-y))*R*(x-y)*(x-y).'*R.';

%% 1 interpolation point U(0,0) = 1

disp('1 interpolation point')

U = [1];
V = [0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0]';  % interpolation points

A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

if takeLimit
    disp('taking limit')
    L = limit(simplify(t),ep,0);
    disp('Limit when ep -> 0')
    pretty(L)
end

%% 2 interpolation points U(0,0) = 1 & U(1,0) = 0

disp('2 interpolation points')

U = [1 0];
V = [0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0]';  % interpolation points

A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

if takeLimit
    disp('taking limit')
    L = limit(simplify(t),ep,0);
    disp('Limit when ep -> 0')
    pretty(L)
end

%% 3 interpolation points U(0,0) = 1, U(1,0) = 0 & U(-1,0) = 0

disp('3 interpolation points')

U = [1 0 0];
V = [0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 1 0; -1 0]';  % interpolation points

A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

if takeLimit
    disp('taking limit')
    L = limit(simplify(t),ep,0);
    disp('Limit when ep -> 0')
    pretty(L)
end

%% 4 interpolation points U(0,0) = 1, U(1,0) = 0, U(-1,0) = 0 & U(1,1) = 0

disp('4 interpolation points')

U = [1 0 0 0];
V = [0 0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [-1 -1; 1 1; -1 1; 1 -1]';  % interpolation points

A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

if takeLimit
    disp('taking limit')
    L = limit(simplify(t),ep,0);
    disp('Limit when ep -> 0')
    pretty(L)
end

%% 5 interpolation points U=1 @ (0,0), U = 0 at the other points

disp('5 interpolation points')

U = [1 0 0 0 0];
V = [0 0 0 0 0];
t = [U(:) V(:)];
d = reshape(t.',1,numel(t)).';

x_i = [0 0; 0 1; 0 -1; -1 0; 1 0]';  % interpolation points

A = sym('A',2*size(x_i,2));
for i = 1:size(x_i,2)
    for j = 1:size(x_i,2)
        A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
    end
end

s = A\d;
t = [0;0];  % initializing
for i = 1:size(x_i,2)
    t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
end

if takeLimit
    disp('taking limit')
    L = limit(simplify(t),ep,0);
    disp('Limit when ep -> 0')
    pretty(L)
end

% %% 9 interpolation points U=1 @ (1,1), U = 0 at the other points
% tic
% disp('9 interpolation points')
% 
% U = [0 0 0 1 0 0 0 0 0];
% V = [0 0 0 0 0 0 0 0 0];
% t = [U(:) V(:)];
% d = reshape(t.',1,numel(t)).';
% 
% x_i = [0 0; 1 0; -1 0; 1 1; 0 1; -1 1; 1 -1; 0 -1; -1 -1]';  % interpolation points
% 
% A = sym('A',2*size(x_i,2));
% for i = 1:size(x_i,2)
%     for j = 1:size(x_i,2)
%         A((i-1)*2+1:2*i, (j-1)*2+1:2*j) = K(ep,x_i(:,i),x_i(:,j));
%     end
% end
% s = A\d;    % takes 700 seconds
% 
% t = [0;0];  % initializing
% for i = 1:size(x_i,2)
%     t = t + K(ep,x_e,x_i(:,i))*s((i-1)*2+1:2*i);
% end
% toc
% 
% tic
% if takeLimit
%     disp('taking limit')
%     L = limit(simplify(t),ep,0);
%     disp('Limit when ep -> 0')
%     pretty(L)
% end
% toc
% 
% t1 = matlabFunction(t(1));
% t2 = matlabFunction(t(2));
% 
% figure(1)
% 
% subplot(1,2,1)
% mesh(xx,yy,t1(xx,yy,e))
% title('U'), xlabel('x'), ylabel('y')
% 
% subplot(1,2,2)
% mesh(xx,yy,t2(xx,yy,e))
% title('V'), xlabel('x'), ylabel('y')
% snapnow
