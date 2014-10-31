%% testFunction3D
% This function generates a divergence-free vector field given some vector
% valued function $F(x,y)$.
%
%  INPUT:
%  P  : list of points where the test function will be evaluated.
%  p  : point where we calculate the derivatives of the test function.
%  F1 : creates divFree vector field (default: sin(x^2 + y^2 + z^2).
%  F2 : creates divFree vector field (default: cos(x^2 + y^2 + z^2).
%  F3 : creates divFree vector field (default: tanh(x^2 + y^2 + z^2).
%
%  OUTPUT:
%  u, v, w    : divergence vector field at the points (X,Y,Z).
%  Du, Dv, Dw : struct containing the derivatives at the grid center.
%%
function [u, v, w, Du, Dv, Dw] = testFunction3d(P,p,F1,F2,F3)

%% Setting up the function
if nargin < 3
    F1 = @(x,y,z) sin(x.^2 + y.^2 + z.^2);
    F2 = @(x,y,z) cos(x.^2 + y.^2 + z.^2);
    F3 = @(x,y,z) tanh(x.^2 + y.^2 + z.^2);
end

Du = struct('x',{[]}, 'y',{[]}, 'z',{[]});
Dv = struct('x',{[]}, 'y',{[]}, 'z',{[]});
Dw = struct('x',{[]}, 'y',{[]}, 'z',{[]});

syms x y z

%% Defining a divergence-free vector field
u = diff(F3(x,y,z),y) - diff(F2(x,y,z),z); 
v = diff(F1(x,y,z),z) - diff(F3(x,y,z),x);
w = diff(F2(x,y,z),x) - diff(F1(x,y,z),y);

%% Calculating derivatives
Du.x = matlabFunction(diff(u,x), 'vars',[x,y,z]);
Du.y = matlabFunction(diff(u,y), 'vars',[x,y,z]);
Du.z = matlabFunction(diff(u,z), 'vars',[x,y,z]);

Dv.x = matlabFunction(diff(v,x), 'vars',[x,y,z]);
Dv.y = matlabFunction(diff(v,y), 'vars',[x,y,z]);
Dv.z = matlabFunction(diff(v,z), 'vars',[x,y,z]);

Dw.x = matlabFunction(diff(w,x), 'vars',[x,y,z]);
Dw.y = matlabFunction(diff(w,y), 'vars',[x,y,z]);
Dw.z = matlabFunction(diff(w,z), 'vars',[x,y,z]);

%% Evaluating derivatives at p
Du.x = Du.x(p(1),p(2),p(3));
Du.y = Du.y(p(1),p(2),p(3));
Du.z = Du.z(p(1),p(2),p(3));

Dv.x = Dv.x(p(1),p(2),p(3));
Dv.y = Dv.y(p(1),p(2),p(3));
Dv.z = Dv.z(p(1),p(2),p(3));

Dw.x = Dw.x(p(1),p(2),p(3));
Dw.y = Dw.y(p(1),p(2),p(3));
Dw.z = Dw.z(p(1),p(2),p(3));

%% Calculating the divergece-free vector field at grid points P
u = matlabFunction(u, 'vars',[x,y,z]);
v = matlabFunction(v, 'vars',[x,y,z]);
w = matlabFunction(w, 'vars',[x,y,z]);

u = u(P(:,1),P(:,2),P(:,3));
v = v(P(:,1),P(:,2),P(:,3));
w = w(P(:,1),P(:,2),P(:,3));