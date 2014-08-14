%% testFunction3D
% This function generates a divergence-free vector field given some vector
% valued function $F(x,y)$.
%
%  INPUT:
%  X,Y,Z : matrices of coordinates where the divergence-free vector 
%          field will be evaluated.
%  F1    : creates divFree vector field (default: sin(x^2 + y^2 + z^2).
%  F2    : creates divFree vector field (default: cos(x^2 + y^2 + z^2).
%  F3    : creates divFree vector field (default: tanh(x^2 + y^2 + z^2).
%
%  OUTPUT:
%  u, v, w    : divergence vector field at the points (X,Y,Z).
%  Du, Dv, Dw : struct containing the derivatives at the grid center.
%  O          : struct containing the indexes of the center of the grid.
%%
function [u, v, w, Du, Dv, Dw, O] = testFunction3d(X,Y,Z,F1,F2,F3)

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
Du.x = Du.x(X,Y,Z).*ones(size(X));
Du.y = Du.y(X,Y,Z).*ones(size(X));
Du.z = Du.z(X,Y,Z).*ones(size(X));

Dv.x = matlabFunction(diff(v,x), 'vars',[x,y,z]);
Dv.y = matlabFunction(diff(v,y), 'vars',[x,y,z]);
Dv.z = matlabFunction(diff(v,z), 'vars',[x,y,z]);
Dv.x = Dv.x(X,Y,Z).*ones(size(X));
Dv.y = Dv.y(X,Y,Z).*ones(size(X));
Dv.z = Dv.z(X,Y,Z).*ones(size(X));

Dw.x = matlabFunction(diff(w,x), 'vars',[x,y,z]);
Dw.y = matlabFunction(diff(w,y), 'vars',[x,y,z]);
Dw.z = matlabFunction(diff(w,z), 'vars',[x,y,z]);
Dw.x = Dw.x(X,Y,Z).*ones(size(X));
Dw.y = Dw.y(X,Y,Z).*ones(size(X));
Dw.z = Dw.z(X,Y,Z).*ones(size(X));

%% Evaluating derivatives at the origin
gridCenterX = ceil(size(X,2)/2);
gridCenterY = ceil(size(Y,1)/2);
gridCenterZ = ceil(size(Z,1)/2);
O = struct('i',gridCenterX, 'j',gridCenterY, 'k',gridCenterZ);

Du.x = Du.x(O.i,O.j,O.k);
Du.y = Du.y(O.i,O.j,O.k);
Du.z = Du.z(O.i,O.j,O.k);


Dv.y = Dv.y(O.i,O.j,O.k);
Dv.z = Dv.z(O.i,O.j,O.k);

Dw.x = Dw.x(O.i,O.j,O.k);
Dw.y = Dw.y(O.i,O.j,O.k);
Dw.z = Dw.z(O.i,O.j,O.k);

%% Calculating the functions at grid points
u = matlabFunction(u, 'vars',[x,y,z]);
v = matlabFunction(v, 'vars',[x,y,z]);
w = matlabFunction(w, 'vars',[x,y,z]);
u = u(X,Y,Z); v = v(X,Y,Z); w = w(X,Y,Z);