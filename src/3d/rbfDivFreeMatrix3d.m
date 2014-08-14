%% rbfDivFreeMatrix3d
% Construct the interpolation matrix *A* for the divergence-free RBF
% interpolant and also returns matrices *Ax*, *Ay* and *Az* used to
% calculate the differentiation matrices (*Dx*, *Dy* and *Dz*).
%
%  INPUT:
%  r    : distance matrix
%  d1   : difference matrix on the first coordinate  (x).
%  d2   : difference matrix on the second coordinate (y).
%  d3   : difference matrix on the third coordinate  (z).
%  t    : shape parameter of the kernel.
%
%  OUTPUT:
%  A          : interpolation matrix. Invert to find interpolant coeffs.
%  Ax, Ay, Az : used to calculate derivative matrices, eg, Dx = Ax*inv(A).

%%
function [A, Ax, Ay, Az] = rbfDivFreeMatrix3d(r, d1, d2, d3, t)

%% Creating temp variables to save time
temp1 = t*r.^2;
temp2 = exp(-temp1);
temp3_1 = temp1 - t*d1.^2 - 1;
temp3_2 = temp1 - t*d2.^2 - 1;
temp3_3 = temp1 - t*d3.^2 - 1;
temp4 = -8*t^3 * d1 .* d2 .* d3 .* temp2;
temp5_1 =  2*t * d1.^2 - 1;
temp5_2 =  2*t * d2.^2 - 1;
temp5_3 =  2*t * d3.^2 - 1;
temp6_1 = 8*t^2 * d1 .* temp2;
temp6_2 = 8*t^2 * d2 .* temp2;
temp6_3 = 8*t^2 * d3 .* temp2;

%% Creating matrices entries
A11 = -4*t * temp2 .* temp3_1;
A22 = -4*t * temp2 .* temp3_2;
A33 = -4*t * temp2 .* temp3_3;
A12 = 4*t^2 * d1 .* d2 .* temp2;
A13 = 4*t^2 * d1 .* d3 .* temp2;
A23 = 4*t^2 * d2 .* d3 .* temp2;
% A21 = A12; A31 = A13; A32 = A23;

Ax11 = temp6_1 .* temp3_1;
Ax22 = temp6_1 .* (temp3_2 - 1);
Ax33 = temp6_1 .* (temp3_3 - 1);
Ax12 = -4*t^2 * d2 .* temp2 .* temp5_1;
Ax13 = -4*t^2 * d3 .* temp2 .* temp5_1;
Ax23 = temp4;
% Ax21 = Ax12; Ax31 = Ax13; Ax32 = Ax23;

Ay11 = temp6_2 .* (temp3_1 - 1);
Ay22 = temp6_2 .* temp3_2;
Ay33 = temp6_2 .* (temp3_3 - 1);
Ay12 = -4*t^2 * d1 .* temp2 .* temp5_2;
Ay13 = temp4;
Ay23 = -4*t^2 * d3 .* temp2 .* temp5_2;
% Ay21 = Ay12; Ay31 = Ay13; Ay32 = Ay23;

Az11 = temp6_3 .* (temp3_1 - 1);
Az22 = temp6_3 .* (temp3_2 - 1);
Az33 = temp6_3 .* temp3_3;
Az12 = temp4;
Az13 = -4*t^2 * d1 .* temp2 .* temp5_3;
Az23 = -4*t^2 * d2 .* temp2 .* temp5_3;
% Az21 = Az12; Az31 = Az13; Az32 = Az23;

%% Forming matrices
A = [A11 A12 A13; ...
     A12 A22 A23; ...
     A13 A23 A33];
 
Ax = [Ax11 Ax12 Ax13; ...
      Ax12 Ax22 Ax23; ...
      Ax13 Ax23 Ax33];
  
Ay = [Ay11 Ay12 Ay13; ...
      Ay12 Ay22 Ay23; ...
      Ay13 Ay23 Ay33];
  
Az = [Az11 Az12 Az13; ...
      Az12 Az22 Az23; ...
      Az13 Az23 Az33];