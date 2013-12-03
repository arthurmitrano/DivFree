%% RBFdivFreeInterpComplex
% Calculates the divergence-free interpolant for a matrix or vector of
% shape parameters |ep| for the evaluation points |ePoints|.
function interp = RBFdivFreeInterpComplex(Aep, d, r, d1, d2, F, G, ep)
% Calculates the divergence-free interpolant for a matrix or vector of
% shape parameters |ep| for the evaluation points |ePoints|.
%
% INPUT:
% Aep  : interpolation matrix depending on the shape parameter
% d    : right hand side vector
% r    : distance matrix
% d1   : difference matrix for 1st coordinate
% d2   : difference matrix for 2nd coordinate
% F, G : functions used to create the kernel and evaluate the 
%        interpolant (annonymous functions)
% ep   : vector or matrix of shape parameters
%
% OUTPUT:
% interp : struct containing the components (u and v) of the interpolant.
%          u and v are cell arrays of the same size of ep, containing the
%          evaluations of the interpolant for each evaluation point.
interp = struct('u', {{}}, 'v', {{}});
interp.u = cell(size(ep));
interp.v = cell(size(ep));

for i = 1:size(ep,1)
    for j = 1:size(ep,2)
%         tic
        coeffs = Aep(ep(i,j))\d;
%         disp('rank')
%         rank(Aep(ep(i,j)))
%         disp('cond')
%         cond(Aep(ep(i,j)))
        coeffs = reshape(coeffs, 2, numel(d)/2).';
        t = RBFdivFreeInterp(coeffs, r, d1, d2, F, G, ep(i,j));
        interp.u{i,j} = t(:,1);
        interp.v{i,j} = t(:,2);
%         disp(['i = ' num2str(i) ', j = ' num2str(j)])
%         toc
    end
end
end