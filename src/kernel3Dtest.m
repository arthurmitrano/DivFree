%% Test if kernel matrix of the 3D case has columns with divergence zero
syms t x1 x2 x3
x = [x1; x2; x3];
h = ((4*t - 4*t^2*(x.'*x))*eye(3) + 4*t^2*(x*x.'))*exp(-t*(x.'*x));
pretty(h)
fprintf('\n')

%% Calculating the divergence of the first column (h1)
div_h1 = diff(h(1,1),x1) + diff(h(2,1),x2) + diff(h(3,1),x3);
disp('div_h1 = ')
pretty(simplify(div_h1))
fprintf('\n')

%% Calculating the divergence of the first column (h2)
div_h2 = diff(h(1,2),x1) + diff(h(2,2),x2) + diff(h(3,2),x3);
disp('div_h2 = ')
pretty(simplify(div_h2))
fprintf('\n')

%% Calculating the divergence of the first column (h3)
div_h3 = diff(h(1,3),x1) + diff(h(2,3),x2) + diff(h(3,3),x3);
disp('div_h3 = ')
pretty(simplify(div_h3))
fprintf('\n')