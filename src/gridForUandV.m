% Plot grid points for u and v

h = 0.05;
xu = -1:2*h:1;
yu = -1:h:1;
xv = -1:h:1;
yv = -1:2*h:1;
[XU, YU] = meshgrid(xu,yu);
[XV, YV] = meshgrid(xv,yv);
plot(XU,YU,'r.', XV,YV,'bo')
title('u dots and v circles')
shg