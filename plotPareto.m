function plotPareto(matrix)

V = scatteredInterpolant(matrix(:,1), matrix(:,2), matrix(:,3));

maxjes = ceil(100*max(matrix))/100;
minjes = floor(100*min(matrix))/100;

xx = minjes(1):0.01:maxjes(1);
yy = minjes(2):0.01:maxjes(2);
[xq,yq] = meshgrid(xx,yy);
zq = V(xq,yq);

figure;
mesh(xq,yq,zq)
xlabel('Track')
ylabel('Effort')
zlabel('Moment')
hold on
scatter3(matrix(:,1), matrix(:,2), matrix(:,3))