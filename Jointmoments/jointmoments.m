function [M,dMdx] = jointmoments(x, model)
% Returns joint moments and derivatives of knee in residual leg

moms = model.gait2d('Jointmoments',x(1:50));
M = moms(2);
dx = 1e-6;

dMdx = spalloc(1,54, 34); %spalloc ??
% Finite differences
for i = 1:34
    xuse = x;
    xuse(i) = x(i)+dx;
    Mdiff = model.gait2d('Jointmoments',xuse(1:50));
    dMdxi = (Mdiff-moms)/dx;
    dMdx(i) = -dMdxi(2); %Only use knee moment
end

dMdx;
