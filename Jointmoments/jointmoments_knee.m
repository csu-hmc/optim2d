function [M,dMdx] = jointmoments_knee(x, model)
% Returns joint moments and derivatives of knee of both legs

moms = model.gait2d('Jointmoments',x(1:50));
M = moms([2,5]);
dx = 1e-6;

dMdx = spalloc(2,54, 68);
% Finite differences
for i = 1:34
    xuse = x;
    xuse(i) = x(i)+dx;
    Mdiff = model.gait2d('Jointmoments',xuse(1:50));
    dMdxi = (Mdiff-moms)/dx;
    dMdx(:,i) = -dMdxi([2,5])'; %Only use knee moment
end

dMdx;
