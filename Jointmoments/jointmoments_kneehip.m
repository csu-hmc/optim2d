function [M,dMdx] = jointmoments_kneehip(x, model)
% Returns joint moments and derivatives of knee and hip for both legs

moms = model.gait2d('Jointmoments',x(1:50));
M = moms([1,2,4,5]);
dx = 1e-7;

dMdx = spalloc(4,50, 68); %spalloc ??
% Finite differences
for i = 1:34
    xuse = x;
    xuse(i) = x(i)+dx;
    Mdiff = model.gait2d('Jointmoments',xuse(1:50));
    dMdxi = (Mdiff-moms)/dx;
    dMdx(:,i) = dMdxi([1,2,4,5])'; %Only use knee moment
end

dMdx;
