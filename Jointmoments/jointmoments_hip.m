function [M,dMdx] = jointmoments_hip(x, model)
% Returns joint moment and derivatives of hip of both legs 

moms = model.gait2d('Jointmoments',x(1:50));
M = moms([1,4]);
dx = 1e-6;

dMdx = spalloc(2,54, 68);

% Finite differences
for i = 1:34
    xuse = x;
    xuse(i) = x(i)+dx;
    Mdiff = model.gait2d('Jointmoments',xuse(1:50));
    dMdxi = (Mdiff-moms)/dx;
    dMdx(:,i) = -dMdxi([1,4])'; 
end

dMdx;
