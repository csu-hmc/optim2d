function [F, g] = easyobj(x)

% f1 = sum(x.^2);
% g(:,1) = 2*transpose(x);
% 
% [f, dfdx, dfdxdot, dfdu] = easydyn1(x);
% 
% F = [f1;transpose(f);diag(dfdx)];
% g = [transpose(g);[dfdx,dfdxdot]]

global problem

[F,g] = problem.easyobj(x);
