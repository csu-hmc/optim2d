function [f, dfdx, dfdxdot, dfdu] = easydyn1(x, xdot)

f = x+xdot;
dfdx = speye(length(x));
dfdxdot = speye(length(x));
dfdu = spalloc(length(x),1,0);