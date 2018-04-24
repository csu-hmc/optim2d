function [f, dfdx, dfdxdot, dfdu] = easydyn(x, xdot)

f = x+xdot;
dfdx = speye(54);
dfdxdot = speye(54);
dfdu = spalloc(54,18,0);