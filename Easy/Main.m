function [X,F,INFO] = Main()
% clear all
% Simple snopt application
global problem
x0 = 100;

% Variables
problem.N = 10; % Nodes
problem.ndof = 5;
problem.nstates = 2*problem.ndof;
problem.ncontrols = 3;
problem.nvarpernode = problem.nstates + problem.ncontrols;		% number of unknowns per time node
problem.nvar = problem.nvarpernode * problem.N;		%Number of variables (controls and states per time node)
problem.ncon = problem.nstates * problem.N;

X0 = x0*randn(problem.nvar,1);
L = -inf*ones(problem.nvar,1);
U = inf*ones(problem.nvar,1);

FL = [-inf;zeros(problem.ncon,1)];
FU = [inf;zeros(problem.ncon,1)];
userfun = 'easyobj';

problem.easyobj = @easyobj;
% Output informative files
snprint   ('probName.out');
snsummary ('prName.sum');
[X,F,INFO] = snopt(X0,L,U,FL,FU,userfun);
snprint    off;
snsummary off;

end

function [F,g]= easyobj(x)

global problem
N = problem.N;

h = 0.01;

g = zeros(1+problem.ncon,problem.nvar);

F(1) = sum(x(problem.ncon+1:end).^2);
g(1,:) = 2*[zeros(1,problem.ncon), transpose(x(problem.nstates*N+1:end))];

for i = 1:N
    if i < N
        x1 = [x(problem.ndof*(i-1)+(1:problem.ndof));x(problem.ndof*(N+(i-1))+(1:problem.ndof))];
        x2 = [x(problem.ndof*i+(1:problem.ndof));x(problem.ndof*(N+i)+(1:problem.ndof))];
    else
        x1 = [x(problem.ndof*(i-1)+(1:problem.ndof));x(problem.ndof*(N+(i-1))+(1:problem.ndof))];
        x2 = [x(1:problem.ndof);x(problem.ndof*N+(1:problem.ndof))];
    end
    [f, dfdx, dfdxdot, dfdu] = easydyn1(x2,(x1-x2)/h); %(x1+x2)/2
    
    F = [F;f];
    
    g(problem.nstates*(i-1)+1+(1:problem.nstates),(i-1)+(1:problem.nstates)) = -dfdxdot/h;
    g(problem.nstates*(i-1)+1+(1:problem.nstates),i+(1:problem.nstates)) = dfdx + dfdxdot/h;
    g(problem.nstates*(i-1)+1+(1:problem.nstates),i+problem.nstates+1) = dfdu;
end
end
