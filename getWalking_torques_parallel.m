clear all
close all
clc

% define which movement we are working on
movement = 'Winter/Winter_normal';

% general settings
problem.Solver = 'IPOPT';
problem.MaxIterations = 5000;
problem.ConstraintTol = .0001;
problem.Tol = .0001;
problem.symmetry = 1;
problem.discretization = 'BE';		% start with backward Euler
problem.checkderivatives = 0;		% set this to 1 to make optim.m check derivatives
problem.debug = 0;                  % only used in debug mode
problem.N = 60;
if problem.symmetry
    problem.N = problem.N/2;
end
problem.Ncycles = 1;
problem.warmstart = 0;
problem.Printinterval = 1;
problem.proscontrol = 0;
   
% define an able-bodied model and the target gait data for the simulation
ablemodel.parameterfile 	= 'gait2d_par.xls';
ablemodel.type 				= 'able';
ablemodel.datafile			= movement;
ablemodel.Wtrack 	= 10;	% weight of tracking term in optimization objective
ablemodel.Weffort 	= 0.01;	% weight of muscle effort term in optimization objective
ablemodel.Wreg = 0;       % weight of the regularization term in optimization objective
ablemodel.effort.fatigue 	= 0;
ablemodel.effort.Fmaxweighted = 0;
ablemodel.effort.exponent	= 2;
ablemodel.Wvalve 	= 0.00;	% weight of valve operating cost in optimization objective
ablemodel.discretization = 'euler';
ablemodel.reducedW = 0;

torquemodel		= ablemodel;
torquemodel.type= 'torque';

Ncycles = 5;
for i = 1:5
    for j = 1:length(Ncycles)
        result = runoptimization(problem, torquemodel, movement, Ncycles(j),i);
    end
end