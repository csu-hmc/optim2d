% script for running the optim.m program
clear all
clear global
clear mex

% define which movement we are working on
movement = 'Winter/Winter_normal';

% general settings
problem.Solver = 'IPOPT';
problem.checkderivatives = 0;
problem.debug = 0;
problem.MaxIterations = 5000;
problem.ConstraintTol = .0001;
problem.Tol = .00001;
problem.Printinterval = 3.0;
problem.N = 60;					% number of collocation nodes in the movement cycle

% define an able-bodied model and the target gait data for the simulation
ablemodel.parameterfile 	= 'gait2d_par.xls';
ablemodel.type 				= 'able';
ablemodel.datafile			= movement;
ablemodel.Wtrack 	= 1;	% weight of tracking term in optimization objective
ablemodel.Weffort 	= 20;	% weight of muscle effort term in optimization objective
ablemodel.effort.fatigue 	= 0;
ablemodel.effort.Fmaxweighted = 0;
ablemodel.effort.exponent	= 3;
ablemodel.Wvalve 	= 0.001;	% weight of valve operating cost in optimization objective
ablemodel.discretization = 'euler';
ablemodel.reducedW = 0;

% define a below knee amputee model, based on the able-bodied model
bkamodel		= ablemodel;
bkamodel.type	= 'bka';
% 450 Nm/rad foot stiffness in able bodied subjects (Hansen et al., J Biomech 37: 1467–1474, 2004)
bkamodel.anklestiffness = 450;		% stiffness (Nm/rad) of prosthetic ankle

% define an above knee amputee model with C-Leg, based on the BKA model
clegmodel		= bkamodel;
clegmodel.type	= 'aka';
clegmodel.RA	= 7;			% knee actuator volume, use zero to disconnect the actuator from the knee
clegmodel.C1max = 0;			% valve constant for valve 1, cm^3 s^-1 MPa^-0.5 (use zero to keep valve closed)
clegmodel.C2max = 50;			% valve constant for valve 2, cm^3 s^-1 MPa^-0.5
clegmodel.B1	= 0.01;			% valve drag, MPa s cm^-3
clegmodel.B2 	= 0.01;			% same for valve 2
clegmodel.L		= 0;			% actuator leakage parameter, cm^3 s^-1 (Nm)^-1
clegmodel.k		= 3.0;			% accumulator stiffness, MPa cm^-3

% define a model with CCF knee, based on C-leg model
ccfmodel 		= clegmodel;
ccfmodel.C1max	= 50;			% now we allow valve 1 to open also

% solve the gait for the able-bodied model
problem.model 			= ablemodel;
problem.initialguess 	= 'mid';
problem.resultfile		= [movement '_result_able.mat'];
optim(problem);
disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
pause;

% solve the gait for the BKA model
disp('Starting optimal control solution process for BKA model...');
problem.model 			= bkamodel;
problem.initialguess 	= problem.resultfile;
problem.resultfile 		= [movement '_result_bka.mat'];
optim(problem);
disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
pause;

% solve the gait for the C-Leg model
disp('Starting optimal control solution process for C-LEG model...');
problem.model 			= clegmodel;
problem.initialguess 	= problem.resultfile;
problem.resultfile 		= [movement '_result_cleg.mat'];
optim(problem);
disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
pause;

% solve the gait for the CCF knee model
disp('Starting optimal control solution process for CCF model...');
problem.model 			= ccfmodel;
problem.initialguess 	= problem.resultfile;
problem.resultfile 		= [movement '_result_ccf.mat'];
optim(problem);

disp('All optimizations completed.')

