% script for running the optim.m program
clear all

% tic
% define which movement we are working on
movement = 'Winter/Winter_slow';

% general settings
problem.Solver = 'IPOPT';
problem.MaxIterations = 5000;
problem.ConstraintTol = .01;
problem.Tol = .0001;
problem.symmetry = 1;
problem.discretization = 'BE';		% start with backward Euler
problem.checkderivatives = 0;		% set this to 1 to make optim.m check derivatives
problem.debug = 0;                  % only used in debug mode
problem.N = 60;						% start with a coarse mesh
problem.Printinterval = 3;
   
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

% Find standing model
% define model used for standing
standmodel = ablemodel;
standmodel.speed = 0;
standmodel.Wtrack = 1e-8;

% optimize for standing
problem.N = 2;
problem.model = standmodel;
problem.initialguess = 'mid';
problem.resultfile		= [movement '_result_stand_2.mat'];
result = optim(problem);

% optimize for standing
problem.N = 30;
problem.model = standmodel;
problem.initialguess = 'Winter/Winter_slow_result_stand_2';
problem.resultfile		= [movement '_result_stand_30.mat'];
result = optim(problem);

% optimize for standing
problem.N = 60;
problem.model = standmodel;
problem.initialguess = 'Winter/Winter_slow_result_stand_30';
problem.resultfile		= [movement '_result_stand.mat'];
result = optim(problem);

% define model used for slowest walking
slowestmodel = ablemodel;
slowestmodel.speed = 0.15;
slowestmodel.Wtrack = 0.05;

% optimize for slow walking
% problem.N = 60;
problem.model = slowestmodel;
problem.initialguess = 'Winter/Winter_slow_result_stand.mat';  % 'mid';%
problem.resultfile		= [movement '_result_slowest.mat'];
result = optim(problem);
% disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
% pause;

% define model used for slow walking
slowmodel = ablemodel;
slowmodel.speed = 0.5;
slowmodel.Wtrack = 0.1;

% optimize for slow walking
% problem.N = 60;
problem.model = slowmodel;
problem.initialguess = 'Winter/Winter_slow_result_slowest.mat';  % 'mid';%
problem.resultfile		= [movement '_result_slow.mat'];
result = optim(problem);
% disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
% pause;

% define model used for slow walking
slow1model = ablemodel;
slow1model.speed = 0.7;
slow1model.Wtrack = 0.3;

% optimize for slow walking
% problem.N = 60;
problem.model = slow1model;
problem.initialguess = 'Winter/Winter_slow_result_slow.mat';  % 'mid';%
problem.resultfile		= [movement '_result_slow1.mat'];
result = optim(problem);
% disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
% pause;

% define model used for slow walking
slow2model = ablemodel;
slow2model.speed = 0.9;
slow2model.Wtrack = 1;

% optimize for slow walking
% problem.N = 60;
problem.model = slow2model;
problem.initialguess = 'Winter/Winter_slow_result_slow1.mat';  % 'mid';%
problem.resultfile		= [movement '_result_slow2.mat'];
result = optim(problem);
% disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
% pause;


% Find result for walking
problem.model = ablemodel;
problem.initialguess = 'Winter/Winter_slow_result_slow2.mat';
problem.resultfile      = [movement '_result_able.mat'];
optim(problem)

% % define a below knee amputee model, based on the able-bodied model
% bkamodel		= ablemodel;
% bkamodel.type	= 'bka';
% % 450 Nm/rad foot stiffness in able bodied subjects (Hansen et al., J Biomech 37: 1467–1474, 2004)
% bkamodel.anklestiffness = 10000;		% stiffness (Nm/rad) of prosthetic ankle
% %  
% % solve the gait for the BKA model
% disp('Starting optimal control solution process for BKA model...');
% problem.model 			= bkamodel;
% problem.initialguess 	= 'Winter/Winter_slow_result_able';
% problem.resultfile 		= [movement '_result_bka_10000.mat'];
% optim(problem);
% 
% % define an above knee amputee model with zero knee moment
% akamodel		= bkamodel;
% akamodel.type	= 'aka';
% akamodel.RA	= 0;			% knee actuator volume, use zero to disconnect the actuator from the knee
% akamodel.C1max = 0;			% valve constant for valve 1, cm^3 s^-1 MPa^-0.5 (use zero to keep valve closed)
% akamodel.C2max = 50;			% valve constant for valve 2, cm^3 s^-1 MPa^-0.5
% akamodel.B1	= 0.01;			% valve drag, MPa s cm^-3
% akamodel.B2 	= 0.01;			% same for valve 2
% akamodel.L		= 0;			% actuator leakage parameter, cm^3 s^-1 (Nm)^-1
% akamodel.k		= 3.0;			% accumulator stiffness, MPa cm^-3
% 
% % solve the gait for the AKA model
% disp('Starting optimal control solution process for AKA model...');
% problem.model 			= akamodel;
% problem.initialguess 	= 'Winter/Winter_slow_result_bka_10000.mat';
% problem.resultfile 		= [movement '_result_aka_10000.mat'];
% optim(problem);
% % disp('Hit ENTER to continue with next optimization, or CTRL-C to quit');
% % pause;
% 
% % 
disp('All optimizations completed.')

