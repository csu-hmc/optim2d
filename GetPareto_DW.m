% script for running the optim.m program
% clear all

% define which movement we are working on
movement = 'Winter/Winter_normal';
   
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
problem.Printinterval = 10;
        
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
ablemodel.kneeconstraint = 1;
ablemodel.hipconstraint = 1;

% define a below knee amputee model, based on the able-bodied model
bkamodel		= ablemodel;
bkamodel.type	= 'bka';
% 450 Nm/rad foot stiffness in able bodied subjects (Hansen et al., J Biomech 37: 1467–1474, 2004)
bkamodel.anklestiffness = 600;		% stiffness (Nm/rad) of prosthetic ankle
bkamodel.ankledamping = 15;         % damping (Nm/rad/s) of prosthetic ankle

bkamodel.Wmom = 5;
bkamodel.Weffort = 2;
bkamodel.Wtrack = 0.1;

problem.model 			= bkamodel;
problem.initialguess 	= 'Winter/Winter_normal_result_able';%'Winter/Pareto_DW/Result_3_mom_1_effort_12_track_0.6.mat';
problem.resultfile 		= ['Winter/Pareto_DW/Result2_3_mom_', num2str(bkamodel.Wmom), '_effort_', num2str(bkamodel.Weffort), '_track_' num2str(bkamodel.Wtrack), '.mat'];
optim_threeobj(problem)



% bkamodel.Wtrack = 1;
% bkamodel.Weffort = 20*bkamodel.Wtrack;
% for j = 1:10
%     bkamodel.Wmom = j/10;
%     problem.model 			= bkamodel;
%     problem.initialguess 	= 'Winter\Winter_normal_result_able.mat';
%     problem.resultfile 		= ['Winter/Pareto_DW/Result_3_mom_', num2str(bkamodel.Wmom), '_effort_', num2str(bkamodel.Weffort), '_track_' num2str(bkamodel.Wtrack), '.mat'];
%     optim_threeobj(problem);
% end
% 
% bkamodel.Wmom = 1;
% for j = 1:9
%     bkamodel.Wtrack = j/10;
%     bkamodel.Weffort = 20*bkamodel.Wtrack;
%     problem.model 			= bkamodel;
%     problem.initialguess 	= 'Winter\Winter_normal_result_able.mat';
%     problem.resultfile 		= ['Winter/Pareto_DW/Result_3_mom_', num2str(bkamodel.Wmom), '_effort_', num2str(bkamodel.Weffort), '_track_' num2str(bkamodel.Wtrack), '.mat'];
%     optim_threeobj(problem);
% end


