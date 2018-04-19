% script for running the optim.m program
% clear all
rng('shuffle')

tic
% define which movement we are working on
movement = 'Winter/Winter_normal';

% general settings
problem.Solver = 'IPOPT';
problem.MaxIterations = 10000;
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
problem.Printinterval = 10;
problem.proscontrol = 0;
   
% define an able-bodied model and the target gait data for the simulation
ablemodel.parameterfile 	= 'gait2d_par.xls';
ablemodel.type 				= 'able';
ablemodel.datafile			= movement;
ablemodel.Wtrack 	= 10;	% weight of tracking term in optimization objective
ablemodel.Weffort 	= 0.01;	% weight of muscle effort term in optimization objective
ablemodel.Wreg = 1;       % weight of the regularization term in optimization objective
ablemodel.effort.fatigue 	= 0;
ablemodel.effort.Fmaxweighted = 0;
ablemodel.effort.exponent	= 2;
ablemodel.Wvalve 	= 0.00;	% weight of valve operating cost in optimization objective
ablemodel.discretization = 'euler';
ablemodel.reducedW = 0;

% Starts with walking model that is found using
% getWalking_IPOPT_deterministic
% 
% define a below knee amputee model, based on the able-bodied model
torquemodel		= ablemodel;
% torquemodel.dur = 1.1;
% torquemodel.speed = torquemodel.dur*1.1627;
torquemodel.type	= 'torque';
% 
%solve the gait for the BKA model
disp('Starting optimal control solution process for torque model...');
problem.model 			= torquemodel;
problem.initialguess 	= [movement '_result_able.mat'];%'mid';%
problem.resultfile 		= [movement '_result_torque_1cycle_effort_test.mat'];
if exist(problem.resultfile, 'file')
    disp('Problem already solved')
else
    optim(problem);
end
% keyboard
% problem.Ncycles = 10;
% torquemodel.rndval = 10*rand(6,problem.N,problem.Ncycles);
% % solve the gait for the BKA model
% disp('Starting optimal control solution process for torque model...');
% problem.model 			= torquemodel;
% problem.initialguess 	= problem.resultfile; % 'Winter/Winter_normal_result_torque_1cycle.mat'; %'Winter/Winter_normal_result_able.mat';%'mid';%
% problem.resultfile 		= [movement '_result_torque_' num2str(problem.Ncycles) 'cycles.mat'];
% optim(problem);


% Ncycles = [2 5 10 20 50];
Ncycles = [6 7 9];%[2 5];
stdev = [0.1 0.5 1 5 10];% 12 15 20 30 40 50];% 60 65 70 75 80 90 100];% 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300];

iniguess = problem.resultfile;
% j = 1;
for k = 1:length(Ncycles)
    for j = 1:5
        problem.Ncycles = Ncycles(k);
        %start with 1 cycle solution
        problem.resultfile = iniguess;
        %new random data
        rndval = (rand(6,problem.N*problem.Ncycles)-0.5)*2; %randn(6,problem.N*problem.Ncycles);% 
        for ii = 1:size(rndval,1)
            for jj = 1:size(rndval,2)
                if abs(rndval(ii,jj)) > 3
                    rndval(ii,jj) = sign(rndval(ii,jj))*3;
                end
            end
        end  
        problem.warmstart = 0;
        for i = 1:length(stdev)
            torquemodel.rndval = stdev(i)*rndval;

            % solve the gait for the BKA model
            disp('Starting optimal control solution process for torque model...');
            problem.model 			= torquemodel;
            problem.initialguess 	= problem.resultfile; % 'Winter/Winter_normal_result_torque_1cycle.mat'; %'Winter/Winter_normal_result_able.mat';%'mid';%
            problem.resultfile 		= [movement '_torque_rand_' num2str(problem.Ncycles) 'cycles_std' num2str(stdev(i)) '_time' num2str(j) '.mat'];
            if exist(problem.resultfile, 'file')
                disp('Problem already solved')
                load(problem.resultfile)
                rndval = result.model.rndval/stdev(i);
            else
                optim(problem);
            end
            problem.warmstart = 1;
        end
        
    end
end
% 
% % 
% disp('All optimizations completed.')
% 
