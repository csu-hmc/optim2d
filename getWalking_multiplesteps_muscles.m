% script for running the optim.m program
% clear all
rng('shuffle')

tic
% define which movement we are working on
movement = 'Winter/Winter_normal';

% general settings
problem.Solver = 'IPOPT';
problem.MaxIterations = 10000;
problem.ConstraintTol = .01;
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
ablemodel.Wtrack 	= 0;	% weight of tracking term in optimization objective
ablemodel.Weffort 	= 10;	% weight of muscle effort term in optimization objective
ablemodel.Wreg = 10;       % weight of the regularization term in optimization objective
ablemodel.effort.fatigue 	= 0;
ablemodel.effort.Fmaxweighted = 0;
ablemodel.effort.exponent	= 2;
ablemodel.Wvalve 	= 0.00;	% weight of valve operating cost in optimization objective
ablemodel.discretization = 'euler';
ablemodel.reducedW = 0;

% Starts with walking model that is found using
% getWalking_IPOPT_deterministic

Ncycles = 10;
stdev = [0.1 1 10 20];
% stdev = [0.1 0.5 1 5 10 20];% 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300];

iniguess = 'Winter/Winter_normal_result_able_notracking';%
%Multiple steps without tta
for k = 1:length(Ncycles)
    for j = 1:1
        problem.Ncycles = Ncycles(k);
        %start with 1 cycle solution
        problem.resultfile = iniguess;
        %new random data
        rndval = (rand(4,problem.N*problem.Ncycles)-0.5)*2; %randn(6,problem.N*problem.Ncycles);%  
        problem.warmstart = 0;
        for i = 1:length(stdev)
            ablemodel.rndval = stdev(i)*rndval;

            % solve the gait for the BKA model
            disp('Starting optimal control solution process for able model...');
            problem.model 			= ablemodel;
            problem.initialguess 	= problem.resultfile; % 'Winter/Winter_normal_result_torque_1cycle.mat'; %'Winter/Winter_normal_result_able.mat';%'mid';%
            problem.resultfile 		= [movement '_effort_notrack_able_fromoldsol_' num2str(problem.Ncycles) 'cycles_std' num2str(stdev(i)) '_time' num2str(j) '.mat'];
            if exist(problem.resultfile, 'file')
                disp('Problem already solved')
                load(problem.resultfile)
                rndval = result.model.rndval/stdev(i);
            else
                optim(problem);
%                 keyboard
            end
            problem.warmstart = 1;
        end
        
    end
end


% problem.N = 60;
% problem.symmetry = 0;
% problem.initialguess = 'Winter/Winter_normal_result_notrack_walk.mat';% 'Winter/Winter_normal_result_stand.mat';
% problem.resultfile = 'Winter/Winter_normal_result_notrack_walk_nosym.mat';
% problem.model = ablemodel;
% if exist(problem.resultfile, 'file')
%     disp('Problem already solved')
% else
%     optim(problem);
% end
% 
% %Get TTA model
% % define a below knee amputee model, based on the able-bodied model
% bkamodel		= ablemodel;
% bkamodel.type	= 'bka';
% % 450 Nm/rad foot stiffness in able bodied subjects (Hansen et al., J Biomech 37: 1467–1474, 2004)
% bkamodel.anklestiffness = 600;		% stiffness (Nm/rad) of prosthetic ankle
% bkamodel.ankledamping = 15;         % damping (Nm/rad/s) of prosthetic ankle
% 
% problem.initialguess = problem.resultfile;
% problem.resultfile = 'Winter/Winter_normal_result_notrack_tta.mat';
% problem.model = bkamodel;
% if exist(problem.resultfile, 'file')
%     disp('Problem already solved')
% else
%     optim(problem);
% end
% 
% Ncycles = 15;
% stdev = [0.1 1 10 20];
% % stdev = [0.1 0.5 1 5 10 20];% 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300];
% 
% iniguess = problem.resultfile;
% 
% % j = 1;
% for k = 1:length(Ncycles)
%     for j = 1:1
%         problem.Ncycles = Ncycles(k);
%         %start with 1 cycle solution
%         problem.resultfile = iniguess;
%         %new random data
%         rndval = (rand(4,problem.N*problem.Ncycles)-0.5)*2; %randn(6,problem.N*problem.Ncycles);%  
%         problem.warmstart = 0;
%         for i = 1:length(stdev)
%             bkamodel.rndval = stdev(i)*rndval;
% 
%             % solve the gait for the BKA model
%             disp('Starting optimal control solution process for able model...');
%             problem.model 			= bkamodel;
%             problem.initialguess 	= problem.resultfile; % 'Winter/Winter_normal_result_torque_1cycle.mat'; %'Winter/Winter_normal_result_able.mat';%'mid';%
%             problem.resultfile 		= [movement '_effort_notrack_TTA_fromstanding_' num2str(problem.Ncycles) 'cycles_std' num2str(stdev(i)) '_time' num2str(j) '.mat'];
%             if exist(problem.resultfile, 'file')
%                 disp('Problem already solved')
%                 load(problem.resultfile)
%                 rndval = result.model.rndval/stdev(i);
%             else
%                 optim(problem);
%                 keyboard
%             end
%             problem.warmstart = 1;
%         end
%         
%     end
% end
% % 
% % % 
% % disp('All optimizations completed.')
% % 
