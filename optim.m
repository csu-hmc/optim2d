function [result] = optim(prob);

% predictive optimization

	tic;
	clear global
	clear mex
	global problem
	
	problem = prob;
	problem.print = 1;						% 0: never, 1: at Printinterval, 2: always
		
	%----------------------------------------------------------------------------------
	% Number of nodes
	N			= problem.N;				% number of collocation nodes
    if ~isfield(problem, 'Ncycles')
        problem.Ncycles = 1;
    end
    Ncycles     = problem.Ncycles;          % number of gait cycles
	if mod(N,2)~=0
		error('Number of nodes must be even.');
	end
	
	%----------------------------------------------------------------------------------
	% Some model related constants
	ndof = 9;
    nprosstates = 0;
    if strcmp(problem.model.type, 'torque')
        nmus = 0;
        ncontrols = 6;
        
    else
        nmus = 16;
        ncontrols = nmus;
        if problem.proscontrol
            nprosstates = 4;
            ncontrols = ncontrols+2;
        end
    end
    nstates = 2*ndof + 2*nmus + nprosstates;
	nvarpernode1 = nstates + ncontrols;                         % number of unknowns per time node, first gait cycle
    nvarpernode = nstates;
	nvar = nvarpernode1*N + nvarpernode*N*(Ncycles-1) + Ncycles + 2;	% total number of unknowns (states, controls, duration, feedback gain)
	ncon = nstates * N * Ncycles;       						% number of constraints due to discretized dynamics and task
    if ~isfield(problem.model,'kneeconstraint')
        problem.model.kneeconstraint = 0;
    end
    if ~isfield(problem.model,'hipconstraint')
        problem.model.hipconstraint = 0;
    end
    if ~isfield(problem.model,'Wreg')
        problem.model.Wreg = 0;
    end
    if ~isfield(problem.model, 'rndval')
        problem.model.rndval = zeros(6,problem.N*problem.Ncycles);
    end
    if problem.model.kneeconstraint == 1
        ncon = ncon + N;
    end
    if problem.model.hipconstraint == 1
        ncon = ncon+N;
    end
    if ~isfield(problem, 'warmstart')
        problem.warmstart = 0;
    end
    %load warm start settings
    if and(strcmp(problem.Solver,'IPOPT'),problem.warmstart);
        load(problem.initialguess);
        if result.problem.Ncycles == 1
            zu = reshape(result.zu(1:end-3),result.problem.nvarpernode1,result.problem.N);
            zl = reshape(result.zl(1:end-3),result.problem.nvarpernode1,result.problem.N);
            zux = zu(1:nstates,:);
            zlx = zl(1:nstates,:);
            problem.zu = [result.zu(1:end-3);repmat(zux(:),Ncycles-1,1);repmat(result.zu(end-2),Ncycles,1);result.zu(end-1:end)];
            problem.zl = [result.zl(1:end-3);repmat(zlx(:),Ncycles-1,1);repmat(result.zl(end-2),Ncycles,1);result.zl(end-1:end)];
            problem.llambda = repmat(result.lambda,problem.Ncycles,1);
        else
            problem.zu = result.zu;
            problem.zl = result.zl;
            problem.llambda = result.lambda;
        end
        clear result
    end

	problem.ndof = ndof;
	problem.nmus = nmus;
	problem.nstates = nstates;
	problem.ncontrols = ncontrols;
	problem.nvarpernode = nvarpernode;
    problem.nvarpernode1 = nvarpernode1;
	problem.nvar = nvar;
	problem.ncon = ncon;
	
	% precompute the indices for controls, states, first derivatives
	iu = [];
    iq = [];
    iqd= [];
    if strcmp(problem.model.type,'torque')
        % First gait cycle
        for i=0:N-1
            iu = [iu  nvarpernode1*i+nstates+(1:ncontrols)];
            iq = [iq  nvarpernode1*i+(4:ndof)];
            iqd= [iqd nvarpernode1*i+ndof+(4:ndof)];
        end
        % other gait cycles
        for j = 0:Ncycles-2
            for i=0:N-1
                iq = [iq  nvarpernode1*N+nvarpernode*N*j+nvarpernode*i+(4:ndof)];
                iqd= [iqd nvarpernode1*N+nvarpernode*N*j+nvarpernode*i+ndof+(4:ndof)];
            end
        end
    else
        for i=0:N-1
            iu = [iu nvarpernode*i+nstates+(1:nmus)];
        end
    end
	problem.iumus = iu;
    problem.iq = iq;
    problem.iqd= iqd;
    
    %make mirror for symmetry
    problem = makemirror(problem);
    
	%----------------------------------------------------------------------------------
	% Initialize one or two models and data
	problem.model.gait2d = @gait2d;				% create function handle for the musculoskeletal model
	problem.model = initialize(problem.model);
	if isfield(problem, 'model1')
		original = which('gait2d.mexw32');				% the original gait2d.mexw32 should be in search path
		system(['copy ' original ' gait2d_1.mexw32']); 	% make a copy in current directory, called gait2d_1.mexw32
		problem.model1.gait2d = @gait2d_1;				% create function handle for the musculoskeletal model
		problem.model1 = initialize(problem.model1);
	else
		problem.lambda = 1;			% to ensure that model1 is not used
	end

	%----------------------------------------------------------------------------------
	% Initialize logging mechanism
	problem.log = [];
	
	%----------------------------------------------------------------------------------
	% Set lower and upper bounds for the optimization variables X
	L = zeros(nvar,1);
	U = zeros(nvar,1);
	% bounds for musculoskeletal controls
	
    Lqqdot = [-1; 0.5; -pi/4; -pi; -pi; -pi; -pi; -pi; -pi; 		... % q
            -200; -200;  -500; -500; -500; -500; -500; -500; -500];  	% qdot
    Uqqdot = [4*Ncycles; 1.5; pi/4; pi; pi; pi; pi; pi; pi; 				... % q
            500; 500; 500; 500; 500; 500; 500; 500; 500];
    if strcmp(problem.model.type, 'torque')
        Lu = -1000+zeros(ncontrols,1);
        Uu = 1000+zeros(ncontrols,1);
%         Lu([3 6]) = -100;
%         Uu([3 6]) = 100;
        Lx = Lqqdot;
        Ux = Uqqdot;
    else
        Lu = zeros(nmus,1);
        Uu = 2*ones(nmus,1);
    
        % bounds for musculoskeletal states
        Lx = [Lqqdot; zeros(nmus,1)-1;									...	% Lce
            Lu]; 											...	% active states

        Ux = [Uqqdot; zeros(nmus,1)+5; 									...	% Lce
            Uu];											... % active states
    end
   
	% bounds for prosthesis states
    if problem.proscontrol
        Lxp = [-100 ; -100; -100; -500];			% s, v1, v2, M
        Uxp = [100 ; 100; 100; 500];
        % bounds for prosthesis controls
        Lup = zeros(2,1)+0.001;
        Uup = ones(2,1);    
    else
        Lxp = []; Uxp = []; Lup = []; Uup = [];
    end
	
    % First gait cycle
    for k = 0:N-1
		L(k*nvarpernode1 + (1:nvarpernode1) ) = [Lx; Lxp; Lu; Lup];
		U(k*nvarpernode1 + (1:nvarpernode1) ) = [Ux; Uxp; Uu; Uup];
    end
    for j = 0:Ncycles-2
        for k = 0:N-1
            L(nvarpernode1*N+nvarpernode*N*j+nvarpernode*k+(1:nvarpernode)) = [Lx; Lxp];
            U(nvarpernode1*N+nvarpernode*N*j+nvarpernode*k+(1:nvarpernode)) = [Ux; Uxp];
        end
    end
    L(end-1-Ncycles:end-2) = 0.2;		% minimum duration of movement cycle, for each gait cycle
    U(end-1-Ncycles:end-2) = 2.0;		% maximum duration of movement cycle    
    if problem.model.rndval ~= 0
        L(end-1:end) = -10000;               % Feedback gain
        U(end-1:end) = 10000;                % Derivative Feedback gain
    else
        % No feedback
        L(end-1:end) = 0;
        U(end-1:end) = 0;
    end
	% constrain X of trunk at first node to be zero (since it does not matter, and this helps convergence)
	L(1) = 0;
	U(1) = 0;
	
	%----------------------------------------------------------------------------------
	% generate constraint scaling factors
	% we suggest divide by 500 for the implicit equations of motion (these are moment imbalances)
	Wc = ones(ncon,1);
	ieom = ndof+(1:ndof);		% index of equations of motion in first node
	for i=1:N
		Wc(ieom) = 1/500;
		ieom = ieom + nstates;
	end
	problem.Wc = spdiags(Wc,0,ncon,ncon);	

	%----------------------------------------------------------------------------------
	% make an initial guess
	if strcmp(problem.initialguess, 'mid')
		X0 = (L + U)/2; %	zeros(size(L))+1e-6;%							% halfway between upper and lower bound 
%         if ~strcmp(problem.model.type,'torque')
%             for k = 1:N
%                 X0((k-1)*72+1) = X0(end);
%             end
%         end
		% X0 = X0 + 0.01*(U-L).*randn(size(X0));	% to add some randomness
% 		X0 = L + (U-L).*rand(size(L));			% to make the initial guess completely random
	else
		% load a previous solution
		load(problem.initialguess);
		Nresult = size(result.u,2);
		t0 = (0:(Nresult-1))'/Nresult;
		x0 = result.x';
		u0 = result.u';
        durchange = 0; %boolean required for symmetry
        if (size(x0,1) ~= Nresult) %A solution with multiple gait cycles is used
            %assume that the initial guess has same number of cycles and
            %time nodes etc as current problem
            X0 = reshape([x0(1:Nresult,:) u0]',nvarpernode1*N,1);
            X1 = reshape(x0(Nresult+1:end,:)',nvarpernode*N*(Ncycles-1),1); %just states for the other gait cycles
            X0 = [X0;X1];
        else
            if and(strcmp(result.model.type,'able'),strcmp(problem.model.type, 'torque'))
                mom = zeros(Nresult,ncontrols);
                result.model = initmodel(result.model);
                for i = 1:Nresult
                    mom(i,:) = result.model.gait2d('Jointmoments',x0(i,1:50)');
                end
                problem.model = initmodel(problem.model);
                u0 = mom;
            end
            if ~isfield(result, 'problem')
                %This is from before symmetry was implemented, so always
                %going to be a full gait cycle
                if problem.symmetry
                    t0 = (0:2:(Nresult-1))'/Nresult;
                    x0 = x0(1:Nresult/2,1:nstates);
                    u0 = u0(1:Nresult/2,1:ncontrols);
                    result.dur = result.dur/2;
                    durchange = 1;
                end
            elseif and(~result.problem.symmetry,problem.symmetry)
                t0 = (0:2:(Nresult-1))'/Nresult;
                x0 = x0(1:Nresult/2,1:nstates);
                u0 = u0(1:Nresult/2,1:ncontrols);
                result.dur = result.dur/2;
                durchange = 1;
                %todo: implement symmetry -> no symmetry
            end
            
            % Do not use muscle states and less control states using torque
            % model
%             if and(strcmp(problem.model.type, 'torque'),~strcmp(result.model.type, 'torque'))
%                 x0 = x0(:,1:ndof*2);
%                 u0 = u0(:,1:ncontrols);
%             end
            
            % duplicate first node node at the end so we can interpolate with periodicity
            t0 = [t0 ; 1];
            if problem.symmetry
                u0 = [u0; u0(1,problem.vmu)];
                x0 = [x0; x0(1,problem.vmx)];
            else
                u0 = [u0 ; u0(1,:)];
                x0 = [x0 ; x0(1,:)];
            end
%             x0(end,1) = x0(end,1) + result.speed * result.dur(1);

            % interpolate states and controls from initial guess to current node times
            times = (0:(N-1))/N;
            x0 = interp1(t0,x0,times,'linear','extrap');
            u0 = interp1(t0,u0,times,'linear','extrap');
            X0 = reshape([x0(:,1:nstates) u0(:,1:ncontrols)]',nvarpernode1*N,1);
            if Ncycles > 1
                for j = 1:Ncycles-1
                    if length(result.dur) > 1
                        dur = result.dur(j);
                    else
                        dur = result.dur;
                    end
                    x0(:,1) = x0(:,1)+result.speed*dur;
                    X1 = reshape(x0(:,1:nstates)',nstates*N,1); %just states for the other gait cycles
                    X0 = [X0;X1];
                end
            end
        end
        if durchange
            result.dur = result.dur*2;
        end
        if length(result.dur) > 1                
            X0 = [X0;result.dur];
            X0 = [X0;result.K]; %feedback gains
        else
            X0 = [X0;repmat(result.dur,Ncycles,1)]; %durations
            X0 = [X0;0;0]; %feedback gains
        end
        
	end
	problem.X0 = X0;			% store initial guess in case we need it later
	
	% determine sparsity structure of Jacobian
	% we have verified that we always get same structure by testing with random X, so do this only once
	for i=1:1
		problem.Jnnz = 1;
		X = L + (U-L).*rand(size(L));		% a random vector of unknowns
		J = conjac_2d(X);
		problem.Jnnz = nnz(J);
		fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n',problem.Jnnz, ncon*nvar, 100*problem.Jnnz/(ncon*nvar));
		problem.Jpattern = double(J~=0);
	end
	
	%----------------------------------------------------------------------------------
	% check the derivatives
	if problem.checkderivatives
		% check the model derivatives
		hh = 1e-7;
		x = randn(nstates,1);
		xdot = randn(nstates,1);
		u = randn(ncontrols,1);
		[f, dfdx, dfdxdot, dfdu] = dyn(problem.model, x, xdot, u);
		dfdx_num = zeros(size(dfdx));
		dfdxdot_num = zeros(size(dfdxdot));
		dfdu_num = zeros(size(dfdu));
		for i=1:nstates
			tmp = x(i);
			x(i) = x(i) + hh;
			fhh = dyn(problem.model, x, xdot, u);
			dfdx_num(:,i) = (fhh-f)/hh;
			x(i) = tmp;
			
			tmp = xdot(i);
			xdot(i) = xdot(i) + hh;
			fhh = dyn(problem.model, x, xdot, u);
			dfdxdot_num(:,i) = (fhh-f)/hh;
			xdot(i) = tmp;
		end
		for i=1:ncontrols
			tmp = u(i);
			u(i) = u(i) + hh;
			fhh = dyn(problem.model, x, xdot, u);
			dfdu_num(:,i) = (fhh-f)/hh;
			u(i) = tmp;
		end
		% report maximal differences between analytical derivatives and numerical results
		fprintf('Max. error in dfdx: ');
		matcompare(dfdx, dfdx_num);
		fprintf('Max. error in dfdxdot: ');
		matcompare(dfdxdot, dfdxdot_num);
		fprintf('Max. error in dfdu: ');
		matcompare(dfdu, dfdu_num);
% 		keyboard
	
		% check the NLP derivatives
		hh = 1e-7;
		X = L + (U-L).*rand(size(L));		% a random vector of unknowns
		f = objfun_2d(X);
		grad = objgrad_2d(X);
		c = confun_2d(X);
		cjac = conjac_2d(X);
		cjac_num = zeros(ncon,nvar);
		grad_num = zeros(nvar,1);
		for i=1:nvar
			fprintf('checking objgrad and conjac for unknown %4d of %4d\n',i,nvar);
			Xisave = X(i);
			X(i) = X(i) + hh;
			cjac_num(:,i) = (confun_2d(X) - c)/hh;
			grad_num(i) =   (objfun_2d(X) - f)/hh;
			X(i) = Xisave;
		end
		
		% report maximal differences between analytical derivatives and numerical results
		fprintf('Max. error in constraint jacobian: ');
		matcompare(cjac, cjac_num);
		fprintf('Max. error in objective gradient: ');
		matcompare(grad, grad_num);
		keyboard
	end
	
	%------------------------------------------------------------------------------------
	screen = get(0,'ScreenSize');
	close all;
	figure(1);		% all plotting during optimization is done in this figure
	clf;
	set(gcf,'OuterPosition',[1 screen(4)-300 screen(3)-200 300]);
	problem.log = [];							% clear the log

	%------------------------------------------------------------------------------------
	% evaluate initial guess
	fprintf('Initial guess evaluation:\n');
	tempprint = problem.print;
    problem.print = 2;
	evaluate_if_needed(X0);
    problem.print = tempprint;
	if (problem.debug)
		disp('Hit ENTER to start optimization...');pause
	end

	%----------------------------------------------------------------------------------------
	% run optimization
	result.solver = problem.Solver;
	if (problem.MaxIterations > 0)	
		fprintf('Starting optimization...\n');
		if strcmp(problem.Solver,'IPOPT')
			funcs.objective = @objfun_2d;
			funcs.gradient  = @objgrad_2d;
			funcs.constraints = @confun_2d;
			funcs.jacobian    = @conjac_2d;
			funcs.jacobianstructure = @conjacstructure_2d;
			options.lb = L;
			options.ub = U;
			options.cl = [zeros(ncon,1)];
			options.cu = [zeros(ncon,1)];
			options.ipopt.max_iter = problem.MaxIterations;
			options.ipopt.hessian_approximation = 'limited-memory';
			options.ipopt.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
			options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
			options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
			options.ipopt.bound_push = options.ipopt.bound_frac;
            options.ipopt.tol = problem.Tol;
            options.ipopt.linear_solver = 'mumps';
            options.ipopt.dual_inf_tol = problem.Tol;
			options.ipopt.constr_viol_tol = problem.ConstraintTol;
			options.ipopt.compl_inf_tol = problem.ConstraintTol;
            options.ipopt.print_level = 0;
            if problem.warmstart
                options.ipopt.warm_start_init_point = 'yes';
%                 options.ipopt.warm_start_bound_frac = 1e-16;
%                 options.ipopt.warm_start_bound_push = 1e-16;
%                 options.ipopt.warm_start_mult_bound_push = 1e-16;
%                 options.ipopt.warm_start_slack_bound_frac = 1e-16;
%                 options.ipopt.warm_start_slack_bound_push = 1e-16;
                options.zl = problem.zl;
                options.zu = problem.zu;
                options.lambda = problem.llambda;
            end 
			[X, info] = ipopt(X0,funcs,options);
            
			result.info = info.status;
			result.message = ipoptmessage(info.status);
            result.obj = objfun_2d(X);
            result.zl = info.zl;
            result.zu = info.zu;
            result.lambda = info.lambda;
            result.K = X(end-1:end);
		elseif strcmp(problem.Solver,'SNOPT')
            problem.objconfun = @objconfun;
            % Change options
            testspec.spc = which('testspec.spc');
            snspec ( testspec.spc );
            % Output informative files
            snprint   ('probName.out');
            snsummary ('prName.sum');
            if problem.warmstart == 1;
                snset ('Warm start')
            end
            FL = [-inf;zeros(problem.ncon,1)];
            FU = [inf;zeros(problem.ncon,1)];
            xmul = zeros(size(L));
            Fmul = zeros(size(FL));
            xstate = 2*ones(size(X0));
            Fstate = 2*ones(size(FL));
            [X,F,INFO] = snopt(X0,L,U,xmul, xstate,FL,FU,Fmul, Fstate, 'objconfun');
%             [X,F,INFO] = snopt(X0,L,U,FL,FU, 'objconfun');
            snprint    off;
            snsummary off;
            result.info = INFO;
            result.obj = F(1);
            result.X = X;
            [result.f, result.g, result.c, result.J] =evaluate(problem.model,X);
            if result.info == 1
                result.message = 'Optimization Solved';
            else
                result.message = 'Check INFO';
            end
        elseif strcmp(problem.Solver,'TOMLAB')
            Prob = conAssign(@objfun_2d, @objgrad_2d, [], [], L, U, 'optim2d', X0, ...
                            [], 0, ...
                            [], [], [], @confun_2d, @conjac_2d, [], problem.Jpattern, ...
							zeros(ncon,1), zeros(ncon,1), ...
                            [], [], [],[]);
			% Prob.SOL.optPar(1)= 1;			% uncomment this to get snoptsum.txt and snoptpri.txt
			Prob.SOL.optPar(9) = problem.ConstraintTol;		% feasibility tolerance
			Prob.SOL.optPar(10) = problem.Tol;				% optimality tolerance
			Prob.SOL.optPar(11) = 1e-6; % Minor feasibility tolerance (1e-6)
			Prob.SOL.optPar(30) = 1000000; % maximal sum of minor iterations (max(10000,20*m))
			Prob.SOL.optPar(35) = problem.MaxIterations;
			Prob.SOL.optPar(36) = 40000; % maximal number of minor iterations in the solution of the QP problem (500)
			Prob.SOL.moremem = 10000000; % increase internal memory
%             Prob.PriLevOpt = 2; % Print every 10 major iterations
			Result = tomRun('snopt',Prob);
			X = Result.x_k;
			result.message = Result.ExitText;
			result.info = Result.ExitFlag;    
		else
			error('Solver name not recognized: %s', problem.Solver);
		end
	else		% skip optimization
		X = X0;
		result.info = 0;		% indicate success even if we did not optimize
    end
%     keyboard
    result.problem = problem;
	disp(['Solver completion status: ' result.message]);
    disp(problem.model.type)
	fprintf('Duration of movement: %8.4f s\n',sum(X(end-1-Ncycles:end-2)));
	
	% save optimization result on file
	savefile(X, problem.resultfile, result);
	disp(['Result was saved in the file ' problem.resultfile]);
	
	% show stick figure of result and optimization log
	disp('Result of optimization:');
	problem.print = 2;
	evaluate_if_needed(X);
	
end						% end of function "optim"
%===========================================================================================
function model = initialize(model);
	global problem

	% Initialize the musculoskeletal model
	model = initmodel(model);
	
	%----------------------------------------------------------------------------------
	% Determine number of nonzeros in model Jacobians
	% we have verified that we always get same structure by testing with random inputs, so do this only once
	model.nnz_dfdx = 1;
	model.nnz_dfdxdot = 1;
	model.nnz_dfdu = 1;
	x 		= randn(problem.nstates,1);
	xdot 	= randn(problem.nstates,1);
	u 		= randn(problem.ncontrols,1);
	[f,dfdx,dfdxdot,dfdu] = dyn(model,x,xdot,u);
	model.nnz_dfdx 		= nnz(dfdx);
	model.nnz_dfdxdot 	= nnz(dfdxdot);
	model.nnz_dfdu 		= nnz(dfdu);
	fprintf('Number of nonzeros in model Jacobians: %d %d %d\n', model.nnz_dfdx, model.nnz_dfdxdot, model.nnz_dfdu);
	
	% Initialize the data for this model
	% read gait data from file created by preproc.m
	fprintf('Reading subject movement data from file %s...\n', model.datafile);
	load(model.datafile);
    
    if isfield(model,'speed');
        gait.speed = model.speed;
    end

    if isfield(model,'dur')
        gait.dur(1) = model.dur;
    end
    
	% copy prescribed speed and target gait cycle duration into objdata struct
	fprintf('Gait speed:    %9.4f m/s\n',gait.speed(1));  
	fprintf('Gait period:   %9.4f  +- %9.4f s\n',gait.dur);
	
	% make gait data periodic and apply appropriate smoothing
	nterms = 6;										% number of fourier terms
	% I have turned this off because it distorts the GRF, try something else later.
	% perdata = makeperiodic(gait.data, nterms);	
	perdata = gait.data;
	
	% interpolate gait data at each collocation node
    N = problem.N;
    if problem.symmetry
        N = 2*N;
    end
	Ncycle = size(gait.data,1)-1;		% number of data samples per gait cycle (data file includes the 100% gait cycle point!)
	tcycle = (0:Ncycle)/Ncycle; 		% time vector for gait data (in fraction of gait cycle)
	tresamp = (0:(N-1))/N;				% we resample the gait cycle into N samples (because N is gait cycle)
	av = interp1(tcycle, perdata, tresamp);
	sd = interp1(tcycle, gait.sd, tresamp);
	
	% Convert kinematic data to our gait2d model generalized coordinates and to radians
	iang = 1:3;							% columns containing the hip, knee, and ankle angles
	av(:,iang) = av(:,iang)*pi/180;			
	sd(:,iang) = sd(:,iang)*pi/180;
	av(:,2) = -av(:,2);					% knee angle must be inverted
	
	% HACK: increase the SD for ankle angle, to relax the tracking there
	% sd(:,3) = 10*sd(:,3);

	% Convert GRF to our coordinate system and to units of body weight
	BW = 2*mean(av(:,5));			% body weight from vertical GRF
	igrf = 4:5;
	av(:,igrf) = av(:,igrf)/BW;
	sd(:,igrf) = sd(:,igrf)/BW;
	
	% prevent very low SD values
	sd(:,iang) = max(sd(:,iang), 2*pi/180);		% 2 deg minimum SD for angles
	sd(:,igrf) = max(sd(:,igrf), 0.02);			% 2% BW minimum SD for forces
	if (gait.dur(2) < 1e-4)
		disp('Warning: gait cycle duration had small SD.  Using SD = 0.5 s instead.');
		gait.dur(2) = 0.5;
	end

	% store average and SD in N x 10 matrix
    if problem.symmetry
        av = [av(1:N/2,:) av(N/2+1:N,:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
        sd = [sd(1:N/2,:) sd(N/2+1:N,:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
    else
        av = [av(1:N,:) av([N/2+1:N 1:N/2],:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
    	sd = [sd(1:N,:) sd([N/2+1:N 1:N/2],:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
    end
	
	% store everything in the data struct within model
	model.data.av = av;
	model.data.sd = sd;
	model.data.dur = gait.dur(1);
	model.data.dursd = gait.dur(2);
	model.data.speed = gait.speed(1);
	
	% the 'zero' model can't move, so we set prescribed speed to zero
	if strcmp(model.type, 'zero')
		model.data.speed = 0.0;
    end

end
%===========================================================================================
function [f, g, c, J] = evaluate(model, X)
	% evaluates the objective, constraints, and their derivatives
	global problem
	
	av = model.data.av;					% target angles and GRF
	sd = model.data.sd;					% and their SD
	dur = model.data.dur;				% target duration
	dursd = model.data.dursd;			% and its SD
	speed = model.data.speed;			% prescribed speed
	
	N = problem.N;
    Ncycles = problem.Ncycles;
	nstates = problem.nstates;
	ncontrols = problem.ncontrols;
	nvarpernode = problem.nvarpernode1;
    ndof = problem.ndof;
	nmus = problem.nmus;
	ncon = problem.ncon;
	nvar = problem.nvar;
	nxg = 2*problem.ndof + 2*nmus;		% number of musculoskeletal states
	iumus = problem.iumus;
    iq = problem.iq;
    iqd= problem.iqd;
	
	BW = 9.81*model.mass;
	gait2d = model.gait2d;
	
	h = X(end-1-Ncycles:end-2)/N;		% time step size, for each gait cycle
    if problem.symmetry
        h = h/2;
    end
	if strcmp(model.discretization, 'midpoint')
		discr = 1;
	elseif strcmp(model.discretization, 'euler')
		discr = 2;
	else
		error('optim:evaluate error: unknown discretization option.');
	end

	%--------------------------------------------------------------
	% Objective function f and its gradient g
	g = zeros(size(X));

	% tracking term for kinematics and GRF
	ixg = 1:nxg;				% index to musculoskeletal states in node 1
	f1 = 0;    
    notrack = size(av,2)*Ncycles+Ncycles; %Ncycles;%size(av,2)*Ncycles;%
    for j = 1:Ncycles
        for i=1:N
            if strcmp(model.type, 'torque')
                x(1:ndof*2,:) = X(ixg);
                x(ndof*2+1:50) = zeros(32,1);
            else
                x = X(ixg);
            end
            [GRF, dGRFdx] = gait2d('GRF',x);
            GRF = GRF/BW;							% convert to units of body weight
            dGRFdx = dGRFdx/BW;						% here also
            if strcmp(model.type, 'torque')
                dGRFdx(:,ndof*2+1:end)= [];
            end
            simdata = [x(4:6); GRF(1:2); x(7:9); GRF(3:4)]';	% the ten tracked variables: right angles, right GRF, left angles, left GRF
            res = (simdata-av(i,:))./sd(i,:);		% the ten residuals, normalized to SD
            f1 = f1 + sum(res.^2)/(notrack*N);
            res = res./sd(i,:);						% divide again by SD, needed for gradient
            g(ixg(4:6)) = g(ixg(4:6)) + 2*res(1:3)'/(notrack*N);		% gradient of objective with respect to right side angles
            g(ixg(7:9)) = g(ixg(7:9)) + 2*res(6:8)'/(notrack*N);		% gradient of objective with respect to left side angles
            g(ixg) = g(ixg) + 2*dGRFdx'*res([4 5 9 10])'/(notrack*N);	% gradient of GRF terms in objective, with respect to all state variables
            
            ixg = ixg + nvarpernode;								% move pointers to next node
            if i == N
                nvarpernode = problem.nvarpernode;
            end
        end
        
        % Always track duration with torque model
%         if strcmp(problem.model.type,'torque')
%             Wdur = max(1,model.Wtrack);
%             if model.Wtrack < 1e-4
%                 notrack = Ncycles;
%             end        
            Wdur = 1; %1 BECAUSE WTRACK COMES AFTER
%         else
%             Wdur = model.Wtrack;
%         end
        	
    end
    
    Wdur = 1;
    %Track complete duration instead of individual
    f1 = f1 + Wdur*((sum(X(end-2-Ncycles+(1:Ncycles))) - dur*Ncycles)/dursd)^2/notrack;				
    g(end-2-Ncycles+(1:Ncycles)) = g(end-2-Ncycles+(1:Ncycles)) + Wdur*2*(sum(X(end-2-Ncycles+(1:Ncycles))) - dur*Ncycles)/dursd^2/notrack;	
    
    % apply weighting to the tracking term
    f1 = f1*model.Wtrack;
    g = g*model.Wtrack;


	% effort term for muscles
	f2 = 0;
	expon = abs(model.effort.exponent);
    nvarpernode = problem.nvarpernode1;
    if strcmp(model.type, 'torque')
        W = model.Weffort;%/100;
        ix = 1:nstates;
        for j = 1:Ncycles
            iu = nstates+(1:ncontrols);
            for i = 1:N %inputs only first gait cycle
                u1 = X(iu);
                x1 = X(ix);
                u1 = u1+x1(4:ndof)*X(end-1)*10+x1(ndof+(4:ndof))*X(end)*10;     % feedback isnot on first three states, or last two controls
                f2 = f2 + 1/(ncontrols*N*Ncycles)*sum(u1.^expon);
                
                dfdu = expon/(ncontrols*N*Ncycles)*W*u1.^(expon-1);

                g(iu) = g(iu) + dfdu;
                g(ix) = g(ix) + [zeros(3,1);dfdu*X(end-1)*10;zeros(3,1);dfdu*X(end)*10];
                g(end-1)= g(end-1)+ 10*sum(dfdu.*x1(4:ndof));	
                g(end)= g(end)+ 10*sum(dfdu.*x1(ndof+(4:ndof)));	
                iu = iu + problem.nvarpernode1;			% move pointers to next node
                ix = ix + nvarpernode;

                if i == N
                    nvarpernode = problem.nvarpernode;
                end
            end
        end
        f2 = f2*W;
    else
        % More than one gait cycle is not implemented for model with muscles
        if Ncycles == 1
            for i=1:nmus
                if model.effort.Fmaxweighted
                    W = model.Weffort * model.Fmax(i) / nmus / sum(model.Fmax);
                else
                    W = model.Weffort / nmus;
                end
                if model.reducedW && (i>5) && (i<9)
                    W = W*0.01;
                end 
                iu = nstates + i + nvarpernode*(0:N-1);
                if model.effort.fatigue
                    % fatigue-like effort calculation
                    meanact = mean(X(iu));
                    f2 = f2 + W * meanact^expon;
                    g(iu) = g(iu) + expon * W * meanact^(expon-1) / N;
                else
                    % simple mean of activation^expon over all muscles and nodes
                    f2 = f2 + W * mean(X(iu).^expon);
                    g(iu) = g(iu) + expon * W * X(iu).^(expon-1) / N;
                end
            end
        else
            W = model.Weffort/nmus;
            ix1 = 1:nstates;
            for iCyc = 1:Ncycles
                iu1 = nstates+(1:ncontrols);
                for iNode = 1:N
                    u = X(iu1);
                    x1 = X(ix1);
                    D = [1	0	0	0	0	0; -1	0	0	0	0	0; -1	-1	0	0	0	0;
                            1	1	0	0	0	0; 0	1	0	0	0	0; 0	-1	-1	0	0	0;
                            0	0	-1	0	0	0; 0	0	1	0	0	0; 0	0	0	1	0	0;
                            0	0	0	-1	0	0; 0	0	0	-1	-1	0; 0	0	0	1	1	0;
                            0	0	0	0	1	0; 0	0	0	0	-1	-1;0	0	0	0	0	-1;
                            0	0	0	0	0	1];
                    %feedback gain is 10x feedback states
                    u1 = u+D*(x1(4:ndof)*X(end-1)*10+x1(ndof+(4:ndof))*X(end)*10);
                    
                    f2 = f2 + W * sum(u1.^expon)/N;
                    g(iu1) = g(iu1) + expon * W * u1.^(expon-1) / N;
                    
                    g(ix1) = g(ix1) + [zeros(nmus,3) D*10*X(end-1) zeros(nmus,3)  D*10*X(end)  zeros(nmus,nmus*2)]'*(expon * W * u1.^(expon-1) / N );
                    g(end-1) = g(end-1) +  (10*D*X(ix1(4:ndof)))'* expon * W * u1.^(expon-1) / N;
                    g(end) = g(end) + (10*D*X(ix1(ndof+(4:ndof))))'* expon * W * u1.^(expon-1) / N;
                
                    iu1 = iu1 + problem.nvarpernode1;
                    ix1 = ix1 + nvarpernode; %ix2;%
                    if i == N
                       nvarpernode = problem.nvarpernode;
                    end
                end
            end
        end
    end

	% cost of valve operation (mean of squared valve speed)
	f3 = 0;
    nvarpernode = problem.nvarpernode1;
    for j = 1:Ncycles
        iu1 = nstates + nmus + (1:2);		% valve controls in node 1
        for i=1:N
            if (i==N)
                iu2 = nstates + nmus + (1:2);
            else
                iu2 = iu1 + nvarpernode;	% valve controls in next node
            end
            v = (X(iu2) - X(iu1))/h(j);
            f3 = f3 + sum(v.^2)/(2*N);
            g(iu1) = g(iu1) - model.Wvalve * v/h(j)/(N); 
            g(iu2) = g(iu2) + model.Wvalve * v/h(j)/(N); 
            g(end) = g(end) - model.Wvalve * sum(v.^2)/(N)^2/h(j);
            iu1 = iu1 + nvarpernode;
        end
    end

	f3 = model.Wvalve*f3;

    % Regularization term
    f4 = 0;	
    nvarpernode = problem.nvarpernode1;
	if model.Wreg ~= 0
        factor = model.Wreg/N^3/(problem.nvarpernode1+problem.nvarpernode*(Ncycles-1));
        % state regularization
        ix = 1:nstates;			% index to states and controls of first node
        % mean square of 1st derivative
        for j = 1:Ncycles
            for i=1:N
                x1 = X(ix);
                ix2 = ix + nvarpernode;
                if i == N
                    if j == Ncycles %start from beginning
                        xstart = 0;
                    else
                        xstart = ix(1)+nvarpernode-1;
                    end
                    if problem.symmetry
                        ix2 = xstart+problem.vmx;     
                    else
                        ix2 = xstart+(1:nstates);
                    end
                end
                x2 = X(ix2);

                xd = x2 - x1;
                f4 = f4 + sum(xd.^2);			% divide by N^2 is done when calculating f4
                % x1 = X(ix);
                g(ix) = g(ix) - 2*factor * xd;			
                % x2 = X(ix+nvarpernode);
                g(ix2) = g(ix2) + 2*factor * xd;		
                
                ix = ix + nvarpernode;			% move pointers to next node
                if i == N
                   nvarpernode = problem.nvarpernode;
                end
            end
        end
        % inputs
        nvarpernode = problem.nvarpernode1;
        iu = nstates+(1:ncontrols);
        for i = 1:N %inputs only first gait cycle
            u1 = X(iu);
            iu2 = iu + nvarpernode;
            if i == N             
                if (problem.symmetry) %use inputs from opposite leg     
                    iu2 = nstates+problem.vmu;
                else % use 'normal order'
                    iu2 = nstates+(1:ncontrols);
                end
            end
            u2 = X(iu2);

            ud = u2 - u1;
            f4 = f4 + sum(ud.^2);			% divide by N^2 is done when calculating f4
            % x1 = X(ix);
            g(iu) = g(iu) - 2*factor * ud;			
            % x2 = X(ix+nvarpernode);
            g(iu2) = g(iu2) + 2*factor * ud;		
            iu = iu + nvarpernode;			% move pointers to next node
        end    
		f4 = f4*factor;								% make it the average
	end
    
	% add up the cost function components and store them all in a row
	f = f1 + f2 + f3 + f4;
	f = [f f1 f2 f3 f4];
	
	%-----------------------------------------------------------------------------
	% constraint violations c and their jacobian J
	c = zeros(ncon,1);
	J = spalloc(ncon, nvar, problem.Jnnz);

    nvarpernode = problem.nvarpernode1;
	% pointers to states and controls in node 1
	ix1 = 1:nstates;
	% pointers to constraints for node 1
	ic = 1:nstates;
    
	% evaluate dynamics at each pair of successive nodes
    for j = 1:Ncycles
        iu1 = nstates+(1:ncontrols);
        for i=1:N
            ix2 = ix1 + nvarpernode;
            iu2 = iu1 + problem.nvarpernode1;
            if i == N
                if j == Ncycles %start from beginning
                    xstart = 0;
                else
                    xstart = ix1(1)+nvarpernode-1;
                end
                if problem.symmetry
                    % use state of other leg
                    ix2 = xstart+problem.vmx;     
                    iu2 = nstates+problem.vmu;
                else
                    ix2 = xstart+(1:nstates);     
                    iu2 = nstates+(1:ncontrols);
                end
            end
            x1 = X(ix1);
            x2 = X(ix2);
            if and(i == N, j == Ncycles)
                if problem.symmetry
                    % Divide by two here, duration state for full gait cycle
                    x2(1) = x2(1) + speed * sum(X(end-1-Ncycles:end-2))/2;		% add horizontal translation for duration
                else
                    x2(1) = x2(1) + speed * sum(X(end-1-Ncycles:end-2));		% add horizontal translation for duration
                end
            end
            u1 = X(iu1);
            u2 = X(iu2);
            if strcmp(problem.model.type, 'torque')
                D = eye(6);
            else 
                D = [1	0	0	0	0	0; -1	0	0	0	0	0; -1	-1	0	0	0	0;
                    1	1	0	0	0	0; 0	1	0	0	0	0; 0	-1	-1	0	0	0;
                    0	0	-1	0	0	0; 0	0	1	0	0	0; 0	0	0	1	0	0;
                    0	0	0	-1	0	0; 0	0	0	-1	-1	0; 0	0	0	1	1	0;
                    0	0	0	0	1	0; 0	0	0	0	-1	-1;0	0	0	0	0	-1;
                    0	0	0	0	0	1];
            end
            %feedback gain is 10x feedback states
            u1 = u1+D*(x1(4:ndof)*X(end-1)*10+x1(ndof+(4:ndof))*X(end)*10);     % feedback isnot on first three states
            u2 = u2+D*(x2(4:ndof)*X(end-1)*10+x2(ndof+(4:ndof))*X(end)*10);

            % evaluate dynamics violation, and derivatives
            if (discr==1)						% midpoint
                [c(ic), dfdx, dfdxdot, dfdu] = dyn(model,(x1+x2)/2,(x2-x1)/h(j),(u1+u2)/2);	%easydyn((x1+x2)/2,(x2-x1)/h);%
                J(ic,ix1) = dfdx/2 - dfdxdot/h(j);
                J(ic,ix2) = dfdx/2 + dfdxdot/h(j);
%                 if strcmp(problem.model.type, 'torque')
                    J(ic,ix1) = J(ic,ix1) + [zeros(nstates,3) dfdu/2*D*X(end-1)*10 zeros(nstates,3) dfdu/2*D*X(end)*10 zeros(nstates,nmus*2)];
                    J(ic,ix2) = J(ic,ix2) + [zeros(nstates,3) dfdu/2*D*X(end-1)*10 zeros(nstates,3) dfdu/2*D*X(end)*10  zeros(nstates,nmus*2)];
%                 end
                J(ic,end-1) = 10*dfdu/2*D*(X(ix1(4:ndof))+X(ix2(4:ndof)));
                J(ic,end) = 10*dfdu/2*D*(X(ix1(ndof+(4:ndof)))+X(ix2(ndof+(4:ndof))));
                J(ic,iu1) = dfdu/2;
                J(ic,iu2) = dfdu/2;
            else								% backward euler
                [c(ic), dfdx, dfdxdot, dfdu] = dyn(model,x2,(x2-x1)/h(j),u2);	%easydyn(x2,(x2-x1)/h);%
                %add noise -> neural means to input, this is the same
%                 c(ic(ndof+4:ndof*2)) = c(ic(ndof+4:ndof*2))+problem.model.rndval((j-1)*N+i);
                % no noise in prosthesis, but extra in to simulate
                % imperfect connection
                c(ic([13 14 16 17])) = c(ic([13 14 16 17]))+problem.model.rndval((j-1)*N+i);
                
                J(ic,ix1) = -dfdxdot/h(j);
                J(ic,ix2) = dfdx + dfdxdot/h(j);
%                 if strcmp(problem.model.type, 'torque')
                    J(ic,ix2) = J(ic,ix2) + [zeros(nstates,3) dfdu*D*X(end-1)*10 zeros(nstates,3) dfdu*D*X(end)*10 zeros(nstates,nmus*2)];
%                 end
                J(ic,end-1) = 10*dfdu*D*x2(4:ndof);
                J(ic,end) = 10*dfdu*D*x2(ndof+(4:ndof));
                J(ic,iu2) = dfdu;
            end

            % and also derivatives of constraint w.r.t duration X(end)
            % df/ddur = df/dxdot * dxdot/dh * dh/ddur + (df/dx * dx/dx2 + df/dxdot * dxdot/dx2) * dx2/ddur
            %         = df/dxdot * (-(x2-x1)/h^2) * (1/N) + ...
            % the last term with dx2/ddur only exists at node N, and remember that f does not depend on x(1)!
            % also remember that df/dx is zero for the x that is horizontal position, on level ground
            J(ic,end-2-Ncycles+j) = -dfdxdot * (x2-x1) / h(j)^2 / (N);	
            
            if problem.symmetry
                %h = dur/N/2
                J(ic,end-2-Ncycles+j) = J(ic,end-2-Ncycles+j)/2;
            end
            % add df/dxdot * dxdot/dx2 * dx2/ddur
            if and(i==N,j==Ncycles)
                if problem.symmetry
                    %x2(1) = x2(1) + speed * sum(X(end-1-Ncycles:end-2))/2;
                    % Divide by two here, duration state for full gait cycle
                    J(ic,end-1-Ncycles:end-2) = J(ic,end-1-Ncycles:end-2) + repmat(dfdxdot(:,1) / h(j) * speed/2,1,Ncycles);
                else
                    J(ic,end-1-Ncycles:end-2) = J(ic,end-1-Ncycles:end-2) + repmat(dfdxdot(:,1) / h(j) * speed,1,Ncycles);
                end
                
            end
                
            % advance the indices
            ix1 = ix1 + nvarpernode; %ix2;%
            iu1 = iu1 + problem.nvarpernode1; %iu2;%
            ic = ic + nstates;		% there are nstates constraints for each node
            
            if i == N
               nvarpernode = problem.nvarpernode;
            end
        end 
    end

	c = problem.Wc * c;
	J = problem.Wc * J;
end
%=====================================================
function evaluate_if_needed(X)
	global optim problem
	% if X was not seen before, evaluate the objective and constraints
	if problem.print == 2 || ~isfield(optim,'X') || (min(X == optim.X) == 0)
		if problem.lambda == 1
			[optim.f, optim.g, optim.c, optim.J] = evaluate(problem.model, X);
		else
			[f1, g1, c1, J1] = evaluate(problem.model1, X);
			[f, g, c, J] = evaluate(problem.model, X);
			optim.f = (1-lambda)*f1 + lambda*f;
			optim.g = (1-lambda)*g1 + lambda*g;
			optim.c = (1-lambda)*c1 + lambda*c;
			optim.J = (1-lambda)*J1 + lambda*J;		
		end
		optim.X = X;

		% log the results of this evaluation
		normc = norm(optim.c(1:end-2*problem.N));
		row = [normc optim.f];
		row(find(row==0)) = NaN;
		problem.log = [problem.log ; row];

		% print something on screen, if it is requested
		if and(problem.print == 0,toc > problem.Printinterval)
            report(X);
            tic;
        elseif or(problem.print == 2,toc > problem.Printinterval)
			report(X);
			fprintf('%d -- Normc: %8.6f  ', size(problem.log,1), normc);
			fprintf('Obj: %8.5f = %8.5f (track) + %8.5f (effort) + %8.5f (valves)+ %8.5f (reg)\n', optim.f);
			savefile(X, 'tmpsave.mat');
			tic;
		end
	end
end
%=====================================================

function [fc, gc] = objconfun(X)

    global optim 
	
	evaluate_if_needed(X);
    fc = [optim.f(1);optim.c];
    gc = [transpose(optim.g);optim.J];
    
end

function c = confun_2d(X)

	global optim
	
	evaluate_if_needed(X);
	c = optim.c;
end
%=====================================================
function J = conjac_2d(X);

	global optim 
	
	evaluate_if_needed(X);
	J = optim.J;
end
%=====================================================
function J = conjacstructure_2d;

	global problem 
	
	% copy result from global struct
	J = problem.Jpattern;
end
%=====================================================
function f = objfun_2d(X);

	global optim 
	
	evaluate_if_needed(X);
	f = optim.f(1);
end
%=====================================================
function g = objgrad_2d(X);

	global optim 
	
	evaluate_if_needed(X);
	g = optim.g;
end
%================================================================================================
function savefile(X, filename, result_input)
	% save the result X in a result structure on a file
	global problem
	
	% first gait cycle
    x = reshape(X(1:problem.nvarpernode1*problem.N),problem.nvarpernode1, problem.N);
    % other gait cycles
    xj = reshape(X(problem.nvarpernode1*problem.N+1:end-2-problem.Ncycles),problem.nvarpernode,problem.N*(problem.Ncycles-1));
	if (nargin > 2)
		result = result_input;
	end
	result.x = [x(1:problem.nstates,:) xj];
	result.u = x(problem.nstates+(1:problem.ncontrols),:);
	result.dur = X(end-1-problem.Ncycles:end-2);
	result.speed = problem.model.data.speed;
	result.model = problem.model;
	result.normc = problem.log(end,1);
	result.f = problem.log(end,2:end);
	save(filename,'result');
end
%===========================================================================================
function [outdata] = makeperiodic(indata, terms);
	% makes movement data periodic by fitting a truncated Fourier series
	% data format: rows represent time, columns represent different variables

	spectrum = fft(indata);
	
	% only keep low frequency terms in spectrum
	n = size(spectrum,1);
	spectrum(2+terms:n-terms,:) = 0;
	
	outdata = ifft(spectrum);
end
%===========================================================================================
function drawstick(x);

	R = [1:6 4];			% right stick points
	L = [2 7:10 8];			% left stick points
	xrange = [min(x(1,:))-0.5  , max(x(1,:))+0.5];
	yrange = [min(x(2,:))-1.2  , max(x(2,:))+0.5];
	for i=1:size(x,2)
		plot(xrange,[0 0],'k');		% draw ground surface as a black line
		hold on
		d = gait2d('Stick',x(:,i));
		plot(d(R,1),d(R,2),d(L,1),d(L,2));
		axis([xrange yrange]);
	end
	hold off;
end
%=========================================================================================
function report(X)
	% displays a progress report on the screen, using current iterate X
	global problem

	% stick figure of current solution
	subplot(1,3,1)
    % first gait cycle
    x = reshape(X(1:problem.nvarpernode1*problem.N),problem.nvarpernode1, problem.N);
    % other gait cycles
    xj = reshape(X(problem.nvarpernode1*problem.N+1:end-2-problem.Ncycles),problem.nvarpernode,problem.N*(problem.Ncycles-1));
	nxg = 2*problem.ndof + 2*problem.nmus;
    x = [x(1:nxg,:) xj(1:nxg,:)]; % only musculoskeletal states
    if strcmp(problem.model.type,'torque')
        x(end+1:end+32,:) = 0;
    end
	drawstick(x);
	title([num2str(problem.N) ' nodes']);
	
	% optimization log
	if size(problem.log,1) > 1
		subplot(1,3,2);
		nlog = size(problem.log,1);
		plot(problem.log(:,[3 4 5 2]));
		% compute a suitable scale for the vertical axis
		maxy = mean(problem.log(:,2)) + std(problem.log(:,3));
		maxy = 2*max(problem.log(round(nlog/2):end,2));
		miny = min(min(problem.log(:,[3 4 5 2])));
		set(gca,'YLim',[miny maxy])
		legend('tracking','effort','valve','total','Location','Best');
		title('Objective');

		subplot(1,3,3);
		semilogy(problem.log(:,1)); 
		title('Norm of constraint violations');
		xlabel('Function evaluations');
	end
	pause(1e-6);

end
%==============================================================================================
function matcompare(a,b);
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('%8.5f at %d %d (%f should be %f)\n', maxerr, irow, icol, full(a(irow,icol)),full(b(irow,icol)));
end
%=================================================================
function [s] = ipoptmessage(info)

 	if info==0;  s = 'solved'; return; end;
 	if info==1;  s = 'solved to acceptable level'; return; end;
 	if info==2;  s = 'infeasible problem detected'; return; end;
 	if info==3;  s = 'search direction becomes too small'; return; end;
 	if info==4;  s = 'diverging iterates'; return; end;
 	if info==5;  s = 'user requested stop'; return; end;
%     
 	if info==-1;  s = 'maximum number of iterations exceeded'; return; end;
 	if info==-2;  s = 'restoration phase failed'; return; end;
 	if info==-3;  s = 'error in step computation'; return; end;
 	if info==-10;  s = 'not enough degrees of freedom'; return; end;
 	if info==-11;  s = 'invalid problem definition'; return; end;
 	if info==-12;  s = 'invalid option'; return; end;
 	if info==-13;  s = 'invalid number detected'; return; end;
%
 	if info==-100;  s = 'unrecoverable exception'; return; end;
 	if info==-101;  s = 'non-IPOPT exception thrown'; return; end;
 	if info==-102;  s = 'insufficient memory'; return; end;
 	if info==-199;  s = 'internal error'; return; end;

end

function [newproblem] = makemirror(problem)	
    if strcmp(problem.model.type, 'torque')
        problem.vmx = [1:3 7:9 4:6 10:12 16:18 13:15];
        problem.vmu = [4:6 1:3];
    elseif problem.proscontrol
        problem.vmx = [1:3 7:9 4:6 10:12 16:18 13:15 27:34 19:26 43:50 35:42 51:54];
    	problem.vmu = [9:16 1:8 17:18];
    else
        problem.vmx = [1:3 7:9 4:6 10:12 16:18 13:15 27:34 19:26 43:50 35:42];
    	problem.vmu = [9:16 1:8];
    end

	newproblem = problem;
	
end
