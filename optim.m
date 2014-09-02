function [result] = optim(prob);

% predictive optimization

	tic;
	clear global
	clear mex
	global problem
	
	problem = prob;
	problem.print = 0;						% 0: never, 1: at Printinterval, 2: always
		
	%----------------------------------------------------------------------------------
	% Number of nodes
	N			= problem.N;				% number of collocation nodes
	if mod(N,2)~=0
		error('Number of nodes must be even.');
	end
	
	%----------------------------------------------------------------------------------
	% Some model related constants
	ndof = 9;
	nmus = 16;
	nstates = 2*ndof + 2*nmus + 4;
	ncontrols = nmus + 2;
	nvarpernode = nstates + ncontrols;		% number of unknowns per time node
	nvar = nvarpernode * N + 1;				% total number of unknowns (states, controls, duration)
	ncon = nstates * N;						% number of constraints due to discretized dynamics and task

	problem.ndof = ndof;
	problem.nmus = nmus;
	problem.nstates = nstates;
	problem.ncontrols = ncontrols;
	problem.nvarpernode = nvarpernode;
	problem.nvar = nvar;
	problem.ncon = ncon;
	
	% precompute the indices for muscle controls
	iu = [];
	for i=0:N-1
		iu = [iu nvarpernode*i+nstates+(1:nmus)];
	end
	problem.iumus = iu;
	
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
	Lu = zeros(nmus,1);
	Uu = 2*ones(nmus,1);
	% bounds for musculoskeletal states
	Lx = [-1; 0.7; -pi/4; -pi; -pi; -pi; -pi; -pi; -pi; 		... % q
			-5; -5;  -20; -20; -20; -20; -20; -20; -20;  	... % qdot
			zeros(16,1)-1;									...	% Lce
			Lu]; 											...	% active states
	Ux = [4; 1.3; pi/4; pi; pi; pi; pi; pi; pi; 				... % q
			5; 5; 20; 20; 20; 20; 20; 20; 20;  			... % qdot
			zeros(16,1)+5; 									...	% Lce
			Uu];											... % active states
	% bounds for prosthesis states
	Lxp = [-100 ; -100; -100; -500];			% s, v1, v2, M
	Uxp = [100 ; 100; 100; 500];
	% bounds for prosthesis controls
	Lup = zeros(2,1)+0.001;
	Uup = ones(2,1);
	for k = 0:N-1
		L(k*nvarpernode + (1:nvarpernode) ) = [Lx; Lxp; Lu; Lup];
		U(k*nvarpernode + (1:nvarpernode) ) = [Ux; Uxp; Uu; Uup];
	end
	L(end) = 0.2;		% minimum duration of movement cycle
	U(end) = 2.0;		% maximum duration of movement cycle
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
		X0 = (L + U)/2;								% halfway between upper and lower bound 
		% X0 = X0 + 0.01*(U-L).*randn(size(X0));	% to add some randomness
		% X0 = L + (U-L).*rand(size(L));			% to make the initial guess completely random
	else
		% load a previous solution
		load(problem.initialguess);
		Nresult = size(result.x,2);
		t0 = (0:(Nresult-1))'/Nresult;
		x0 = result.x';
		u0 = result.u';
		
		% initial guess for valve controls is random (to avoid getting stuck in local optimum)
		% u0(:,17:18) = rand(size(u0(:,17:18)));
		% u0(:,17:18) = ones(size(u0(:,17:18)));
		
		% duplicate first node node at the end so we can interpolate with periodicity
		t0 = [t0 ; 1];
		u0 = [u0 ; u0(1,:)];
		x0 = [x0 ; x0(1,:)];
		x0(end,1) = x0(end,1) + result.speed * result.dur;
		
		% interpolate states and controls from initial guess to current node times
		times = (0:(N-1))/N;
		x0 = interp1(t0,x0,times,'linear','extrap');
		u0 = interp1(t0,u0,times,'linear','extrap');
		X0 = reshape([x0 u0]',nvar-1,1);
		X0 = [X0 ; result.dur];
	end
	problem.X0 = X0;			% store initial guess in case we need it later
	
	% determine sparsity structure of Jacobian
	% we have verified that we always get same structure by testing with random X, so do this only once
	for i=1:1
		problem.Jnnz = 1;
		X = L + (U-L).*rand(size(L));		% a random vector of unknowns
		J = conjac(X);
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
		keyboard
	
		% check the NLP derivatives
		hh = 1e-7;
		X = L + (U-L).*rand(size(L));		% a random vector of unknowns
		f = objfun(X);
		grad = objgrad(X);
		c = confun(X);
		cjac = conjac(X);
		cjac_num = zeros(ncon,nvar);
		grad_num = zeros(nvar,1);
		for i=1:nvar
			fprintf('checking objgrad and conjac for unknown %4d of %4d\n',i,nvar);
			Xisave = X(i);
			X(i) = X(i) + hh;
			cjac_num(:,i) = (confun(X) - c)/hh;
			grad_num(i) =   (objfun(X) - f)/hh;
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
	problem.print = 2;
	evaluate_if_needed(X0);
	if (problem.debug)
		disp('Hit ENTER to start optimization...');pause
	end

	%----------------------------------------------------------------------------------------
	% run optimization
	result.solver = problem.Solver;
	if (problem.MaxIterations > 0)	
		problem.print = 1;
		fprintf('Starting optimization...\n');
		if strcmp(problem.Solver,'IPOPT')
			funcs.objective = @objfun;
			funcs.gradient  = @objgrad;
			funcs.constraints = @confun;
			funcs.jacobian    = @conjac;
			funcs.jacobianstructure = @conjacstructure;
			options.lb = L;
			options.ub = U;
			options.cl = zeros(ncon,1);
			options.cu = zeros(ncon,1);	
			options.ipopt.max_iter = problem.MaxIterations;
			options.ipopt.hessian_approximation = 'limited-memory';
			% options.ipopt.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
			% options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
			% options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
			% options.ipopt.bound_push = options.ipopt.bound_frac;
			options.ipopt.tol = problem.Tol;
			options.ipopt.constr_viol_tol = problem.ConstraintTol;
			options.ipopt.compl_inf_tol = problem.ConstraintTol;
			% options.ipopt.linear_solver = 'ma57';        % default is ma27, options: ma57, mumps (requires library install)
			options.ipopt.print_level = 0;
			[X, info] = ipopt(X0,funcs,options);
			result.info = info.status;
			result.message = ipoptmessage(info.status);
		elseif strcmp(problem.Solver,'SNOPT')	
			Prob = conAssign(@objfun, @objgrad, [], [], L, U, 'optim2d', X0, ...
                            [], 0, ...
                            [], [], [], @confun, @conjac, [], problem.Jpattern, ...
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
	disp(['Solver completion status: ' result.message]);
	fprintf('Duration of movement: %8.4f s\n',X(end));
	
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
		disp('Warning: gait cycle duration had small SD.  Using SD = 5 ms instead.');
		gait.dur(2) = 0.005;
	end

	% store average and SD in N x 10 matrix
	av = [av(1:N,:) av([N/2+1:N 1:N/2],:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
	sd = [sd(1:N,:) sd([N/2+1:N 1:N/2],:)];		% right side data in columns 1-5, left side in columns 6-10, phase shifted	
	
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
	nstates = problem.nstates;
	ncontrols = problem.ncontrols;
	nvarpernode = problem.nvarpernode;
	nmus = problem.nmus;
	ncon = problem.ncon;
	nvar = problem.nvar;
	nxg = 2*problem.ndof + 2*nmus;		% number of musculoskeletal states
	iumus = problem.iumus;
	
	BW = 9.81*model.mass;
	gait2d = model.gait2d;
	
	h = X(end)/N;		% time step size
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
	for i=1:N
		x = X(ixg);
		[GRF, dGRFdx] = gait2d('GRF',x);
		GRF = GRF/BW;							% convert to units of body weight
		dGRFdx = dGRFdx/BW;						% here also
		simdata = [x(4:6); GRF(1:2); x(7:9); GRF(3:4)]';	% the ten tracked variables: right angles, right GRF, left angles, left GRF
		res = (simdata-av(i,:))./sd(i,:);		% the ten residuals, normalized to SD
		f1 = f1 + sum(res.^2)/(11*N);
		res = res./sd(i,:);						% divide again by SD, needed for gradient
		g(ixg(4:6)) = g(ixg(4:6)) + 2*res(1:3)'/(11*N);		% gradient of objective with respect to right side angles
		g(ixg(7:9)) = g(ixg(7:9)) + 2*res(6:8)'/(11*N);		% gradient of objective with respect to left side angles
		g(ixg) = g(ixg) + 2*dGRFdx'*res([4 5 9 10])'/(11*N);	% gradient of GRF terms in objective, with respect to all state variables
		ixg = ixg + nvarpernode;								% move pointers to next node
	end
	
	% add tracking error for movement duration
	f1 = f1 + ((X(end) - dur)/dursd)^2/11;				
	g(end) = g(end) + 2*(X(end) - dur)/dursd^2/11;		
	
	% apply weighting to the tracking term
	f1 = f1*model.Wtrack;
	g = g*model.Wtrack;
	
	% effort term for muscles
	f2 = 0;
	expon = abs(model.effort.exponent);
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

	% cost of valve operation (mean of squared valve speed)
	iu1 = nstates + nmus + (1:2);		% valve controls in node 1
	f3 = 0;
	for i=1:N
		if (i==N)
			iu2 = nstates + nmus + (1:2);
		else
			iu2 = iu1 + nvarpernode;	% valve controls in next node
		end
		v = (X(iu2) - X(iu1))/h;
		f3 = f3 + sum(v.^2)/(2*N);
		g(iu1) = g(iu1) - model.Wvalve * v/h/N; 
		g(iu2) = g(iu2) + model.Wvalve * v/h/N; 
		g(end) = g(end) - model.Wvalve * sum(v.^2)/N^2/h;
		iu1 = iu1 + nvarpernode;
	end
	f3 = model.Wvalve*f3;

	% add up the cost function components and store them all in a row
	f = f1 + f2 + f3;
	f = [f f1 f2 f3];
	
	%-----------------------------------------------------------------------------
	% constraint violations c and their jacobian J
	c = zeros(ncon,1);
	J = spalloc(ncon, nvar, problem.Jnnz);

	% pointers to states and controls in node 1
	ix1 = 1:nstates;
	iu1 = nstates+(1:ncontrols);
	% pointers to constraints for node 1
	ic = 1:nstates;
	% evaluate dynamics at each pair of successive nodes
	for i=1:N
		if (i < N)
			ix2 = ix1 + nvarpernode;
			iu2 = iu1 + nvarpernode;
		else
			ix2 = 1:nstates;							% use node 1
			iu2 = nstates+(1:ncontrols);				% use node 1
		end
		x1 = X(ix1);
		x2 = X(ix2);
		if (i == N)
			x2(1) = x2(1) + speed * X(end);		% add horizontal translation for duration
		end
		u1 = X(iu1);
		u2 = X(iu2);
		
		% evaluate dynamics violation, and derivatives
		if (discr==1)						% midpoint
			[c(ic), dfdx, dfdxdot, dfdu] = dyn(model,(x1+x2)/2,(x2-x1)/h,(u1+u2)/2);	
			J(ic,ix1) = dfdx/2 - dfdxdot/h;
			J(ic,ix2) = dfdx/2 + dfdxdot/h;
			J(ic,iu1) = dfdu/2;
			J(ic,iu2) = dfdu/2;
		else								% euler
			[c(ic), dfdx, dfdxdot, dfdu] = dyn(model,x2,(x2-x1)/h,u2);	
			J(ic,ix1) = -dfdxdot/h;
			J(ic,ix2) = dfdx + dfdxdot/h;
			J(ic,iu2) = dfdu;
		end
		
		% and also derivatives of constraint w.r.t duration X(end)
		% df/ddur = df/dxdot * dxdot/dh * dh/ddur + (df/dx * dx/dx2 + df/dxdot * dxdot/dx2) * dx2/ddur
		%         = df/dxdot * (-(x2-x1)/h^2) * (1/N) + ...
		% the last term with dx2/ddur only exists at node N, and remember that f does not depend on x(1)!
		% also remember that df/dx is zero for the x that is horizontal position, on level ground
		J(ic,end) = -dfdxdot * (x2-x1) / h^2 / N;	
		% add df/dxdot * dxdot/dx2 * dx2/ddur
		if (i==N)
			J(ic,end) = J(ic,end) + dfdxdot(:,1) / h * speed;
		end
				
		% advance the indices
		ix1 = ix1 + nvarpernode;
		iu1 = iu1 + nvarpernode;
		ic = ic + nstates;		% there are nstates constraints for each node
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
		normc = norm(optim.c);
		row = [normc optim.f];
		row(find(row==0)) = NaN;
		problem.log = [problem.log ; row];

		% print something on screen, if it is requested
		if problem.print == 0
			return
		end
		if problem.print == 2 || toc > problem.Printinterval
			report(X);
			fprintf('%d -- Normc: %8.6f  ', size(problem.log,1), normc);
			fprintf('Obj: %8.5f = %8.5f (track) + %8.5f (effort) + %8.5f (valves)\n', optim.f);
			savefile(X, 'tmpsave.mat');
			tic;
		end
	end
end
%=====================================================
function c = confun(X);

	global optim
	
	evaluate_if_needed(X);
	c = optim.c;
end
%=====================================================
function J = conjac(X);

	global optim 
	
	evaluate_if_needed(X);
	J = optim.J;
end
%=====================================================
function J = conjacstructure;

	global problem 
	
	% copy result from global struct
	J = problem.Jpattern;
end
%=====================================================
function f = objfun(X);

	global optim 
	
	evaluate_if_needed(X);
	f = optim.f(1);
end
%=====================================================
function g = objgrad(X);

	global optim 
	
	evaluate_if_needed(X);
	g = optim.g;
end
%================================================================================================
function savefile(X, filename, result_input)
	% save the result X in a result structure on a file
	global problem
	
	x = reshape(X(1:end-1), problem.nvarpernode, problem.N);
	if (nargin > 2)
		result = result_input;
	end
	result.x = x(1:problem.nstates,:);
	result.u = x(problem.nstates+(1:problem.ncontrols),:);
	result.dur = X(end);
	result.speed = problem.model.data.speed;
	result.model = problem.model;
	result.normc = problem.log(end,1);
	result.f = problem.log(end,2);
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
		axis('equal');
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
	x = reshape(X(1:end-1), problem.nvarpernode, problem.N);
	nxg = 2*problem.ndof + 2*problem.nmus;
	x = x(1:nxg,:);			% use only the musculoskeletal states
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
	pause(0.1);

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
