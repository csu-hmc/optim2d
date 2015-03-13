function simulate(duration,resultfile)
	% simulates the model with initial conditions and controls from resultfile
	
	% use global variable to store the controller parameters
	global controlpar

	% Load the file with simulation result
	if (nargin == 1)
		[file,parentfolder] = uigetfile('*.mat');
		resultfile = [parentfolder file];
	end
	result = load(resultfile);
	result = result.result;
	dur = result.dur;
	speed = result.speed;
	N = size(result.x,2);
	h = dur/N;				% time step for simulation
	
	% store the open loop control pattern in the controller parameters
	controlpar.t = (0:N-1)'*h;		% sample times
	controlpar.T = result.dur;		% repetition period
	controlpar.u = result.u';		% controls at sample times

	% Initialize the model
	model = initmodel(result.model);

	% Simulate the model with either implicit midpoint or implicit euler method
	nsteps = round(duration/h);
	xsim = zeros(size(result.x,1),nsteps+1);
	usim = zeros(size(result.u,1),nsteps+1);
	xsim(:,1) = result.x(:,1);				
	usim(:,1) = result.u(:,1);				
	euler = strcmp(model.discretization, 'euler');
	options.maxiterations = 100;
	options.ftol = 1e-6;
	t = 0;
	tic
	for i = 1:nsteps
		i1 = mod(i-1,N)+1;			% index to first node of time interval
		i2 = mod(i,N)+1;			% index to the next node
		
		% solve the time step
		x1 = xsim(:,i);
		xsim(:,i+1) = fsolve1(@simfun, x1, options);
		usim(:,i+1) = result.u(:,i2);
		t = t + h;
	end
	fprintf('Simulated %.3f seconds (%d integration steps) in %.3f seconds real time.\n', duration, nsteps, toc);
	
	% save the result
	result.x = xsim;
	result.u = usim;
	result.dur = duration;
	result.normc = NaN;
	result.f = NaN;
	result.solver = 'simulate';
	result.message = 'Forward simulation';
	outputfile = [parentfolder 'sim_' file];
	save(outputfile, 'result');
	fprintf('Result was saved on file:\n%s\n', outputfile);
	
	%================================================================
	function [f,J] = simfun(x);
		% function to solve x at the end of the time step using f(x) = 0
		if (euler)
			% Euler discretization uses controls at t+h, x(t+h)
			u = controller(t+h, x);
			[f, dfdx, dfdxdot] = dyn(model, x, (x-x1)/h, u);
			J = dfdx + dfdxdot/h;
		else
			% Midpoint discretization uses controls at t+h/2, (x(t) + x(t+h))/2
			u = controller(t+h/2, (x+x1)/2);
			[f, dfdx, dfdxdot] = dyn(model, (x+x1)/2, (x-x1)/h, u);
			J = dfdx/2 + dfdxdot/h;		
		end
	end
	%==================================================================
end
%===============================================================================
function [u] = controller(t,x)
	% here is where we put our controller (open loop or closed loop)
	
	global controlpar

	% present version only has open loop (u depends only on time t, not on state x)
	
	% open loop controls are periodic with period T, so map t into the range of zero to T
	tper = mod(t, controlpar.T);
	
	% interpolate the stored control patterns, to obtain controls at any time t
	u = interp1(controlpar.t, controlpar.u, tper, 'linear', 'extrap')';
	
	% add your feedback corrections here...
	
end
%===============================================================================
function [x,f,info] = fsolve1(fun, x0, options)
	% solves f(x)=0 using Newton's method

	% initialize
	info = 0;
	iterations = 0;
	x = x0;
	[f, J] = fun(x);
	bestfnorm = norm(f);
		
	% Start loop
	while (1)
		iterations = iterations + 1;
		dx = -J\f;						% Newton direction
		
		% do the step dx
		x = x + dx;
		
		% check if number of iterations exceeded
		if (iterations > options.maxiterations)
			info = 1;
			return
		end
		
		% Make sure we have improvement in the merit function norm(f)
		[f,J] = fun(x);
		fnorm = norm(f);
		while (fnorm > bestfnorm)
			dx = dx/2;
			x = x-dx;
			[f,J] = fun(x);
			fnorm = norm(f);
		end

		% are we done?
		if fnorm < options.ftol
			info = 0;
			return;
		end

	end
	
end
