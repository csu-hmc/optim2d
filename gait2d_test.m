function [result] = gait2d_test(command)

% This program runs various tests on the gait2d model

	close all

	if (nargin ~= 1)
		error('gait2d_test needs one argument (command string)');
	end
		
	% Some model related variables
	ndof = 9;
	nmus = 16;
	nstates = 2*ndof + 2*nmus;
	joints = {'Hip' 'Knee' 'Ankle' 'LHip' 'LKnee' 'LAnkle'};
	njoints = size(joints,2);

	% make sure we are using latest version of MEX function
	
	% Initialize the model
	warning('OFF','MATLAB:xlsread:Mode');
	parameters = xlsread('gait2d_par.xls','','','basic');
	warning('ON','MATLAB:xlsread:Mode');
	gait2d('Initialize', parameters);
	
	% to test the amputee options
	% gait2d('Set','Right AKA', 450.0);				% 450 Nm/rad based on Hansen et al., J Biomech 37: 1467–1474, 2004.
	
	% Extract Lceopt values, we need those in the isokinetic simulations
	Lceopt = gait2d('Get', 'Lceopt');
	
	dofnames = {'trunk x','trunk y','trunk angle', ...
			'Rhip angle','Rknee angle','Rankle angle', ...
			'Lhip angle','Lknee angle','Lankle angle'};

	% construct a free fall state
	xff = zeros(nstates,1);
	xff(2) = 1.2;					% put model some distance above origin
	xff(4) = 15*pi/180;				% right hip 5 deg flexed
	xff(5) = -15*pi/180;			% right knee 5 deg flexed
	xff(7) = -15*pi/180;			% left hip 5 deg extended
	xff(8) = -25*pi/180;			% left knee 5 deg flexed
	xff(19:34) = 1.4;				% Lce long enough to get small muscle force
	
	% xff = xff + 0.01*randn(size(xff));

	% ======================= do the stickfigure test
	if strcmp(command, 'stick')
		disp('Stick figure test...');
		clf;
		gait2dstick(xff);
	end
	
	% ======================= do the speed test
	if strcmp(command, 'speed')
		disp('Speed test...');
		tic
		Neval = 10000;
		for i=1:Neval
			x = rand(nstates,1);
			xdot = rand(nstates,1);
			u = rand(nmus,1);
			[f, dfdx, dfdxdot, dfdu] = gait2d('Dynamics',x,xdot,u);
		end
		fprintf('Computation time for each model evaluation: %8.5f ms\n',1000*toc/Neval);		
		tic
		for i=1:Neval
			x = rand(nstates,1);
			gait2d('Jointmoments',x);
		end
		fprintf('Computation time without Autolev code:      %8.5f ms\n',1000*toc/Neval);		
	end
	
	% ======================= do the GRF test
	if strcmp(command, 'grf')
		disp('GRF test...');
		x = xff;
		x(1:9) = zeros(9,1);		% put model in upright position
		x([6 9]) = 0.3;				% tilt the feet so only heel is on ground
		y = 0.945:0.0001:1.0;
		ny = size(y,2);
		vy = [-1 -.5 0 .5 1];
		nvy = size(vy,2);
		Fy = zeros(size(y,2),size(vy,2));
		for i = 1:nvy
			x(11) = vy(i);				% vertical speed
			for j = 1:ny
				x(2) = y(j);			% vertical position
				GRF = gait2d('GRF',x);
				Fy(j,i) = GRF(2);		% vertical GRF of right foot (heel)
			end
		end
		subplot(1,2,1);
		plot(y,Fy);
		legend([repmat('speed ',nvy,1) num2str(vy')]);
		xlabel('hip height (m)');
		ylabel('right VGRF (N)');

		% horizontal force-velocity curves
		x(11) = 0;		% reset vertical velocity to zero
		vx = -.2:0.001:.2;
		nvx = size(vx,2);
		y = 0.945:0.008:0.961;
		ny = size(y,2);
		Fx = zeros(nvx,ny);
		for i = 1:ny
			x(2) = y(i);				% vertical position
			for j = 1:nvx
				x(10) = vx(j);			% horizontal speed
				GRF = gait2d('GRF',x);
				Fx(j,i) = GRF(1);		% horizontal GRF of right foot (heel)
			end
		end
		subplot(1,2,2);
		plot(vx,Fx);
		legend([repmat('hip height ',ny,1) num2str(y')]);
		xlabel('horizontal speed (m/s)');
		ylabel('Fx (N)');
	end
	
	% ======================= do the free fall dynamics test
	if strcmp(command, 'freefall')	
		disp('Freefall dynamics test...');
		
		% we use the freefall state x constructed above
		% construct state derivatives that should be nearly correct for this state
		xdot = zeros(nstates,1);
		xdot(ndof+2) = -9.81;			% free fall acceleration, f should now be close to zero
		
		% test the dynamics
		u = zeros(nmus,1);
		[f] = gait2d('Dynamics',xff,xdot,u);
		
		% print the residuals
		fprintf('--------------------------------------------------\n');
		fprintf('velocity    multibody    contraction    activation\n');
		fprintf('--------------------------------------------------\n');
		for i=1:16
			if (i<10)
				fprintf('%9.4f %11.4f',f(i),f(9+i));
			else
				fprintf('                     ');
			end
			fprintf('%13.4f %13.4f\n',f(18+i),f(34+i));
		end
		fprintf('--------------------------------------------------\n');
		
		% get the GRF
		GRF = gait2d('GRF',xff);
		fprintf('GRF_heel = %8.3f %8.3f\nGRF_toe  = %8.3f %8.3f\n', GRF);

		fprintf('\nCheck that dynamics residuals in the table, and GRF, are all small.\n');
	end
	
	%======================== do the accelerometer test
	if strcmp(command, 'accelerometer')
		disp('Testing accelerometer model...');
		
		% orientation test
		angles = (0:pi/20:2*pi)';
		aa = [];
		for i = 1:numel(angles)
			q = [0 0 angles(i) 0 0 0 0 0 0]';
			qd = zeros(ndof,1);
			qdd = zeros(ndof,1);
			a = gait2d('Accelerometer',q,qd,qdd);
			aa = [aa ; a'];
		end
		accplot(angles,aa,'orientation (deg)');
		
		% horizontal acceleration test
		accelerations = -10:10;
		aa = [];
		for i=1:numel(accelerations)
			q = zeros(ndof,1);
			qd = zeros(ndof,1);
			qdd = zeros(ndof,1);
			qdd(1) = accelerations(i);
			a = gait2d('Accelerometer',q,qd,qdd);
			aa = [aa ; a'];
		end
		accplot(accelerations,aa,'hor acc (m s^-2)');

		% vertical acceleration test
		accelerations = -10:10;
		aa = [];
		for i = 1:numel(accelerations)
			q = zeros(ndof,1);
			qd = zeros(ndof,1);
			qdd = zeros(ndof,1);
			qdd(2) = accelerations(i);
			a = gait2d('Accelerometer',q,qd,qdd);
			aa = [aa ; a'];
		end
		accplot(accelerations,aa,'vert acc (m s^-2)');

		% centrifugal test
		angvels = -10:10;
		aa = [];
		for i = 1:numel(angvels)
			q = zeros(ndof,1);
			qd = zeros(ndof,1);
			qdd = zeros(ndof,1);
			qd(3) = angvels(i);
			a = gait2d('Accelerometer',q,qd,qdd);
			aa = [aa ; a'];
		end
		accplot(accelerations,aa,'ang vel (rad s^-1)');
		
	end
	
	% ======================= do the derivatives test
	if strcmp(command, 'derivatives')
		disp('Checking derivatives...');
		for i=1:5					% edit this line to do multiple tests
			x = rand(nstates,1);
			xdot = .1*randn(nstates,1);
			u = rand(nmus,1);
			M = rand(6,1);
			[f,dfdx,dfdxdot,dfdu,dfdM] = gait2d('Dynamics',x,xdot,u,M);
			[GRF, dGRFdx] = gait2d('GRF',x);
			for j=1:nstates
				% compute df/dx
				h = 1e-7;
				saved = x(j);
				x(j) = x(j)+h;
				[fnew] = gait2d('Dynamics',x,xdot,u,M);
				[GRFnew] = gait2d('GRF',x);
				dfdx_num(:,j) = (fnew-f)/h;
				dGRFdx_num(:,j) = (GRFnew-GRF)/h;
				x(j) = saved;
				
				% compute df/dxdot
				saved = xdot(j);
				if j>ndof && j<=2*ndof
					h = 1e10;		% for mass matrix, use large finite difference (it is not dependent on xdot!) to avoid roundoff error
				else
					h = 1e-5;		% for the other elements, this seems to be a good finite difference
				end
				xdot(j) = xdot(j)+h;
				[fnew] = gait2d('Dynamics',x,xdot,u,M);
				dfdxdot_num(:,j) = (fnew-f)/h;				
				xdot(j) = saved;
			end
			
			% compute df/du
			h = 1e-7;
			for j=1:nmus
				saved = u(j);
				u(j) = u(j)+h;
				[fnew] = gait2d('Dynamics',x,xdot,u,M);
				dfdu_num(:,j) = (fnew-f)/h;
				u(j) = saved;
			end
			
			% compute df/dM
			h = 1e7;		  % large finite difference because we know f is linear in M
			for j=1:6
				saved = M(j);
				M(j) = M(j)+h;
				[fnew] = gait2d('Dynamics',x,xdot,u,M);
				dfdM_num(:,j) = (fnew-f)/h;
				M(j) = saved;
			end
			
			% check the accelerometer derivatives
			q = rand(ndof,1);
			qd = rand(ndof,1);
			qdd = rand(ndof,1);
			[a,dadq,dadqd,dadqdd] = gait2d('Accelerometer',q,qd,qdd);
			h = 1e-7;
			for j=1:ndof
				saved = q(j);
				q(j) = q(j)+h;
				[anew] = gait2d('Accelerometer',q,qd,qdd);
				dadq_num(:,j) = (anew-a)/h;
				q(j) = saved;
				
				saved = qd(j);
				qd(j) = qd(j)+h;
				[anew] = gait2d('Accelerometer',q,qd,qdd);
				dadqd_num(:,j) = (anew-a)/h;
				qd(j) = saved;
				
				saved = qdd(j);
				qdd(j) = qdd(j)+h;
				[anew] = gait2d('Accelerometer',q,qd,qdd);
				dadqdd_num(:,j) = (anew-a)/h;
				qdd(j) = saved;
			end
			
			% report the differences
			fprintf('Checking df/dx...\n'); 		matcompare(dfdx, dfdx_num);
			fprintf('Checking df/dxdot...\n'); 		matcompare(dfdxdot, dfdxdot_num);
			fprintf('Checking df/du...\n'); 		matcompare(dfdu, dfdu_num);
			fprintf('Checking df/dM...\n'); 		matcompare(dfdM, dfdM_num);
			fprintf('Checking dGRF/dx...\n'); 		matcompare(dGRFdx, dGRFdx_num);
			fprintf('Checking da/dq...\n'); 		matcompare(dadq, dadq_num);
			fprintf('Checking da/dqd...\n'); 		matcompare(dadqd, dadqd_num);
			fprintf('Checking da/dqdd...\n'); 		matcompare(dadqdd, dadqdd_num);
			fprintf('===== finished with test %d\n',i);
			fprintf('Hit ENTER to continue.\n');
			pause
		end
	end
	
	% ======================= do the isometric muscle tests
	if strcmp(command, 'isometric')
		disp('Calculating isometric strength curves...');
		clf;
	   
		% hip moment-angle curves, at 30-deg intervals in knee angle
		isometric_curves(1,[-30:5:90],2,[-120:30:0]);
		% knee moment-angle curves, at 30-deg intervals in hip angle
		isometric_curves(2,[-120:5:0],1,[-30:30:90]);
		% ankle moment-angle curves, at 30-deg intervals in knee angle
		isometric_curves(3,[-60:5:60],2,[-120:30:0]);   
	end
	
	% ======================= do the isokinetic muscle tests
	if strcmp(command, 'isokinetic')
		disp('Calculating isokinetic strength curves...');
		clf;
	   
		% moment-angular velocity curves
		isokinetic_curves(1,[-1000:25:1000]);
		isokinetic_curves(2,[-1000:25:1000]);
		isokinetic_curves(3,[-1000:25:1000]);
	end
		
	% ====================== simulate freefall movement, first passive then active ==============
	if strcmp(command, 'simulate')
		disp('Simulating freefall movements...');
		clf;
				
		% set simulation parameters
		tend = 0.5;
		h = 0.001;
		nsteps = round(tend/h);
		options.tol = 1e-6;
		options.maxiterations = 100; 
		
		% run passive simulation
		disp('Running passive freefall simulation...');
		eval1 = gait2d('Get','Evaluations');
		[tout, xout, info] = IMstep(@nostim, xff,[0 tend],nsteps,options);
		eval2 = gait2d('Get','Evaluations');
		fprintf('Number of function evaluations: %d\n', eval2-eval1);
		figure(1);clf;
		gait2dstick(xout);
		title('Passive simulation');
		if (info < 0)
			error('IMstep did not converge');
		end

		% run simulation with 50% muscle stim
		disp('Running active freefall simulation...');
		[tout, xout, info] = IMstep(@halfstim, xff,[0 tend],nsteps,options);	
		eval3 = gait2d('Get','Evaluations');
		fprintf('Number of function evaluations: %d\n', eval3-eval2);
		figure(2);clf;
		gait2dstick(xout);
		title('50% muscle excitation');
		if (info < 0)
			error('IMstep did not converge');
		end
				
		% run the same simulation with 1/10th step size
		disp('Running active freefall simulation with smaller steps...');
		[tout3, xout3, info] = IMstep(@halfstim, xff,[0 tend],10*nsteps,options);	
				
		% run the same simulation with explicit ODE solver
		disp('Running active freefall simulation with ODE23...');
		[tout2, xout2] = ode23(@odefun, [0 tend], xff);
		xout2 = xout2';
		eval4 = gait2d('Get','Evaluations');
		fprintf('Number of integration steps:    %d\n', size(tout2,1));
		fprintf('Number of function evaluations: %d\n', eval4-eval3);
		figure(3);
		gait2dstick(xout2)
		title('50% muscle stimulation - ODE23');
		
		% compare the implicit solution to ODE23
		figure(4)
		for i=1:9
			subplot(3,3,i)
			plot(tout,xout(i,:)',tout3,xout3(i,:)', tout2,xout2(i,:));
			set(gca,'XLim',[0 tend])
			title(dofnames(i));
		end
		legend(['h = ' num2str(tend/nsteps)],['h = ' num2str(tend/nsteps/10)],'ODE23');
					
	end
	
	%===================== Start of embedded functions ============================================
	
	% nostim: muscle stimulation function for passive simulation
	function [u] = nostim(t);
		u = zeros(nmus,1);
	end
	
	% halfstim: muscle stimulation function for simulation with 50% muscle stim
	function [u] = halfstim(t);
		u = 0.5+zeros(nmus,1);
	end
	
	%=============================================================================================================
	function isometric_curves(joint1, range1, joint2, range2)
		% produces isometric moment-angle curves for a joint, one for each value in a range of angles in second joint
		% third joint angle is kept zero
		
		angles = zeros(njoints,1);
		angvel = zeros(njoints,1);
		pascurves = [];
		poscurves = [];
		negcurves = [];
		legends = {};
		for angle2 = range2
			angles(joint2) = angle2;
			pasmoments = [];
			posmoments = [];
			negmoments = [];
			for angle1 = range1
				angles(joint1) = angle1;
				pasmoments = [pasmoments maxmoment(joint1, angles, angvel, 0)];
				posmoments = [posmoments maxmoment(joint1, angles, angvel, 1)];
				negmoments = [negmoments maxmoment(joint1, angles, angvel, -1)];
			end
			pascurves = [pascurves  pasmoments'];
			poscurves = [poscurves  posmoments'];
			negcurves = [negcurves  negmoments'];
			legends = [legends ; [char(joints(joint2)) ' ' num2str(angle2)] ];
		end
		
		% plot total moments on left side of figure
		subplot(3,3,3*joint1-2);
		plot(range1, poscurves);hold on;
		plot(range1, negcurves);
		labels;
		title([char(joints(joint1)) ': total moment']);
		
		% plot passive moments in middle column of figure
		subplot(3,3,3*joint1-1);
		plot(range1, pascurves);hold on;
		labels;
		title([char(joints(joint1)) ': passive moment']);
		
		% subtract passive moments and plot in rightmost column of figure
		subplot(3,3,3*joint1);
		plot(range1, poscurves-pascurves);hold on;
		plot(range1, negcurves-pascurves);
		labels;
		title([char(joints(joint1)) ': active = total -- passive']);
		
		legend(legends);
		
		function labels
			a = get(gca);
			axis([min(range1) max(range1) a.YLim]);
			plot([0 0],a.YLim,'k:');
			plot(a.XLim,[0 0],'k:');
			hold off;
			xlabel('angle (deg)');
			ylabel('moment (Nm)');
		end
	end

	%=============================================================================================================
	function isokinetic_curves(joint, range)
		% produces isokinetic moment-angular velocity curves for a joint
		
		angles = zeros(njoints,1);
		angles([2 5]) = -5;			% knees 5 deg flexed
		angvel = zeros(njoints,1);
		pasmoments = [];
		posmoments = [];
		negmoments = [];
		for vel = range
			angvel(joint) = vel;
			pasmoments = [pasmoments maxmoment(joint, angles, angvel, 0)];
			posmoments = [posmoments maxmoment(joint, angles, angvel, 1)];
			negmoments = [negmoments maxmoment(joint, angles, angvel, -1)];
		end
		
		% plot total moments on left side of figure
		subplot(3,3,3*joint-2);
		plot(range, posmoments);hold on;
		plot(range, negmoments);
		labels;
		title([char(joints(joint)) ': total moment']);
		
		% plot passive moments in middle column of figure
		subplot(3,3,3*joint-1);
		plot(range, pasmoments);hold on;
		labels;
		title([char(joints(joint)) ': passive moment']);
		
		% subtract passive moments and plot in rightmost column of figure
		subplot(3,3,3*joint);
		plot(range, posmoments-pasmoments);hold on;
		plot(range, negmoments-pasmoments);
		labels;
		title([char(joints(joint)) ': active = total -- passive']);
		
		function labels
			a = get(gca);
			axis([min(range) max(range) a.YLim]);
			plot([0 0],a.YLim,'r:');
			plot(a.XLim,[0 0],'r:');
			hold off;
			xlabel('angular velocity (deg/s)');
			ylabel('moment (Nm)');
		end
		
	end

	%=============================================================================================================
	function [mom] = maxmoment(joint, angles, angvel, sign)
		% simulate maximum moment at one joint, as function of all joint angles and angular velocities
		% joint		for which joint we will calculate moment
		% angles	the six joint angles (deg)
		% angvel	the six angular velocities (deg/s)
		% sign		0: passive, 1: max positive moment, -1: max negative moment

		angles = angles*pi/180;		% convert to radians
		angvel = angvel*pi/180;
		
		% determine moment arms so we know which muscles to activate
		% here (2D model) we have constant moment arms.
		% we should in a later version ask the MEX function what the moment arms are at this posture
		momentarms_oneside = [	0.05	0		0		; ...
							-0.062	0		0		; ...
							-0.072	-0.034	0		; ...
							0.034	0.05	0		; ...
							0		0.042	0		; ...
							0		-0.02	-0.053	; ...
							0		0		-0.053	; ...
							0		0		0.037	];
		% make a matrix that contains moment arms for both sides (in case needed...)
		momentarms = [	momentarms_oneside 		zeros(nmus/2,njoints/2) ; ...
						zeros(nmus/2,njoints/2) 	momentarms_oneside];
		
		Act =  sign*momentarms(:,joint) > 0;	% vector that has a 1 for all muscles we want to activate
		
		% determine lengthening velocities of the muscle-tendon complexes, normalize to Lceopt
		Vmuscle = -(momentarms * angvel) ./ Lceopt;
		
		% determine the Lce's at which there is contraction equilibrium (dF/dt = 0, or Lcedot = Vmuscle)
		% use Newton's method
		xdot = [zeros(2*ndof,1) ; Vmuscle ; zeros(nmus,1)];	% we want these state derivatives
		u = zeros(nmus,1);				% no stim, we don't care about activation dynamics here
		tol = 1e-6;						% tolerance of iterative Newton solver
		indexLce = 2*ndof + (1:nmus);	% index to Lce variables within the state vector x
		x = [0 0 0 angles' 0 0 0 angvel' zeros(1,nmus) Act']';
		
		for imus=1:nmus;		% max number of Newton iterations
			Lce = 1.0;		% initial guess for this muscle's Lce			
			[Lce, Fval, Flag] = fzero(@contraction_equilibrium, Lce);				
		end
		
		if (flag < 0)
			fprintf('maxmoment: muscle contraction equilibrium not found within max number of iterations.\n');
			keyboard
		end
		
		% now determine the joint moments at this state of the system
		moments = gait2d('Jointmoments',x);
		mom = moments(joint);
		
		function [F] = contraction_equilibrium(Lce)
			x(2*ndof+imus) = Lce;
			f = gait2d('Dynamics',x,xdot,u);
			F = f(2*ndof+imus);
		end

	end

	%=====================================================
	% the following function solves the state derivatives from state and control,
	% so we can simulate using a Matlab ODE solver
	function [xdot] = odefun(t,x)
		u = 0.5+zeros(16,1);				% all muscles at 50% stim
		if exist('fsolve.m') > 1
			% use matlab's fsolve if optimization toolbox is installed
			options = optimset(optimset('fsolve'),'Jacobian','on','Display','off');
			xdot = fsolve(@equilib2,zeros(nstates,1),options);
		else
			% use a simple Newton-Raphson method
			% this actually turns out to be much faster than fsolve and equally accurate
			% this converges for the standard test (freefall with 50% stim) but this is not guaranteed for other simulations!
			xdot = zeros(size(x));
			resnorm = 1e10;
			while resnorm > 1e-6
				[f,J] = equilib2(xdot);
				xdot = xdot - J\f;		% do one full Newton step
				resnorm = norm(f);			% calculate norm of residuals
			end
		end
		
		function [f,J] = equilib2(xdot)
			[f,dfdx,J] = gait2d('Dynamics',x, xdot, u);	
		end
		
	end
	
end
%=====================================================
function matcompare(a,b);
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('Max. difference: %8.6f at %d %d (%8.6f vs. %8.6f)\n', ...
		maxerr, irow, icol, a(irow,icol),b(irow,icol));
end
%=========================================================
function accplot(x,y,xstring)
	for iseg = 1:7
		for ia = 1:6
			ivar = 6*(iseg-1)+ia;
			subplot(7,6,ivar);
			plot(x,y(:,ivar));
			if (iseg==1)
				title(['a' num2str(ia)]);
			end
			if (iseg==7)
				xlabel(xstring);
			end
			if (ia==1)
				ylabel(['segment ' num2str(iseg)]);
			end
		end
	end
	disp('Hit ENTER to continue...');pause;
end
