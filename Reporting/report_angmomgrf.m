function report_angmomgrf(resultfile)
	% make report for R3FAS 2D model results

	global linecolors
	
	% Constants
	musclenames = {'Iliopsoas' 'Vastus'};
	nmus = size(musclenames,2);

	% Load the file with simulation result
	if (nargin == 0)
		[resultfile,path] = uigetfile('./*.mat');
		load([path resultfile]);
	else
		path = [];
		load(resultfile);
	end
	x = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
	N = size(x,2);
	symm = 0;							% no symmetry is assumed in the R3FAS project

	% Vector of time instants for nodes
	T = (1:N)*dur/N; 

	% Initialize the model (we need to run the model to compute GRFs and muscle forces
	model = initmodel(result.model);

	% Computing GRFs, muscle forces
	GRF = zeros(4,N);
	mfor = zeros(16,N);
	mom = zeros(6,N);
	tmp = zeros(16,N);
	ycontact = zeros(4,N);
	for i=1:N
		[GRF(:,i), dGRFdx, tmp(:,i)] = model.gait2d('GRF',x(1:50,i));
		mfor(:,i) = model.gait2d('Muscleforces',x(1:50,i));
		mom(:,i) = model.gait2d('Jointmoments',x(1:50,i));
		d = model.gait2d('Stick',x(1:50,i));
		ycontact(:,i) = d([5 6 9 10],2)';		% Y coordinates of right heel and toe, left heel and toe
    end
    [R_JRF,L_JRF] = getJointReactionForces(model, x);
    JRF = getJRF(R_JRF, L_JRF);
    
    % determine body weight
	BW = 9.81 * model.mass;

	% rearrange the matrices and unit conversion
	simang = 180/pi*[x([4:6],:) x([7:9],:)]';
	simgrf = [GRF([1 2],:) GRF([3 4],:)]'/BW;
	simmom = [mom([1:3],:) mom([4:6],:) ]';
	simMprosth = [x(54,:) zeros(1,N)]';
    simfor = [mfor([1:8],:) mfor([9:16],:) ]';
    simjrf = [JRF(1:3,:) JRF(4:6,:)]';
	
	% add the prosthesis moment to the knee moment from the musculoskeletal model
	simmom(:,2) = simmom(:,2) + simMprosth;

	% create Figure window
	figure(	'NumberTitle', 		'off',	...
			'Name',				['Report of: ' path resultfile], ...
			'PaperOrientation',	'landscape', ...
			'PaperSize',[11 8.5]);%, ...
% 			'OuterPosition',	[1 1 1024 1024]);
	linecolors = ['r' 'r--'];
    
    % Plot joint angles
    subplot(3,3,1)
    set(gca,'fontsize',12)
	plotvar(simang(:,1), symm);
	ylabel('Angle [deg]', 'Fontsize',12)
	title('Hip angle');
	subplot(3,3,2)
    set(gca,'fontsize',12)
	plotvar(-simang(:,2), symm);
	ylabel('Angle [deg]', 'Fontsize',12)
	title('Knee angle');
    
	% Plot joint moments
	subplot(3,3,4)
    set(gca,'fontsize',12)
	plotvar(-simmom(:,1), symm);
%     xlabel('% of Gait Cycle', 'Fontsize', 12);
	ylabel('Moment [N m]', 'Fontsize',12)
	title('Hip moment');
	subplot(3,3,5)
    set(gca,'fontsize',12)
	plotvar(simmom(:,2), symm);
%     xlabel('% of Gait Cycle', 'Fontsize', 12);
	ylabel('Moment [N m]', 'Fontsize',12)
	title('Knee moment');
    
    subplot(3,3,7)
    set(gca,'fontsize',12)
    plotvar(simjrf(:,1),symm);
    ylabel('Force [N]', 'Fontsize', 12)
    title('Hip Joint Reaction Force')
    
    subplot(3,3,8)
    set(gca,'fontsize',12)
    plotvar(simjrf(:,2),symm);
    ylabel('Force [N]', 'Fontsize', 12)
    title('Knee Joint Reaction Force')
    
    subplot(3,3,9)
    set(gca,'fontsize',12)
    plotvar(simjrf(:,3),symm);
    ylabel('Force [N]', 'Fontsize', 12)
    title('Ankle Joint Reaction Force')
    
% 
% 	% Plot muscle forces
% 	for i = 1:nmus
% 		subplot(4,2,i+4)
%         set(gca,'fontsize',12)
% 		plotvar(simfor(:,i), [], symm);
% 		ylabel('Force [N]', 'Fontsize',12);
%         xlabel('% of Gait Cycle', 'Fontsize', 12);
% 		title(musclenames(i));
%     end
end
%============================================================================
function plotvar(sim,bsymm,av,sd);
	
	global linecolors

	hold on	
	% do we have symmetry or do we plot left and right separately?
	if (bsymm)
		N = size(sim,1);
		t = 100*(0:N)/N;
		if (nargin>3)
			x = [t  t(end:-1:1)];
			Rav = [av(1:N); av(1)];
			Rsd = [sd(1:N); sd(1)];		
			y = [Rav-Rsd ; Rav(end:-1:1)+Rsd(end:-1:1)];	
			fill(x,y,[0.9 0.9 1]);
		end
		plot(t, [sim; sim(1)],linecolors(1), 'LineWidth', 2);
	else
		N = size(sim,1)/2;
		t = 100*(0:N)/N;
		Rsim = [sim(1:N); sim(1)];
		Lsim = [sim(N+1:2*N); sim(N+1)];
		if (nargin>3)
			x = [t  t(end:-1:1)];
			Rav = [av(1:N); av(1)];
			Lav = [av(N+1:2*N); av(N+1)];
			Rsd = [sd(1:N); sd(1)];
			Lsd = [sd(N+1:2*N); sd(N+1)];			
			Ry = [Rav-Rsd ; Rav(end:-1:1)+Rsd(end:-1:1)];	
			Ly = [Lav-Lsd ; Lav(end:-1:1)+Lsd(end:-1:1)];	
			fill(x,Ry,[0.9 0.9 1]);
			fill(x,Ly,[1 0.9 0.9]);
		end
		plot(t, Rsim,linecolors(1),'LineWidth', 2);
		plot(t, Lsim,linecolors(2),'LineStyle', '--','LineWidth', 2);
	end
	box on
	hold off
	if (nargin>3)
		ymin = min([av-sd;sim]);
		ymax = max([av+sd;sim]);
	else
		ymin = min(sim);
		ymax = max(sim);
	end
	margin = 0.05*(ymax-ymin)+ 1e-6;
	ymin = ymin - margin;
	ymax = ymax + margin;
	axis([0 100 ymin ymax]);
	box on
end
