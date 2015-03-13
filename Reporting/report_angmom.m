function report_angmom(resultfile)
	% make report for R3FAS 2D model results

	global linecolors

	% Load the file with simulation result
	if (nargin == 0)
		[resultfile,path] = uigetfile('./*.mat');
		load([path resultfile]);
	else
		path = [];
		load(resultfile);
	end
	x = result.x;
	dur = result.dur;
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
		[GRF(:,i), ~, tmp(:,i)] = model.gait2d('GRF',x(1:50,i));
		mfor(:,i) = model.gait2d('Muscleforces',x(1:50,i));
		mom(:,i) = model.gait2d('Jointmoments',x(1:50,i));
		d = model.gait2d('Stick',x(1:50,i));
		ycontact(:,i) = d([5 6 9 10],2)';		% Y coordinates of right heel and toe, left heel and toe
    end

	% rearrange the matrices and unit conversion
	simang = 180/pi*[x([4:6],:) x([7:9],:)]';
	simmom = [mom([1:3],:) mom([4:6],:) ]';
	simMprosth = [x(54,:) zeros(1,N)]';
	
	% add the prosthesis moment to the knee moment from the musculoskeletal model
	simmom(:,2) = simmom(:,2) + simMprosth;
	
	% Extract stored gait data that was tracked
	av = model.data.av;
	sd = model.data.sd;
	
	% Replicate the data as needed to match the number of simulated nodes (N)
	nrep = ceil(N/size(av,1));
	av = repmat(av,nrep,1);
	sd = repmat(sd,nrep,1);
	av = av(1:N,:);
	sd = sd(1:N,:);
	
	% rearrange data and unit conversion, put left after right side in same columns
	avang = 180/pi*[av(:,[1:3]) ; av(:,[6:8])];
	sdang = 180/pi*[sd(:,[1:3]) ; sd(:,[6:8])];

	% create Figure window
% 	close all				% closes all existing Figure windows
	figure(	'NumberTitle', 		'off',	...
			'Name',				['Report of: ' path resultfile])
	linecolors = ['r' 'r--'];

	% Plot joint angles
		% Plot joint angles
	subplot(2,3,1)
	plotvar(simang(:,1), symm);%, avang(:,1), sdang(:,1));
	ylabel('Angle [deg]', 'Fontsize',12')
	title('Hip angle', 'Fontsize',12');
	subplot(2,3,2)
	plotvar(-simang(:,2), symm);%, -avang(:,2), sdang(:,2));
	ylabel('Angle [deg]', 'Fontsize',12')
	title('Knee angle', 'Fontsize',12');
	subplot(2,3,3)
	plotvar(simang(:,3), symm);%, avang(:,3), sdang(:,3));
	ylabel('Angle [deg]', 'Fontsize',12')
	title('Ankle angle', 'Fontsize',12');

	% Plot joint moments
	subplot(2,3,4)
	plotvar(-simmom(:,1), symm);
	ylabel('Moment [N m]', 'Fontsize',12')
    xlabel('Percentage of Gait Cycle', 'Fontsize',12')
	title('Hip moment', 'Fontsize',12');
	subplot(2,3,5)
	plotvar(simmom(:,2), symm);
	ylabel('Moment [N m]', 'Fontsize',12')
	title('Knee moment', 'Fontsize',12');
    xlabel('Percentage of Gait Cycle', 'Fontsize',12')
	subplot(2,3,6)
	plotvar(-simmom(:,3), symm);
	ylabel('Moment [N m]', 'Fontsize',12')
    xlabel('Percentage of Gait Cycle', 'Fontsize',12')
	title('Ankle moment', 'Fontsize',12');

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
		plot(t, [sim; sim(1)],linecolors(1));
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
		plot(t, Rsim,linecolors(1), 'Linewidth', 1.5);
		plot(t, Lsim,linecolors(2), 'Linestyle', '--', 'Linewidth', 1.5);
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
