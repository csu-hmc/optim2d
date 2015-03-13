function report(resultfile)
	% make report for R3FAS 2D model results

	global linecolors
	
	% Constants
	musclenames = {'Iliopsoas' 'Gluteals' 'Hamstrings' 'Rectus Femoris' 'Vastus' 'Gastrocnemius' 'Soleus' 'Tibialis Anterior'};
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

	% determine body weight
	BW = 9.81 * model.mass;

	% rearrange the matrices and unit conversion
	simang = 180/pi*[x([4:6],:) x([7:9],:)]';
	simangvel = [x([13:15],:) x([16:18],:)]';
	simgrf = [GRF([1 2],:) GRF([3 4],:)]'/BW;
	simact = [u([1:8],:) u([9:16],:) ]';	
	simfor = [mfor([1:8],:) mfor([9:16],:) ]';
	simmom = [mom([1:3],:) mom([4:6],:) ]';
	simaccum = [x(51,:) NaN(1,N)]';
	simvalveu = [u(17,:) u(18,:)]';
	simvalvev = [x(52,:) x(53,:)]';
	simMprosth = [x(54,:) zeros(1,N)]';
	
	% add the prosthesis moment to the knee moment from the musculoskeletal model
	simmom(:,2) = simmom(:,2) + simMprosth;

	% calculate joint power
	simpwr = simmom .* simangvel;
	
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
	avgrf = [av(:,[4 5]) ; av(:,[9 10])];
	sdgrf = [sd(:,[4 5]) ; sd(:,[9 10])];

	% create Figure window
% 	close all				% closes all existing Figure windows
	figure(	'NumberTitle', 		'off',	...
			'Name',				['Report of: ' path resultfile], ...
			'PaperOrientation',	'landscape', ...
			'PaperPosition',	[0 0 10 8.5], ...		% fill the page
			'OuterPosition',	[1 1 1024 1024]);
	linecolors = ['b' 'r'];

	% Optimization Information
	subplot(6,4,1)
	xtext1 = -0.4;
	title(strrep(resultfile,'_','\_'));
	text(xtext1,1,['# Nodes: ' num2str(N)],'FontSize',7)
	text(xtext1,0.85,['Speed: ' num2str(speed,'%8.3f') ' m/s'],'FontSize',7)
	text(xtext1,0.7,['Gait cycle: ', num2str(dur,'%8.3f'), ' s'],'FontSize',7);
	text(xtext1,0.55,['Wtrack: ', num2str(model.Wtrack,'%8.3f')],'FontSize',7);
	text(xtext1,0.4,['Weffort: ', num2str(model.Weffort,'%8.3f')],'FontSize',7);
	text(xtext1,0.25,['Wvalve: ', num2str(model.Wvalve,'%8.3f')],'FontSize',7);
	text(xtext1,0.1,['Datafile: ', strrep(model.datafile,'_','\_')],'FontSize',7);
	
	xtext2 = 0.4;
	text(xtext2,1,['Objective: ' num2str(result.f,'%8.5f')],'FontSize',7)
	text(xtext2,0.85,['Norm(c): ' num2str(result.normc,'%8.5f')],'FontSize',7)
	
	axis off
	box on

	% Horizontal GRF
	subplot(6,4,5)
	plotvar(simgrf(:,1), symm, avgrf(:,1), sdgrf(:,1));
	ylabel('Force [BW]', 'Fontsize',8');
	title('Horizontal GRF');

	% Vertical Ground reaction force
	subplot(6,4,9)
	plotvar(simgrf(:,2), symm, avgrf(:,2), sdgrf(:,2));
	ylabel('Force [BW]', 'Fontsize',8');
	title('Vertical GRF');
	
	% Plot joint angles
	subplot(6,4,2)
	plotvar(simang(:,1), symm, avang(:,1), sdang(:,1));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Hip angle');
	subplot(6,4,3)
	plotvar(-simang(:,2), symm, -avang(:,2), sdang(:,2));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Knee angle');
	subplot(6,4,4)
	plotvar(simang(:,3), symm, avang(:,3), sdang(:,3));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Ankle angle');

	% Plot joint moments
	subplot(6,4,6)
	plotvar(-simmom(:,1), symm);
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Hip moment');
	subplot(6,4,7)
	plotvar(simmom(:,2), symm);
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Knee moment');
	subplot(6,4,8)
	plotvar(-simmom(:,3), symm);
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Ankle moment');
	
	% Plot joint powers
	subplot(6,4,10);
	plotvar(simpwr(:,1), symm);
	ylabel('Power [W]', 'Fontsize',8')
	title('Hip power');
	subplot(6,4,11)
	plotvar(simpwr(:,2), symm);
	ylabel('Power [W]', 'Fontsize',8')
	title('Knee power');
	subplot(6,4,12)
	plotvar(simpwr(:,3), symm);
	ylabel('Power [W]', 'Fontsize',8')
	title('Ankle power');

	a = legend('Right','Left');
	p = get(a,'Position');
	p(1) = 0.72;
	p(2) = 0.13;
	set(a,'Position',p);
	
	% Plot muscle forces
	for i = 1:nmus
		col = mod(i-1,4);
		row = (i-1-col)/4;
		k = 13+col+4*row;
		subplot(6,4,k)
		plotvar(simfor(:,i), [], symm);
		ylabel('Force [N]', 'Fontsize',8');
		title(musclenames(i));
	end
	
	% Plot prosthesis variables
	subplot(6,4,21);
	plotvar(simaccum, symm);
	title('accumulator');
	ylabel('Volume [cm^3]', 'Fontsize',8');
	xlabel('Time [% of movement]', 'Fontsize',8');
	
	linecolors = ['m' 'k'];
	subplot(6,4,22);
	plotvar(simvalveu, symm);
	title('valve controls');
	xlabel('Time [% of movement]', 'Fontsize',8');
	
	subplot(6,4,23);
	plotvar(simvalvev, symm);
	title('flow rates');
	ylabel('Flow [cm^3 s^{-1}]', 'Fontsize',8');
	xlabel('Time [% of movement]', 'Fontsize',8');
	
	a = legend('Valve 1', 'Valve 2');
	p = get(a,'Position');
	p(1) = 0.82;
	p(2) = 0.13;
	set(a,'Position',p);

	% compute net joint power (for looking at net power sources and sinks)
	% Rmeanpwr = mean(simpwr(1:N,:))
	% Lmeanpwr = mean(simpwr(N+(1:N),:))
	% Totalpwr = sum([Rmeanpwr Lmeanpwr])
	
	% process the contact speed,force,power variables (for looking at contact power loss)
	% v = tmp(1:8,:);
	% F = tmp(9:16,:);
	% P = v.*F;
	
	% uncomment the following if you want to look at anything interactively
	% keyboard

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
		plot(t, Rsim,linecolors(1));
		plot(t, Lsim,linecolors(2));
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
