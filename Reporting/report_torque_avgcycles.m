function report_torque_avgcycles(resultfile)
	% make report for R3FAS 2D model results

	global linecolors
	

	% Load the file with simulation result
	load(resultfile)
	x = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
	N = result.problem.N;
    Ncycles = result.problem.Ncycles;
    Ntot = N*Ncycles;
    symm = result.problem.symmetry;
    
    % rearrange with multiple gait cycles and symmetry
    if symm
        xnew = [];
        dist = 0;
        for j = 0:Ncycles-1
            x1 = x(:,N*j+(1:N));
            x1(1,:) = x1(1,:) + dist;
            x2 = x1(result.problem.vmx,:);
            x2(1,:) = x2(1,:) + speed*dur(j+1)/2;
            xnew = [xnew x1 x2];
            dist = dist + speed*dur(j+1)/2;
        end
        %also repeat input
        u1 = u(:,1:N);
        u2 = u1(result.problem.vmu,:);
        u = [u1 u2];
        x = xnew;
        Ntot = Ntot*2;
        N = N*2;
    end
    x = [x(1:18,:); zeros(32,Ntot)];
    xnew = reshape(x,50,60,Ncycles);
%     xmean = x;
    xmean = mean(xnew,3);
    xstd = std(xnew,[],3);
    K_mat = transpose(repmat(result.K,1,6));
    K_mat = [zeros(6,3) diag(K_mat(:,1)) zeros(6,3) diag(K_mat(:,2))];
    
    % Calculate torque with feedback
    utot = zeros(6,N,Ncycles);
    for j = 1:size(xnew,3)
        for i = 1:size(xnew,2)
            utot(:,i,j) = u(:,i) + K_mat*xnew(1:18,i,j)*10;
        end
    end
    
    utot_avg = mean(utot,3);
    
    % Vector of time instants for nodes
	T = (1:N)*sum(dur)/N; 

	% Initialize the model (we need to run the model to compute GRFs and muscle forces
	model = initmodel(result.model);

	% Computing GRFs, muscle forces
	GRF = zeros(4,N);
	ycontact = zeros(4,N);
	for i=1:N
		GRF(:,i) = model.gait2d('GRF',xmean(1:50,i));
		d = model.gait2d('Stick',xmean(1:50,i));
		ycontact(:,i) = d([5 6 9 10],2)';		% Y coordinates of right heel and toe, left heel and toe
	end

	% determine body weight
	BW = 9.81 * model.mass;

	% rearrange the matrices and unit conversion
	simang = 180/pi*[xmean([4:6],:) xmean([7:9],:)]';
    simangstd = 180/pi*[xstd([4:6],:) xstd([7:9],:)]';
	simangvel = [xmean([13:15],:) xmean([16:18],:)]';
	simgrf = [GRF([1 2],:) GRF([3 4],:)]'/BW;
	simmom = [utot_avg([1:3],:) utot_avg([4:6],:)]';
    simutot = [utot_avg([1:3],:) utot_avg([4:6],:)]';
    simycon = [ycontact(2,:) ycontact(4,:)]';
	
	% calculate joint power
    simpwr = simutot.*simangvel;
	 
	% create Figure window
% 	close all				% closes all existing Figure windows
	figure(	'NumberTitle', 		'off',	...
			'Name',				['Report of: ' path resultfile], ...
			'PaperOrientation',	'landscape', ...
			'PaperPosition',	[0 0 10 8.5], ...		% fill the page
			'OuterPosition',	[1 1 1024 1024]);
	linecolors = ['b' 'r'];

	% Optimization Information
	subplot(3,4,1)
	xtext1 = -0.4;
	title(strrep(resultfile,'_','\_'));
	text(xtext1,1,['# Nodes: ' num2str(N)],'FontSize',7)
	text(xtext1,0.85,['Speed: ' num2str(speed,'%8.3f') ' m/s'],'FontSize',7)
	text(xtext1,0.7,['Gait cycles: ', num2str(sum(dur),'%8.3f'), ' s'],'FontSize',7);
	text(xtext1,0.55,['Wtrack: ', num2str(model.Wtrack,'%8.3f')],'FontSize',7);
	text(xtext1,0.4,['Weffort: ', num2str(model.Weffort,'%8.3f')],'FontSize',7);
	text(xtext1,0.25,['Wreg: ', num2str(model.Wreg,'%8.3f')],'FontSize',7);
	text(xtext1,0.1,['Datafile: ', strrep(model.datafile,'_','\_')],'FontSize',7);
	
	xtext2 = 0.4;
	text(xtext2,1,['Objective: ' num2str(result.f,'%8.5f')],'FontSize',7)
	text(xtext2,0.85,['Norm(c): ' num2str(result.normc,'%8.5f')],'FontSize',7)
	
	axis off
	box on

	% Horizontal GRF
	subplot(3,4,5)
	plotvar(simgrf(:,1), 0,length(dur));
	ylabel('Force [BW]', 'Fontsize',8');
	title('Horizontal GRF');

	% Vertical Ground reaction force
	subplot(3,4,9)
	plotvar(simgrf(:,2), 0,length(dur));
	ylabel('Force [BW]', 'Fontsize',8');
	title('Vertical GRF');
	
	% Plot joint angles
	subplot(3,4,2)
	plotvar(simang(:,1), 0,length(dur),simang(:,1),simangstd(:,1));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Hip angle');
	subplot(3,4,3)
	plotvar(-simang(:,2), 0,length(dur),-simang(:,2),simangstd(:,2));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Knee angle');
	subplot(3,4,4)
	plotvar(simang(:,3), 0,length(dur),simang(:,3),simangstd(:,3));
	ylabel('Angle [deg]', 'Fontsize',8')
	title('Ankle angle');

	% Plot joint moments
	subplot(3,4,6)
	plotvar(-simmom(:,1), 0,length(dur));
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Hip moment');
	subplot(3,4,7)
	plotvar(simmom(:,2), 0,length(dur));
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Knee moment');
	subplot(3,4,8)
	plotvar(-simmom(:,3), 0,length(dur));
	ylabel('Moment [N m]', 'Fontsize',8')
	title('Ankle moment');
	
	% Plot joint powers
	subplot(3,4,10);
	plotvar(simpwr(:,1), 0,length(dur));
	ylabel('Power [W]', 'Fontsize',8')
	title('Hip power');
	subplot(3,4,11)
	plotvar(simpwr(:,2), 0,length(dur));
	ylabel('Power [W]', 'Fontsize',8')
	title('Knee power');
	subplot(3,4,12)
	plotvar(simpwr(:,3), 0,length(dur));
	ylabel('Power [W]', 'Fontsize',8')
	title('Ankle power');

	a = legend('Right','Left');
	p = get(a,'Position');
	p(1) = 0.72;
	p(2) = 0.13;
	set(a,'Position',p);
	
    figure
    plotvar(simycon*100,0,length(dur))
    xlabel('% of Gait Cycle')
    ylabel('Foot Clearance [cm]')
end
%============================================================================
function plotvar(sim,bsymm,Ncycles,av,sd);
	
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
