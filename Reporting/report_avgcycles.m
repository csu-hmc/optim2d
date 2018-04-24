function report_avgcycles(resultfile)
	% make report for R3FAS 2D model results

	global linecolors
	
    musclenames = {'Iliopsoas' 'Gluteals' 'Hamstrings' 'Rectus Femoris' 'Vastus' 'Gastrocnemius' 'Soleus' 'Tibialis Anterior'};
	
	% Load the file with simulation result
	load(resultfile)
	x = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
    if ~isfield(result,'problem')
        N = 60;
        Ncycles = 1;
        symm = 0;
        x = x(1:50,:);
    else
        N = result.problem.N;
        Ncycles = result.problem.Ncycles;
        symm = result.problem.symmetry;
    end
    Ntot = N*Ncycles;
    h = sum(dur)/Ntot;
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
    xnew = reshape(x,50,N,Ncycles);
%     xmean = x;

    % Initialize the model (we need to run the model to compute GRFs and muscle forces
	model = initmodel(result.model);
    GRF = zeros(4,N,Ncycles);
    mfor = zeros(16,N,Ncycles);
    mom = zeros(6,N,Ncycles);
    ycontact = zeros(4,N,Ncycles);
    % Calculate moments and ground reaction force for each cycle
    for j = 1:Ncycles
        for i=1:N
            GRF(:,i,j) = model.gait2d('GRF',xnew(1:50,i,j));
            mfor(:,i,j) = model.gait2d('Muscleforces',xnew(1:50,i,j));
            mom(:,i,j) = model.gait2d('Jointmoments',xnew(1:50,i,j));
            d = model.gait2d('Stick',xnew(1:50,i,j));
            ycontact(:,i,j) = d([5 6 9 10],2)';		% Y coordinates of right heel and toe, left heel and toe
        end
    end
        
    xmean = mean(xnew,3);     xstd = std(xnew,[],3);
    GRFmean = mean(GRF,3);     GRFstd = std(GRF,[],3);
    mformean = mean(mfor,3);     mforstd = std(mfor,[],3);
    mommean = mean(mom,3);     momstd = std(mom,[],3);
    ycontactmean = mean(ycontact,3);     ycontactstd = std(ycontact,[],3);
    
    %Rearrange to start at right heel strike
    indstart = find(GRFmean(2,:) > 50,1)-2;
    if indstart < 1
        indstart = 1;
    end
    xmean = xmean(:,[indstart:end 1:indstart-1]);             xstd = xstd(:,[indstart:end 1:indstart-1]);
    GRFmean = GRFmean(:,[indstart:end 1:indstart-1]);         GRFstd = GRFstd(:,[indstart:end 1:indstart-1]);
    mformean = mformean(:,[indstart:end 1:indstart-1]);       mforstd = mforstd(:,[indstart:end 1:indstart-1]);
    mommean = mommean(:,[indstart:end 1:indstart-1]);         momstd = momstd(:,[indstart:end 1:indstart-1]);
    ycontactmean = ycontactmean(:,[indstart:end 1:indstart-1]); ycontactstd = ycontactstd(:,[indstart:end 1:indstart-1]);
    u = u(:,[indstart:end 1:indstart-1]);
    
    %compute stimulation with feedback
    D = [1	0	0	0	0	0; -1	0	0	0	0	0; -1	-1	0	0	0	0;
            1	1	0	0	0	0; 0	1	0	0	0	0; 0	-1	-1	0	0	0;
            0	0	-1	0	0	0; 0	0	1	0	0	0; 0	0	0	1	0	0;
            0	0	0	-1	0	0; 0	0	0	-1	-1	0; 0	0	0	1	1	0;
            0	0	0	0	1	0; 0	0	0	0	-1	-1;0	0	0	0	0	-1;
            0	0	0	0	0	1];
        
    if ~isfield(result,'problem')
        K_mat = zeros(16,18);
        utot_avg = u;
        result.problem.nmus = 16;
    else
        K_mat = [zeros(result.problem.nmus,3) D*result.K(1)*10 zeros(result.problem.nmus,3) D*result.K(2)*10];
        utot = zeros(16,N,Ncycles);
        for j = 1:size(xnew,3)
            for i = 1:size(xnew,2)
                utot(:,i,j) = u(:,i) + K_mat*xnew(1:18,i,j);
            end
        end 
        utot_avg = mean(utot,3);
    end
    
    
    %Calculate CI index between Vastus and Hamstrings
    uVas = utot_avg(5,:);
    uHam = utot_avg(3,:);
    
    istart = 1/5*N;
    iend = 2/5*N;
    if uVas(istart) > uHam(istart)
        iHam = 1;
        iVas = [];
    else
        iVas = 1;
        iHam = [];
    end
    for i = istart+1:iend%length(uVas)
        if and(uVas(i)>0,uHam(i)>0)
            if uVas(i) > uHam(i)
                iHam = [iHam i];
            else
                iVas = [iVas i];
            end
        end
    end
        
    Iant = sum(uVas(iVas))*h +sum(uHam(iHam))*h;
    Itot = sum(uVas)*h+sum(uHam)*h;
    
    CI = 2*Iant/Itot
    
    %Isakov way
    istance = N*0.6;
    
    ST1_ham = sum(uHam(1:istance/2));
    ST2_ham = sum(uHam(istance/2+1:istance));
    ST1_vas = sum(uVas(1:istance/2));
    ST2_vas = sum(uVas(istance/2+1:istance));
    
    CI1 = ST1_ham/ST1_vas;
    CI2 = ST2_ham/ST2_vas;
    % Vector of time instants for nodes
	T = (1:N)*sum(dur)/N; 
    
	% determine body weight
	BW = 9.81 * model.mass;

	% rearrange the matrices and unit conversion
	simang = 180/pi*[xmean([4:6],:) xmean([7:9],:)]';
    simangstd = 180/pi*[xstd([4:6],:) xstd([7:9],:)]';
	simangvel = [xmean([13:15],:) xmean([16:18],:)]';
	simgrf = [GRFmean([1 2],:) GRFmean([3 4],:)]'/BW;
    simstim = [utot_avg([1:8],:) utot_avg([9:16],:) ]';	
    simfor = [mformean([1:8],:) mformean([9:16],:) ]';
    simmom = [mommean([1:3],:) mommean([4:6],:) ]';
    simycon = [ycontactmean(2,:) ycontactmean(4,:)]';
	
	% calculate joint power
    simpwr = simmom.*simangvel;
	 
	% create Figure window
% 	close all				% closes all existing Figure windows
	figure(	'NumberTitle', 		'off',	...
			'Name',				['Report of: ' path resultfile], ...
			'PaperOrientation',	'landscape', ...
			'PaperPosition',	[0 0 10 8.5], ...		% fill the page
			'OuterPosition',	[1 1 1024 1024]);
	linecolors = ['b' 'r'];

%	 Optimization Information
       

        % Horizontal GRF
        subplot(3,5,2)
        plotvar(simgrf(:,1), []);%, avgrf(:,1), sdgrf(:,1));
        ylabel('Force [BW]', 'Fontsize',8');
        title('Horizontal GRF');

        % Vertical Ground reaction force
        subplot(3,5,3)
        plotvar(simgrf(:,2), []);%, avgrf(:,2), sdgrf(:,2));
        ylabel('Force [BW]', 'Fontsize',8');
        title('Vertical GRF');

        % Plot joint angles
        subplot(3,5,6)
        plotvar(simang(:,1), [], simang(:,1), simangstd(:,1));
        ylabel('Angle [deg]', 'Fontsize',8')
        title('Hip angle');
        subplot(3,5,7)
        plotvar(-simang(:,2), [], -simang(:,2), simangstd(:,2));
        ylabel('Angle [deg]', 'Fontsize',8')
        title('Knee angle');
        subplot(3,5,8)
        plotvar(simang(:,3), [], simang(:,3), simangstd(:,3));
        ylabel('Angle [deg]', 'Fontsize',8')
        title('Ankle angle');

        % Plot joint moments
        subplot(3,5,11)
        plotvar(-simmom(:,1), []);
        ylabel('Moment [N m]', 'Fontsize',8')
        title('Hip moment');
        subplot(3,5,12)
        plotvar(simmom(:,2), []);
        ylabel('Moment [N m]', 'Fontsize',8')
        title('Knee moment');
        subplot(3,5,13)
        plotvar(-simmom(:,3), []);
        ylabel('Moment [N m]', 'Fontsize',8')
        title('Ankle moment');

        a = legend('Right','Left');
        p = get(a,'OuterPosition');
        p(1) = 0.72;
        p(2) = 0.13;
        set(a,'OuterPosition',p);

        % Plot muscle forces
        for i = 1:result.problem.nmus/2
            
            subplot(8,5,5*(i-1)+4) % 4 9 14 19 
            plotvar(simfor(:,i), []);
            ylabel('Force [N]', 'Fontsize',8');
            title(musclenames(i));
        end
        
        % Plot muscle stimulations
        for i = 1:result.problem.nmus/2
            subplot(8,5,5*(i-1)+5) % 4 9 14 19 
            plotvar(simstim(:,i), []);
            ylabel('Stimulation [-]', 'Fontsize',8');
            title(musclenames(i));
        end
	
	
    figure
    plotvar(simycon*100,0,length(dur))
    xlabel('% of Gait Cycle')
    ylabel('Foot Clearance [cm]')
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
			fill(x,Ry,[0.9 0.9 1], 'linestyle','none');
			fill(x,Ly,[1 0.9 0.9], 'linestyle','none');
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
