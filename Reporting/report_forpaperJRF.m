function report_forpaperJRF(resultfile)
	% make report for R3FAS 2D model results

    if nargin == 0
        resultfile1 = 'Winter\Winter_normal_result_able.mat';
        resultfile2 = 'Winter/Paper/Result_mom_1e-08_effort_10_track_1.mat';
        resultfile3 = 'Winter\Paper\Result_mom_0.3_effort_10_track_1.mat';
        resultfile4 = 'Winter\Paper\Result_mom_1_effort_1_track_0.1.mat';
    end
    
	global lincolors1
	
	% Constants
	musclenames = {'Iliopsoas' 'Vastus'};
	nmus = size(musclenames,2);

    resultfile = char(resultfile1, resultfile2, resultfile3, resultfile4);
    

    % create Figure window
	figure(	'NumberTitle', 		'off',	...
			'Name',				'Report of Multiple Files', ...
			'PaperOrientation',	'landscape', ...
			'PaperPosition',	[0 0 10 8.5], ...		% fill the page
			'OuterPosition',	[1 1 1024 1024]);
        hold on
	lincolors1 = ['b' 'k', 'r', 'g'];

 
    for ii = 1:size(resultfile,1)
        load(resultfile(ii,:));
        x = result.x;
        u = result.u;
        dur = result.dur;
        speed = result.speed;
        N = size(x,2);
        t = 100*(0:N)/N;
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
        simjrf = JRF;
        simjrf(4:6,:) = [simjrf(4:6,N/2+1:end) simjrf(4:6,1:N/2)];
        simjrf = [simjrf simjrf(:,1)]/BW;
        simgrf = GRF;
        simgrf(3:4,:) = [simgrf(3:4,N/2+1:end) simgrf(3:4,1:N/2)];
        simgrf = [simgrf simgrf(:,1)]/BW;
        
        % Plot Joint Reaction forces

        subplot(3,2,1)
        set(gca,'fontsize',12)
        hold on
        plot(t,simjrf(1,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 2.4])
        if ii == 1
            ylabel('Force [BW]', 'Fontsize', 12)
            title('Prosthesis Side', 'Fontsize', 14)
    		text(98,2,'Hip Contact Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
        
        subplot(3,2,2)
        set(gca,'fontsize',12)
        hold on
        plot(t,simjrf(4,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 2.4])
        if ii == 1
            ylabel('Force [BW]', 'Fontsize', 12)
            title('Intact Side', 'Fontsize', 14)
            text(45,2,'Hip Contact Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
        
        subplot(3,2,3)
        set(gca,'fontsize',12)
        hold on
        plot(t,simjrf(2,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 3.9])
        if ii == 1
            ylabel('Force [BW]', 'Fontsize', 12)
            text(98,3.5,'Knee Contact Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
        
        subplot(3,2,4)
        set(gca,'fontsize',12)
        hold on
        plot(t,simjrf(5,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 3.9])
        if ii == 1
            ylabel('Force [BW]', 'Fontsize', 12)
            text(50,3.5,'Knee Contact Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
         
        subplot(3,2,5)
        set(gca,'fontsize',12)
        hold on
        plot(t,simgrf(2,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 1.4])
        if ii == 1
            xlabel('% of Gait Cycle', 'Fontsize', 12);
            ylabel('Force [BW]', 'Fontsize',12);
            text(98,1.28,'Ground Reaction Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
        subplot(3,2,6)
        set(gca,'fontsize',12)
        hold on
        plot(t,simgrf(4,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-0.01 1.4])
%         ylim = 
        if ii == 1
            xlabel('% of Gait Cycle', 'Fontsize', 12);
            ylabel('Force [BW]', 'Fontsize',12);
            text(55,1.28,'Ground Reaction Force','HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12);
        end
        
    end
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
