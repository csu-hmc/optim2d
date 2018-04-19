function report_musclesforpaper(resultfile1, resultfile2, resultfile3, resultfile4)
	% make report for R3FAS 2D model results

    if nargin == 0
        resultfile1 = 'Winter\Winter_normal_result_able.mat';
        resultfile2 = 'Winter/Paper/Result_mom_1e-08_effort_10_track_1.mat';
        resultfile3 = 'Winter\Paper\Result_mom_0.3_effort_10_track_1.mat';
        resultfile4 = 'Winter\Paper\Result_mom_1_effort_1_track_0.1.mat';
    end
    
	global linecolors
	
	% Constants
	musclenames = {'Iliopsoas' 'Gluteals' 'Hamstrings' 'Rectus Femoris' 'Vastus' 'Gastrocnemius' 'Soleus' 'Tibialis Anterior'};
	nmus = size(musclenames,2);

    resultfile = char(resultfile1, resultfile2, resultfile3, resultfile4);
    
    global lincolors1
	

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

        mfor1 = [mfor mfor(:,1)];
        mfor1(nmus+1:end,:) = [mfor1(nmus+1:end,N/2+1:end) mfor1(nmus+1:end,1:N/2)];
        mfor = mfor1([1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16],:)';

        % Plot muscle forces
        for i = 1:nmus*2  
            subplot(8,2,i)
            hold on
            set(gca,'fontsize',12)
            plot(t,mfor(:,i),lincolors1(ii), 'LineWidth', 1.5)
            if and(i==1,ii == 1)
                title('Rectus Femoris', 'FontSize', 14);
            end
            if and(i==2,ii == 1)
                title('Hamstrings', 'FontSize', 14);
            end
            if or(and(i==15,ii == 1),and(i==16,ii == 1))
                xlabel('% of Gait Cycle', 'Fontsize', 12);
            end
            if ii == 1
                ylabel('Force [N]', 'Fontsize',12);
                text('FontSize',12,'String',musclenames(ceil(i/2)),'Position',[50.3857566765578 954.805194805194 0]);
            end
            if and(i==4,ii == 4)
                legend('ABLE','TTA1','TTA2','TTA3','Orientation', 'horizontal','Position',[0.864252646611359 0.259368839642474 0.0942460305929657 0.123274158385143]);
            end
        end
    end
end

