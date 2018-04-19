function report_forpaper(resultfile1, resultfile2, resultfile3, resultfile4)
	% make report for R3FAS 2D model results

    if nargin == 0
        resultfile1 = 'Winter\Winter_normal_result_able.mat';
        resultfile2 = 'Winter/Paper/Result_mom_1e-08_effort_10_track_1.mat';
        resultfile3 = 'Winter\Paper\Result_mom_0.3_effort_10_track_1.mat';
        resultfile4 = 'Winter\Paper\Result_mom_1_effort_1_track_0.1.mat';
    end
    
	global linecolors
	
	% Constants
	musclenames = {'Iliopsoas' 'Hamstrings'};
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

        % determine body weight
        BW = 9.81 * model.mass;

        %rearrange
        x = [x x(:,1)];
        x(7:9,:) = [x(7:9,N/2+1:end) x(7:9,1:N/2)];
        mom = [mom mom(:,1)];
        mom(4:6,:) = [mom(4:6,N/2+1:end) mom(4:6,1:N/2)];
        mfor = [mfor([1,9],:) mfor([1,9],1); mfor([3,11],:) mfor([3,11],1)];
        mfor([2,4],:) = [mfor([2,4],N/2+1:end) mfor([2,4],1:N/2)];
        mfor = mfor';

        % Plot joint angles
        subplot(4,2,1)
        hold on
        set(gca,'fontsize',12)
        plot(t,180/pi*x(4,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-20 50])
        if ii == 1
            text('FontSize',12,'String','Hip','Position',[46.8249258160237 36.3636363636364 0]);
            ylabel('Angle [deg]', 'Fontsize',12)
            title('Prosthesis Side', 'FontSize', 14);
        end
        
        subplot(4,2,2)
        hold on
        set(gca,'fontsize',12)
        plot(t,180/pi*x(7,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-20 50])
        if ii == 1
            text('FontSize',12,'String','Hip','Position',[46.2314540059347 35.6578947368421 0]);
            ylabel('Angle [deg]', 'Fontsize',12)
            title('Intact Side', 'FontSize', 14);
        end
        
        subplot(4,2,3)
        hold on
        set(gca,'fontsize',12)
        plot(t,-180/pi*x(5,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-5 70])
        if ii == 1
            ylabel('Angle [deg]', 'Fontsize',12)
            text('FontSize',12,'String','Knee','Position',[44.1543026706231 54.5394736842105 0]);
        end

        subplot(4,2,4)
        hold on
        set(gca,'fontsize',12)
        plot(t,-180/pi*x(8,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-5 70])
        if ii == 1
            ylabel('Angle [deg]', 'Fontsize',12)
            text('FontSize',12,'String','Knee','Position',[43.5608308605341 53.5526315789474 0]);
        end
        
        % Plot joint moments
        subplot(4,2,5,'YTick',[-20 0 20 50 80])
        hold on
        set(gca,'fontsize',12)
        plot(t,-mom(1,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-30 95])
        if ii == 1
            ylabel('Moment [N m]', 'Fontsize',12)
            text('FontSize',12,'String','Hip','Position',[46.2314540059347 73.7662337662338 0]);
        end
        subplot(4,2,6,'YTick',[-20 0 20 50 80])
        hold on
        set(gca,'fontsize',12)
        plot(t,-mom(4,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-30 95])
        if ii == 1
            ylabel('Moment [N m]', 'Fontsize',12)
            text('FontSize',12,'String','Hip','Position',[46.2314540059347 74.5454545454545 0]);
        end

        subplot(4,2,7)
        hold on
        set(gca,'fontsize',12)
        plot(t,mom(2,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-30 50])
        if ii == 1
            ylabel('Moment [N m]', 'Fontsize',12)
            xlabel('% of Gait Cycle', 'Fontsize', 12);
            text('FontSize',12,'String','Knee','Position',[44.1543026706231 38.025974025974 0]);
        end
        subplot(4,2,8)
        hold on
        set(gca,'fontsize',12)
        plot(t,mom(5,:),lincolors1(ii), 'LineWidth', 1.5)
        ylim([-30 50])
        if ii == 1
            ylabel('Moment [N m]', 'Fontsize',12)
            xlabel('% of Gait Cycle', 'Fontsize', 12);
            text('FontSize',12,'String','Knee','Position',[44.4510385756676 36.9480519480519 0]);
        end

%         % Plot muscle forces
%         for i = 1:4
%             subplot(4,2,i+8)
%             hold on
%             set(gca,'fontsize',12)
%             plot(t,mfor(:,i),lincolors1(ii), 'LineWidth', 1.5)
%             if ii == 1
%                 ylabel('Force [N]', 'Fontsize',12);
%                 xlabel('% of Gait Cycle', 'Fontsize', 12);
%                 if i == 1
%                     text('FontSize',12,'String','Iliopsoas','Position',[50.3857566765578 954.805194805194 0]);
%                 elseif i ==2
% %                     text('FontSize',12,'String','Iliopsoas','Position',[40.8902077151335 938.961038961039 0]);
%                 elseif i ==3
% %                     text('FontSize',12,'String','Hamstrings','Position',[39.406528189911 601.56015037594 0]);
%                 else
%                     text('FontSize',12,'String','Hamstrings','Position',[39.1097922848665 597.5 0]);
%                 end
%             end
%             if and(i==4,ii == 4)
%                 legend('ABLE','AMP1','AMP2','AMP3','Position',[0.864252646611359 0.259368839642474 0.0942460305929657 0.123274158385143]);
%             end
%         end
    end
    
saveas(gcf,'C:\Users\2625262\Documents\Conferences Papers\2015\Paper1\Pictures\AngMomupdated','png')
end

