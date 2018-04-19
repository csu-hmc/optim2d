function report_torque_footclearance_Dayton(movement, stdev, Ncycles, repeats)
	
    global linecolors
    
    % Create filenames
    for j = 1:repeats
        for k = 1:length(Ncycles)
            result(j,k,1).result = load([movement '_result_torque_1cycle_effort_test.mat']);
            result(j,k,1).data = getSimData(result(j,k,1).result.result);
            dur1(:,1) = repmat(result(j,k,1).result.result.dur,Ncycles,1);
            for i = 2:length(stdev)+1
                filename = [movement '_torque_truncated_' num2str(Ncycles(k)) 'cycles_effort_std' num2str(stdev(i-1)) '_time' num2str(j) '_test1.mat'];
                result(j,k,i).result = load(filename);
                result(j,k,i).data = getSimData(result(j,k,i).result.result);
                dur1(:,i) = result(j,k,i).result.result.dur;
                K(:,i) = result(j,k,i).result.result.K;
                footclearance(:,:,i) = result(j,k,i).data.simycon;
            end

        end
    end

    %Same in all
	dur = result(1,1,1).result.result.dur;
    speed = result(1,1,1).result.result.speed;
	N = result(1,1,1).result.result.problem.N;
    
    % Maybe add here to adjust peaks

    % Vector of time instants for nodes
	T = (1:N)*sum(dur)/N; 
    
	% create Figure window
	figure(	'NumberTitle', 		'off',	...
			'Name',				'Foot Clearance Report', ...
			'PaperOrientation',	'landscape');%, ...
% 			'PaperPosition',	[0 0 10 8.5], ...		% fill the page
% 			'OuterPosition',	[1 1 1024 1024]);
	linecolors = [1 0 0;1 0 1;0 0 1;0 0 0;0.494117647409439 0.184313729405403 0.556862771511078;...
        0 1 1;1 1 0;0.600000023841858 0.200000002980232 0; 1 0.600000023841858 0.7843137383461;...
        0.164705887436867 0.384313732385635 0.274509817361832;0 1 0];

    subplot(2,2,1)
    set(gca,'fontsize',12)
    hold on
    for j = 1:repeats
        for k = 1:length(Ncycles)
            for i = 1:(length(stdev)+1)
                plotvar(j*k*i,result(j,k,i).data.simycon(:,1)*100,0,Ncycles(k))
            end
        end
    end
    hold on
    plot([87 100], [0 0], 'k', 'Linewidth', 2)
    plot([87 100], [2.5 2.5], 'k', 'Linewidth', 2)
    plot([87 87], [0 2.5], 'k', 'Linewidth', 2)
    plot([100 100], [0 2.5], 'k', 'Linewidth', 2)
    xlabel('% of Gait Cycle')
    ylabel('Heel Clearance [cm]')
    subplot(2,2,3)
    set(gca,'fontsize',12)
    for j = 1:repeats
        for k = 1:length(Ncycles)
            for i = 1:(length(stdev)+1)
                plotvar(j*k*i,result(j,k,i).data.simycon(:,2)*100,0,Ncycles(k))
            end
        end
    end
    hold on
    plot([75 90], [-0.2 -0.2], 'k', 'Linewidth', 2)
    plot([75 90], [1 1], 'k', 'Linewidth', 2)
    plot([75 75], [-0.2 1], 'k', 'Linewidth', 2)
    plot([90 90], [-0.2 1], 'k', 'Linewidth', 2)
    xlabel('% of Gait Cycle')
    ylabel('Toe Clearance [cm]')
    legendCell = cellstr(num2str([0 stdev]', 'stdev = %-d Nm'));
    legend(legendCell)
    
    subplot(2,2,2)
    set(gca,'fontsize',12)
    hold on
    for j = 1:repeats
        for k = 1:length(Ncycles)
            for i = 1:(length(stdev)+1)
                plotvar(j*k*i,result(j,k,i).data.simycon(:,1)*100,0,Ncycles(k))
            end
        end
    end
    xlabel('% of Gait Cycle')
    ylabel('Heel Clearance [cm]')
    xlim([87 100])
    ylim([0 2.5])
    
    subplot(2,2,4)
    set(gca,'fontsize',12)
    for j = 1:repeats
        for k = 1:length(Ncycles)
            for i = 1:(length(stdev)+1)
                plotvar(j*k*i,result(j,k,i).data.simycon(:,2)*100,0,Ncycles(k))
            end
        end
    end
    xlabel('% of Gait Cycle')
    ylabel('Toe Clearance [cm]')
    xlim([75 90])
    ylim([-0.2 1])
end
%============================================================================
function data = getSimData(result);
    x = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
	N = result.problem.N;
    Ncycles = result.problem.Ncycles;
    Ntot = N*Ncycles;
    symm = result.problem.symmetry;
    ndof = result.problem.ndof;
    
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
    x = reshape(x,size(x,1),N,Ncycles);
    ufull = zeros([size(u) Ncycles]);
    for i = 1:Ncycles
        ufull(:,:,i) = u+x(4:ndof,:,i)*result.K(1)*10+x(ndof+(4:ndof),:,i)*result.K(2)*10;
    end
    
    xmean = mean(x,3);
    xstd = std(x,[],3);
    umean = mean(ufull,3);
    ustd = std(ufull,[],3);

	% Initialize the model (we need to run the model to compute GRFs and muscle forces
	model = initmodel(result.model);

	% Computing GRFs, muscle forces
	GRF = zeros(4,N);
	tmp = zeros(16,N);
	ycontact = zeros(4,N,Ncycles);
    for j = 1:Ncycles
        for i=1:N
            [GRF(:,i), ~, tmp(:,i)] = model.gait2d('GRF',x(1:50,i,j));
            d = model.gait2d('Stick',x(1:50,i,j));
            ycontact(:,i,j) = d([5 6 9 10],2)';		% Y coordinates of right heel and toe, left heel and toe
        end
    end
    yconmean = mean(ycontact,3);
    yconstd = std(ycontact,[],3);

	% determine body weight
	BW = 9.81 * model.mass;

	% rearrange the matrices and unit conversion
	data.simang = 180/pi*[xmean([4:6],:) xmean([7:9],:)]';
	data.simangvel = [xmean([13:15],:) xmean([16:18],:)]';
	data.simgrf = [GRF([1 2],:) GRF([3 4],:)]'/BW;
	data.simmom = [umean([1:3],:) umean([4:6],:)]';
    data.simycon = [ yconmean([3 4],:) yconmean([1 2],:)]';
    data.stdycon = [ yconstd([3 4],:) yconstd([1 2],:)]';
    
    % standard deviations
    data.stdang = 180/pi*[xstd([4:6],:) xstd([7:9],:)]';
	data.stdangvel = [xstd([13:15],:) xstd([16:18],:)]';
    data.stdmom = [ustd([1:3],:) ustd([4:6],:)]';
    
	% calculate joint power
	data.simpwr = data.simmom.* data.simangvel;
end
%============================================================================

function plotvar(ii,sim,bsymm,Ncycles,std);
	
	global linecolors

	hold on	
	% do we have symmetry or do we plot left and right separately?
	if (bsymm)
		N = size(sim,1);
		t = 100*(0:N)/N;
		if (nargin>4)
			x = [t  t(end:-1:1)];
			Rav = [sim(1:N); sim(1)];
			Rsd = [std(1:N); std(1)];		
			y = [Rav-Rsd ; Rav(end:-1:1)+Rsd(end:-1:1)];	
			fill(x,y,[0.9 0.9 1]);
		end
		plot(t, [sim; sim(1)],'Color', linecolors(:,ii),'linewidth',2);
	else
		N = size(sim,1)/2;
		t = 100*(0:N)/N;
		Rsim = [sim(1:N); sim(1)];
		Lsim = [sim(N+1:2*N); sim(N+1)];
		if (nargin>4)
			x = [t  t(end:-1:1)];
			Rav = [sim(1:N); sim(1)];
			Lav = [sim(N+1:2*N); sim(N+1)];
			Rsd = [std(1:N); std(1)];
			Lsd = [std(N+1:2*N); std(N+1)];			
			Ry = [Rav-Rsd ; Rav(end:-1:1)+Rsd(end:-1:1)];	
			Ly = [Lav-Lsd ; Lav(end:-1:1)+Lsd(end:-1:1)];	
			fill(x,Ry,[0.9 0.9 1]);
			fill(x,Ly,[0.9 0.9 1]);
		end
% 		plot(t, Rsim,'Color', linecolors(:,i));
		plot(t, Lsim,'Color', linecolors(ii,:),'linewidth',1.5);
	end
	box on
	hold off
	if (nargin>4)
		ymin = min([sim-std;sim]);
		ymax = max([sim+std;sim]);
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

%     % Adjust using max knee angle
%     [~,idx] = min(data.simang(:,2));
%     [~,idx1] = min(data1.simang(:,2));
%     if idx1 > idx
%         idxs = idx1-idx;
%         data1.simang = data1.simang([idxs:end 1:idxs-1] ,:);
%         data1.simangvel = data1.simangvel([idxs:end 1:idxs-1] ,:);
%         data1.simgrf = data1.simgrf([idxs:end 1:idxs-1] ,:);
%         data1.simmom = data1.simmom([idxs:end 1:idxs-1] ,:);
%         data1.simycon = data1.simycon([idxs:end 1:idxs-1] ,:);
%         data1.stdmang = data1.stdang([idxs:end 1:idxs-1] ,:);
%         data1.stdangvel = data1.stdangvel([idxs:end 1:idxs-1] ,:);        
%         data1.stdmom = data1.stdmom([idxs:end 1:idxs-1] ,:);
%         data1.stdycon = data1.stdycon([idxs:end 1:idxs-1] ,:);
%         data1.simpwr = data1.simpwr([idxs:end 1:idxs-1] ,:);
%     elseif idx1 < idx
%         idxs = idx-idx1;
%         data1.simang = data1.simang([end-idxs:end 1:end-idxs-1] ,:);
%         data1.simangvel = data1.simangvel([end-idxs:end 1:end-idxs-1] ,:);
%         data1.simgrf = data1.simgrf([end-idxs:end 1:end-idxs-1] ,:);
%         data1.simmom = data1.simmom([end-idxs:end 1:end-idxs-1] ,:);
%         data1.simycon = data1.simycon([end-idxs:end 1:end-idxs-1] ,:);
%         data1.stdmang = data1.stdang([end-idxs:end 1:end-idxs-1] ,:);
%         data1.stdangvel = data1.stdangvel([end-idxs:end 1:end-idxs-1] ,:);
%         data1.stdmom = data1.stdmom([end-idxs:end 1:end-idxs-1] ,:);
%         data1.stdycon = data1.stdycon([end-idxs:end 1:end-idxs-1] ,:);
%         data1.simpwr = data1.simpwr([end-idxs:end 1:end-idxs-1] ,:);
%     end
