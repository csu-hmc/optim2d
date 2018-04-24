function anim(resultfile);
% anim.m: make a movie of the full gait cycle

	% Load the file with simulation result
	if (nargin == 0)
		[file,path] = uigetfile('*.mat');
		resultfile = [path file];
	end
	load(resultfile);
	x = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
	N = size(x,2)/length(dur);
    Ncycles = length(dur);
	gait2d = result.model.gait2d;		% function handle for musculoskeletal model
    if ~isfield(result, 'problem')
        result.problem.symmetry = 0;
    end
	symm = result.problem.symmetry;
    
    % number of repetitions, for smoother playback
	Nrep = ceil(5/Ncycles);	
    

	% Vector of time instants for nodes
	T = (1:N)*sum(dur)/N; 

	% Initialize the model (we need to run the model to compute GRFs and muscle forces
	result.model = initmodel(result.model);

	% if needed, create full gait cycle
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
        x = xnew;
        N = N*2;
    end
    
    if strcmp(result.model.type, 'torque')
        x = [x(1:18,:); zeros(32,N*Ncycles)];
    end
    
	actualfps = size(x,2)/sum(dur);
	
	% initialize movie file
	moviefile = strrep(resultfile,'.mat','.avi');
	avi = VideoWriter(moviefile);%, 'fps', fps, 'compression', 'Cinepak');

	% initialize figure window
	close all
	figure(1);
	clf;
	set(gcf,'Position',[5 5 550 550]);
	set(gcf, 'color', 'white');
	
	% create ground points (random dots)
	np = 5000*max(Nrep,Ncycles);
	xg = rand(np,1);
	yg = rand(np,1);
	xg = -1 + 3*Ncycles*[xg ; 2-xg];
	yg = -0.15*[yg ; yg];
	
	% make the movie
    open(avi)
	R = [1:6 4];			% right stick points
	L = [2 7:10 8];			% left stick points
	nframes = size(x,2);
	shift = speed*sum(dur)/nframes;
	for k = 1:Nrep
		for i=0:nframes-1
			plot(xg,yg,'.','Color',[0.7 0.7 0.7],'MarkerSize',4);
			hold on
			d = result.model.gait2d('Stick',x(1:50,i+1));
			plot(d(R,1),d(R,2),'b',d(L,1),d(L,2),'r','LineWidth',2);
			axis('equal');
			axis([ [-1 1]+shift*i -0.4 2]);
			axis('off');
			if (i==0)
				F = getframe(gca);
				frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
			else
				F = getframe(gca,frame);
			end
            writeVideo(avi,F);
			cla;
		end
	end
	close(avi);
	hold off;
	close all
end