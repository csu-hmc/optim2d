function anim(resultfile);
% anim.m: make a movie of the full gait cycle

	% number of repetitions, for smoother playback
	Nrep = 5;		

	% Load the file with simulation result
	if (nargin == 0)
		[file,path] = uigetfile('*.mat');
		resultfile = [path file];
	end
	load(resultfile);
	x1 = result.x;
	u = result.u;
	dur = result.dur;
	speed = result.speed;
	N = size(x1,2);
	gait2d = result.model.gait2d;		% function handle for musculoskeletal model
	symm = 0;							% no symmetry is assumed in the R3FAS project

	% Vector of time instants for nodes
	T = (1:N)*dur/N; 

	% Initialize the model (we need to run the model to compute GRFs and muscle forces
	result.model = initmodel(result.model);

	% if needed, create full gait cycle
	if (symm)
		x2 = x1(result.model.vmx,:);
		x2(1,:) = x2(1,:) + speed*dur;
		x = [x1 x2];
		dur = 2*dur;
	else
		x = x1;
	end
	actualfps = size(x,2)/dur;
	
	% initialize movie file
	moviefile = strrep(resultfile,'.mat','.avi');
	avi = avifile(moviefile, 'fps', 30, 'compression', 'Cinepak','quality',95);

	% initialize figure window
	close all
	figure(1);
	clf;
	set(gcf,'Position',[5 5 550 550]);
	set(gcf, 'color', 'white');
	
	% create ground points (random dots)
	np = 15000;
	xg = rand(np,1);
	yg = rand(np,1);
	xg = -1 + 3*[xg ; 2-xg];
	yg = -0.15*[yg ; yg];
	
	% make the movie
	R = [1:6 4];			% right stick points
	L = [2 7:10 8];			% left stick points
	nframes = size(x,2);
	shift = speed*dur/nframes;
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
			avi = addframe(avi,F);
			cla;
		end
	end
	avi = close(avi);
	hold off;
	close all
end