function stick(resultfile)

% plot a stick figure of a gait cycle

% Load the result
if (nargin == 0)
	[file,path] = uigetfile('*.mat');
	resultfile = [path file];
end
load(resultfile);
x = result.x;
u = result.u;
dur = result.dur;
speed = result.speed;
N = size(x,2);
gait2d = result.model.gait2d;		% function handle for musculoskeletal model

% Initialize the model (we need to run the model to compute GRFs and muscle forces
result.model = initmodel(result.model);

R = [1:6 4];			% right stick points
L = [2 7:10 8];			% left stick points

xmin = 1e6;
xmax = -1e6;
figure(1);clf;hold on
for i=1:5:size(x,2)
	d = result.model.gait2d('Stick',x(1:50,i));
	d(:,1) = d(:,1) + 0.1*(i-1);
	xmin = min(xmin,d(:,1));
	xmax = max(xmax,d(:,1));
	plot(d(R,1),d(R,2),'b',d(L,1),d(L,2),'r');
end
plot([xmin xmax],[0 0],'k');		% draw ground surface as a black line
axis('equal');
axis('off');
hold off;
