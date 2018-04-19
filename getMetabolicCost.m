%======================================================================
%> @file getMetabolicCost.m
%> @brief Matlab function to calculate Metabolic cost of a gait cycle
%>
%> @author Anne Koelewijn
%> @date March 23, 2015
%>
%> @details
%> Function to calculate the metabolic cost of a gait cycle using
%Umberger's model. Ross Miller's code was also used as a reference
%> 
%>Input is a Matlab structure (\p problem) and double (\p X). If the
%function is called as postprocessing, you need the third input (\p init)
%to be anything
%> @param problem describes the problem that is optimized
%> @param X describes the current state of the problem
%> @param init is used to initialize the problem if necessary

%>
%> The output of getMetabolicCost.m are three values and their derivatives with respect to the state.
%/p MetCost is equal to the Metabolic cost %of the single gait cycle in J/kg/m units, /p MetRate is
%equal to the Metabolic rate of the single gait cycle in W/kg and /p CoT is equal to the cost of transport
%of the gait cycle, CoT = MetCost/g
%> @retval MetCost Metabolic Cost of gait cycle
%> @retval dMetCostdx Derivative of MetCost wrt state
%> @retval MetRate Metabolic rate during gait cycle
%> @retval dMetRatedx Derivative of MetRate wrt state
%> @retval CoT Cost of transport during gait cycle
%> @retval dCotdx Derivative of CoT wrt state
%======================================================================

function [MetCost,MetRate, CoT] = getMetabolicCost(result, init) % 

model = result.model;

if nargin == 2
    % Initialize the model (we need to run the model to compute GRFs and muscle forces)
    model = initmodel(model);
    problem.ndof = 9;
	problem.nmus = 16;
	problem.nstates = 2*problem.ndof + 2*problem.nmus + 4;
	problem.ncontrols = problem.nmus + 2;
    problem.N = size(result.x,2);
end

% Variables
speed = result.speed;                  % Current walking speed
mass = gait2d('Get', 'Total Mass');    % Mass of the subject
g = 9.81;     % Gravitational constant
T = result.dur; % Duration of gait cycle
noMuscles = problem.nmus;        % Number of muscles
N = problem.N;                   % Number of collocation points
h = T/N;                         % Duration of time step

% Initialize paremeters
MetRateperm = zeros(noMuscles,1);    
Edot = zeros(noMuscles,N);

% Get muscle variables
l_ceopt =  gait2d('Get','Lceopt');             % Optimal fiber length
FT = gait2d('Get','FT');                       % Percentage of fast twitch fibers
W = gait2d('Get', 'Width');                    % Width of hill curve
F_max = gait2d('Get', 'Fmax');                 % Max muscle force
Ahill = gait2d('Get', 'Ahill');                % Hill constant
Gmax = gait2d('Get', 'Gmax');                  % Maximum eccentric force
kPEE = gait2d('Get', 'kPEE');                  % PEE stiffness parameter
PEEslack = gait2d('Get', 'PEEslack').*l_ceopt; % PEE slack length

% dMetRatedx = zeros(size(result.x));
ilce = (problem.ndof*2+1):(problem.ndof*2+noMuscles); % Lce comes after velocity
iact = (problem.ndof*2+noMuscles+1):(problem.ndof*2+noMuscles*2); % Activation comes after Lce
for i = 1:N
    x = result.x(:,i);
    STIM = result.u(1:noMuscles,i); % Input is stimulatiom
    ACT = x(iact);
    l_ce = x(ilce).*l_ceopt; %Lce is given as percentage of L_ceopt
    F_ce = gait2d('Muscleforces', x(1:50)); % Force in the muscles

    %calculate v_ce
    if i > 1 % Take node i+1
        l_ceprev = result.x(ilce,i-1).*l_ceopt;
        v_ce = (l_ce-l_ceprev)/h;
    else
        l_ceprev = result.x(ilce,N).*l_ceopt;
        v_ce = (l_ce-l_ceprev)/h;
    end
    % Calculate energy rate in Watts for all muscles at once, use
    % discontinuous model in postprocessing
    if nargin == 2 %initialization means postprocessing
        [Edot(:,i)] = getErate_discont(l_ceopt, FT, W, F_max, F_ce, STIM, ACT, l_ce, Ahill, Gmax, kPEE, PEEslack, v_ce);
        dEdot = zeros(problem.nmus,problem.ndof*2+4); %just here to make the program work
    else
        [Edot(:,i),dEdot] = getErate(l_ceopt, FT, W, F_max, F_ce, STIM, ACT, l_ce, Ahill, Gmax, kPEE, PEEslack, dFdx, v_ce,h,T);
    end
            
    % No Edot in bka muscles
    if strcmp(result.model.type, 'bka')
        % No energy expenditure in lower right leg in amputee
        Edot([6,7,8],i) = 0;
    end
    
%     %Put derivatives in matrix correctly
%     dMetRatedx(ix(1:problem.ndof*2)) = dMetRatedx(ix(1:problem.ndof*2))+sum(dEdot(:,1:problem.ndof*2))'; % dang
%     dMetRatedx(ilce) = dMetRatedx(ilce) + dEdot(:,problem.ndof*2+1); % dlce
%     dMetRatedx(iact) = dMetRatedx(iact) + dEdot(:,problem.ndof*2+2); % dact
%     dMetRatedx(iu) = dMetRatedx(iu) + dEdot(:,problem.ndof*2+4); %du
%     dMetRatedx(end-1) = dMetRatedx(end-1)+sum(dEdot(:,end)); %dT
%     if i > 1
%         dMetRatedx(ilce-nvarpernode) = dMetRatedx(ilce-nvarpernode)+dEdot(:,problem.ndof*2+3);
%     else
%         dMetRatedx(ilce+N*nvarpernode) = dMetRatedx(ilce+N*nvarpernode)+dEdot(:,problem.ndof*2+3);
%     end
end

for i = 1:noMuscles
    % Calculate metabolic rate in Watts
    MetRateperm(i) = 1/N*sum(Edot(i,:));
end

% Metabolic Rate in W/kg
MetRate = 1/mass*sum(MetRateperm);
% dMetRatedx = 1/mass*dMetRatedx;
MetRate = MetRate+1.20; %1.20 W/kg is quiet standing for a men, 1.13 for women

%Determine metabolic cost and cost of transport
MetCost = MetRate/speed;
% dMetCostdx = 1/speed*dMetRatedx;
% dMetCostdx(end) = MetRate*(-1/speed^2);
CoT = MetCost/g;
% dCotdx = dMetCostdx/g;