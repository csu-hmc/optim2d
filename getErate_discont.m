%======================================================================
%> @file getErate_discont.m
%> @brief Matlab function to calculate the energy rate of a single time
%step for post processing. This is the discontinuous version of the code
%>
%> @author Anne Koelewijn
%> @date March 23, 2015
%>
%> @details
%> Function to calculate the energy rate at a single time step using
%Umberger's model. Ross Miller's code was also used as a reference
%> 
%>Input is are muscle parameters of the specific muscle that is considered.
%> @param l_ceopt Optimal fiber length
%> @param FT Percentage of fast twitch fiber in the muscle
%> @param W Width of the force-length curve
%> @param F_max Maximum muscle force
%> @param F_ce Force in the contractile element
%> @param STIM stimulation of the muscle (between 0 and 1)
%> @param ACT activation level of the muscle (between 0 and 1)
%> @param l_ce current lenght of the contractile element
%> @param Ahill Hill constant in the force-velocity curve
%> @param gmax Maximum muscle eccentric velocity
%>
%> The output of getErate.m is a double which is equal to the energy rate
%of the time step in W/kg
%> @retval Edot is the energy expenditure during the current time step
%======================================================================

function [Edot,w_ce, h_sl, h_am] = getErate_discont(l_ceopt, FT, W, F_max, F_cepee, STIM, ACT, l_ce, Ahill, gmax, kPEE, PEEslack, v_ce)

%Variables
rho = 1059.7;   % Muscle density
sigma = 250e3;   % Max muscle stress
S = 1.5;        % For aerobic activity, S = 1.0 for anaerobic activity

%Calculate F_ce by substracting Fpee
k1 = 0.01*F_max;
k2 = kPEE.*F_max./l_ceopt.^2;
F_PEE = k1.*(l_ce-PEEslack)+k2.*(l_ce-PEEslack).^2.*(l_ce > PEEslack);
F_ce = F_cepee-F_PEE;

% Calculate other variables
A = STIM.*(STIM>ACT)+1/2.*(STIM+ACT).*(STIM <= ACT);
F_iso = exp(-(l_ce-l_ceopt).^2./(W.*l_ceopt).^2); % Force length relationship
mmass = (F_max/sigma*rho).*l_ceopt;

% Calculate velocity of contractile element
vbar_cemaxft = 12;  % Maximum shortening velocity
vbar_cemaxst = vbar_cemaxft/2.5; 
v_max = vbar_cemaxft * l_ceopt;
if isempty(F_ce)
    if v_ce < 0
        g = (v_max+v_ce)/(v_max-v_ce./Ahill);
    else
        c = v_max.*Ahill.*(gmax-1)./(Ahill+1);
        g = (gmax*v_ce+c)./(v_ce+c);
    end
    F_ce = F_max.*F_iso.*ACT.*g;
end
vbar_ce = v_ce./l_ceopt;

% Activation and maintenance
A_am = A.^0.6;

% Nominal value
Nh_am = 25.*(ACT<=(1-FT))+(128*FT+25).*(ACT>(1-FT));
h_am = (Nh_am.*A_am*S).*(l_ce <= l_ceopt) + ((0.4*Nh_am+0.6*Nh_am.*F_iso).*A_am*S).*(l_ce > l_ceopt);

% Shortening-Lengthening energy
alpha_ST = 100/vbar_cemaxst;
alpha_FT = 153/vbar_cemaxft*(ACT > 1-FT); %Zero when activation is less than %ST fibers
alpha_L = 4*alpha_ST;

% Nominal value
Nh_sl = alpha_L.*vbar_ce.*(vbar_ce > 0) + (100*(alpha_ST*vbar_cemaxst<-alpha_ST*vbar_ce.*(1-FT))-alpha_ST*vbar_ce.*(1-FT).*(alpha_ST*vbar_cemaxst > -alpha_ST*vbar_ce.*(1-FT))-alpha_FT.*vbar_ce.*FT).*(vbar_ce <= 0);

A_sl = A.*(vbar_ce > 0)+ A.^(2.0).*(vbar_ce <=0);
h_sl = Nh_sl.*A_sl*S;
h_sl = h_sl.*F_iso.*(l_ce > l_ceopt) + h_sl.*(l_ce <= l_ceopt);

% Contractile element work
w_ce = max(0,-F_ce.*v_ce./mmass);

Edot = (1*(1> h_am+h_sl+w_ce)+ (h_am+h_sl+w_ce).*(1 <= h_am+h_sl+w_ce)).*mmass; %Multiply by muscle mass to get W