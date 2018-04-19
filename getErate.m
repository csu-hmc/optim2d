%======================================================================
%> @file getErate.m
%> @brief Matlab function to calculate the energy rate of a single time step
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

function [Edot, dEdot,w_ce, h_sl, h_am] = getErate(l_ceopt, FT, W, F_max, F_cepee, STIM, ACT, l_ce, Ahill, gmax, kPEE, PEEslack, v_ce)

%Variables
rho = 1059.7;   % Muscle density
sigma = 250e3;   % Max muscle stress
S = 1.5;        % For aerobic activity, S = 1.0 for anaerobic activity

%Calculate F_ce by substracting Fpee
k1 = 0.01*F_max;
k2 = kPEE*F_max/l_ceopt^2;
F_PEE = k1*(l_ce-PEEslack)+k2*(l_ce-PEEslack)^2*(l_ce > PEEslack);
F_ce = F_cepee-F_PEE;

% Calculate other variables
A = STIM*(STIM>ACT)+1/2*(STIM+ACT)*(STIM <= ACT);
F_iso = exp(-(l_ce-l_ceopt)^2/(W*l_ceopt)^2); % Force length relationship
g_vce = F_ce/ACT/F_max/F_iso; % Force-velocity relationship
mmass = F_max/sigma*rho*l_ceopt;

% Calculate velocity of contractile element
vbar_cemaxft = 12;  % Maximum shortening velocity
vbar_cemaxst = vbar_cemaxft/2.5; 
v_max = vbar_cemaxft * l_ceopt;
% Hack for negative muscle forces
if g_vce < 0;
    g_vce = 0;
end
if nargin < 13
    if g_vce < 1 %shortening
        v_ce = v_max*(g_vce-1)/(1+g_vce/Ahill);         % Use Hill curve
    else
        c = v_max*Ahill*(gmax-1)/(Ahill+1);
        v_ce = c*(1-g_vce)/(g_vce-gmax);                % Use Katz curve
    end
else
    if v_ce < 0
        g = (v_max+v_ce)/(v_max-v_ce/Ahill);
    else
        c = v_max*Ahill*(gmax-1)/(Ahill+1);
        g = (gmax*v_ce+c)/(v_ce+c);
    end
    F_ce = F_max*F_iso*ACT*g;
end
vbar_ce = v_ce/l_ceopt;

% Activation and maintenance
A_am = A^0.6;

% Nominal value
Nh_am = 25*(ACT<=(1-FT))+(128*FT+25)*(ACT>(1-FT));

if l_ce <= l_ceopt
    h_am = Nh_am*A_am*S;
else
    h_am = (0.4*Nh_am+0.6*Nh_am*F_iso)*A_am*S;
end

% Shortening-Lengthening energy
alpha_ST = 100/vbar_cemaxst;
alpha_FT = 153/vbar_cemaxft*(ACT > 1-FT); %Zero when activation is less than %ST fibers
alpha_L = 4*alpha_ST;

% Nominal value
if vbar_ce > 0
    Nh_sl = alpha_L*vbar_ce;
else
    Nh_sl = min([alpha_ST*vbar_cemaxst,-alpha_ST*vbar_ce*(1-FT)])-alpha_FT*vbar_ce*FT;
end
A_sl = A;
if vbar_ce <= 0
    A_sl = A_sl^2.0;
end
h_sl = Nh_sl*A_sl*S;
if l_ce > l_ceopt
    h_sl = h_sl*F_iso;
end

% Contractile element work
w_ce = -F_ce*v_ce/mmass;

Edot = max([1, h_am+h_sl+w_ce])*mmass; %Multiply by muscle mass to get W

% Derivatives
dFiso_dxi = zeros(1,8);
dFiso_dlce = -(l_ce-l_ceopt)^2/(W*l_ceopt)^2*F_iso*-2*l_ce/(W*l_ceopt)^2;
dFiso_dact = 0;
dA_dxi = zeros(1,8);
dA_dlce = 0;
dA_dact = 1/2*(STIM <= ACT);
dAam_dxi = zeros(1,8);
dAam_dlce = 0;
dAam_dact = 0.6*dA_dact^(-0.4);
if l_ce > l_ceopt
    dham_dxi = A_am*S*0.6*Nh_am*dFiso_dxi+(0.4*Nh_am+0.6*F_iso*Nh_am)*S*dAam_dxi;
    dham_dlce = A_am*S*0.6*Nh_am*dFiso_dlce+(0.4*Nh_am+0.6*F_iso*Nh_am)*S*dAam_dlce;
    dham_dact = A_am*S*0.6*Nh_am*dFiso_dact+(0.4*Nh_am+0.6*F_iso*Nh_am)*S*dAam_dact;  
else
    dham_dxi = Nh_am*S*dAam_dxi;
    dham_dact = Nh_am*S*dAam_dact;
    dham_dlce = Nh_am*S*dAam_dlce;
end
dForce_dxi = ones(1,8)*F_max*k1*-0.02;
dForce_dlce = -F_max*k1*l_ceopt;
dForce_dact = 0;
dFce_dxi = dForce_dxi;
dFce_dlce = dForce_dlce-k1-2*l_ce*k2*(l_ce>PEEslack);
dFce_dact = dForce_dact;
dgvce_dxi = 1/(ACT*F_max*F_iso)*dFce_dxi;
dgvce_dlce = 1/(ACT*F_max*F_iso)*dFce_dlce-F_ce/(ACT*F_max*F_iso^2)*dFiso_dlce;
dgvce_dact = 1/(ACT*F_max*F_iso)*dFce_dact-F_ce/(ACT^2*F_max*F_iso);
% derivative of vbar
if v_ce < 0
    dvce_dxi = 1/l_ceopt*(v_max/(1+g_vce/Ahill)*dgvce_dxi-(g_vce-1)/(1+g_vce/Ahill)^2*v_max/A*dgvce_dxi);
    dvce_dact = 1/l_ceopt*(v_max/(1+g_vce/Ahill)*dgvce_dact-(g_vce-1)/(1+g_vce/Ahill)^2*v_max/A*dgvce_dact);
    dvce_dlce = 1/l_ceopt*(v_max/(1+g_vce/Ahill)*dgvce_dlce-(g_vce-1)/(1+g_vce/Ahill)^2*v_max/A*dgvce_dlce);
else
    dvce_dxi = 1/l_ceopt*(-c/(g_vce-gmax)*dgvce_dxi-c*(1-g_vce)/(g_vce-gmax)^2*dgvce_dxi);
    dvce_dact = 1/l_ceopt*(-c/(g_vce-gmax)*dgvce_dact-c*(1-g_vce)/(g_vce-gmax)^2*dgvce_dact);
    dvce_dlce = 1/l_ceopt*(-c/(g_vce-gmax)*dgvce_dlce-c*(1-g_vce)/(g_vce-gmax)^2*dgvce_dlce);
end
if vbar_ce > 0
    dNhsl_dxi = alpha_L*dvce_dxi;
    dNhsl_dact = alpha_L*dvce_dact;
    dNhsl_dlce = alpha_L*dvce_dlce;
else
    dNhsl_dxi = -alpha_FT*FT*dvce_dxi;
    dNhsl_dact = -alpha_FT*FT*dvce_dact;
    dNhsl_dlce = -alpha_FT*FT*dvce_dlce;
    if alpha_ST*vbar_cemaxst < -alpha_ST*vbar_ce*(1-FT)
        dNhsl_dxi = dNhsl_dxi - alpha_ST*(1-FT)*dvce_dxi;
        dNhsl_dact = dNhsl_dact - alpha_ST*(1-FT)*dvce_dact;
        dNhsl_dlce = dNhsl_dlce - alpha_ST*(1-FT)*dvce_dlce;
    end
end
dAsl_dxi = dA_dxi*(vbar_ce>0)+2*dA_dxi*A*(vbar_ce<=0);
dAsl_dact = dA_dlce*(vbar_ce>0)+2*dA_dact*A*(vbar_ce<=0);
dAsl_dlce = dA_dact*(vbar_ce>0)+2*dA_dlce*A*(vbar_ce<=0);
if l_ce <= l_ceopt
    dhsl_dxi = dNhsl_dxi*A_sl*S+Nh_sl*dAsl_dxi*S;
    dhsl_dact = dNhsl_dact*A_sl*S+Nh_sl*dAsl_dact*S;
    dhsl_dlce = dNhsl_dlce*A_sl*S+Nh_sl*dAsl_dlce*S;
else
    dhsl_dxi = dNhsl_dxi*A_sl*S*F_iso+Nh_sl*dAsl_dxi*S*F_iso+Nh_sl*A_sl*S*dFiso_dxi;
    dhsl_dact = dNhsl_dact*A_sl*S*F_iso+Nh_sl*dAsl_dact*S*F_iso+Nh_sl*A_sl*S*dFiso_dact;
    dhsl_dlce = dNhsl_dlce*A_sl*S*F_iso+Nh_sl*dAsl_dlce*S*F_iso+Nh_sl*A_sl*S*dFiso_dlce;
end
dwce_dxi = -F_ce/mmass*dvce_dxi-vbar_ce/mmass*dFce_dxi;
dwce_dact = -F_ce/mmass*dvce_dact-vbar_ce/mmass*dFce_dact;
dwce_dlce = -F_ce/mmass*dvce_dlce-vbar_ce/mmass*dFce_dlce;

if Edot == 1
    dEdot_dxi = 0;
else
    dEdot_dxi = mmass*(dwce_dxi+dhsl_dxi+dham_dxi);
    dEdot_dact = mmass*(dwce_dact+dhsl_dact+dham_dact);
    dEdot_dlce = mmass*(dwce_dlce+dhsl_dlce+dham_dlce);
end
dEdot = [dEdot_dxi dEdot_dlce dEdot_dact]';
end