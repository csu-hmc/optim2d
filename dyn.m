%==================================================================
function [f,dfdx,dfdxdot,dfdu] = dyn(model,x,xdot,u)
	% evaluates the implicit dynamics equation: f(x,xdot,u)=0
	%
	% state x is the state of the gait2d model, augmented by four state variables
	% of the CCF prosthetic knee model: s, v1, v2, M
	% control u is the control of the gait2d model, plus two valve controls u1,u2

	% this function can evaluate the following versions of the model (model.type):
	% 'zero'		A trivial model that has zeros as its solution
	% 'able'		The gait2d model, with actuator disconnected
	% 'bka'			As 'able', but on right side, no ankle muscles and a passive prosthetic foot
	% 'free'		As 'able', but on right side, no knee/ankle muscles
	% 'cleg'		As 'free', but actuator connected and C1max=0 (valve 1 closed)
	% 'ccf'			As 'cleg', but all valves operating

	nxg = 50;		% number of states for musculoskeletal model
	nug = 16;		% number of controls for musculoskeletal model
	nstates = nxg + 4;
	ncontrols = nug + 2;
	ixg = 1:nxg;
	iug = 1:nug;

	% error checking
	if size(x,1) ~= nstates
		error('optim:dyn error: incorrect size for x, size = %d %d', size(x,1), size(x,2));
	end
	if size(xdot,1) ~= nstates
		error('optim:dyn error: incorrect size for xdot, size = %d %d', size(xdot,1), size(xdot,2));
	end
	if size(u,1) ~= ncontrols
		error('optim:dyn error: incorrect size for u, size = %d %d', size(u,1), size(u,2));
	end

	% x and u indices for the prosthesis variables
	ikneevel = 14;		% knee angular velocity
	is = nxg+1;
	iv1 = nxg+2;
	iv2 = nxg+3;
	iM = nxg+4;
	iu1 = nug+1;
	iu2 = nug+2;
	
	% preallocate the outputs
	nx = size(x,1);
	nu = size(u,1);
	f = zeros(nx,1);
	dfdx 	= spalloc(nx, nx, model.nnz_dfdx);
	dfdxdot = spalloc(nx, nx, model.nnz_dfdxdot);
	dfdu 	= spalloc(nx, nu, model.nnz_dfdu);
	
	if strcmp(model.type, 'zero')
		f = x;
		f(2) = x(2) - 1.0;			% this drives the vertical hip position to 1, avoiding high GRF
		dfdx = speye(nx);
	else
		% first evaluate the musculoskeletal part of the model
		if strcmp(model.type, 'AFO')
            % Add stiffness to ankle
            M = [0 x(iM) (-x(6)*model.anklestiffness-x(15)*model.ankledamping) 0 0 0]';			% apply the right knee moment from the prosthesis
%             dMdx6 = [0 0 (x(6)<0)*(-model.anklestiffness) 0 0 0]; (x(6)<0)*
%             dMdx15 = [0 0 (x(6)<0)*(-model.ankledamping) 0 0 0];
        else
            M = [0 x(iM) 0 0 0 0]';			% apply the right knee moment from the prosthesis
        end
		[f(ixg), dfdx(ixg,ixg), dfdxdot(ixg,ixg), dfdu(ixg,iug), dfdM] = model.gait2d('Dynamics',x(1:nxg),xdot(1:nxg),u(1:nug),M);
		dfdx(ixg,iM) = dfdM(:,2);
%         dfdx(ixg,6) = dfdM*dMdx6';
%         dfdx(ixg,15) = dfdM*dMdx15';
	
		% if there is no prosthetic knee in the model, use simple equations for the prosthetic variables
		if strcmp(model.type, 'able') || strcmp(model.type, 'bka') || strcmp(model.type, 'AFO')
			% s = 0
			f(nxg+1) = x(is);
			dfdx(nxg+1,is) = 1;
			% v1 = 0;
			f(nxg+2) = x(iv1);
			dfdx(nxg+2,iv1) = 1;
			% v2 = 0;
			f(nxg+3) = x(iv2);
			dfdx(nxg+3,iv2) = 1;
			% M = 0;
			f(nxg+4) = x(iM);
			dfdx(nxg+4,iM) = 1;
			
		% if there is a prosthetic knee in the model, use the four hydraulic system equations
		elseif strcmp(model.type, 'aka')
			scale = 1/2500;			% valve equations must be scaled down because Cmax^2 is a large number
			if (model.C1max ~= 0.0)
				% these are the actual equations for the CCF knee energy storage, controlled by valve 1
				
				% state equation for accumulator
				f(nxg+1) = xdot(is) + x(iv1);
				dfdxdot(nxg+1,is) = 1;
				dfdx(nxg+1,iv1) = 1;
				
				% valve equation for valve 1: u1^2 * C1max^2*(M - RA*k*s + RA*B1*v1) + RA * |v1| * v1 = 0;
				f(nxg+2) = scale * ( u(iu1)^2 * model.C1max^2 * ( x(iM) - model.RA*model.k*x(is) + model.RA*model.B1*x(iv1) ) ...
					+ model.RA * abs(x(iv1)) * x(iv1) );
				dfdx(nxg+2,iM) = scale * u(iu1)^2 * model.C1max^2;
				dfdx(nxg+2,is) = -scale * u(iu1)^2 * model.C1max^2 * model.RA * model.k;
				dfdx(nxg+2,iv1) = scale * ( u(iu1)^2 * model.C1max^2 * model.RA * model.B1 + 2*model.RA*x(iv1)*sign(x(iv1)) );
				dfdu(nxg+2,iu1) = scale * 2*u(iu1) * model.C1max^2 * ( x(iM) - model.RA*model.k*x(is) + model.RA*model.B1*x(iv1) );
			else
				% The following just drives accumulator volume (s) and valve flow (v1) to zero, which is OK
				% for C-leg and others in which there is no valve 1.  In principle, we could use the
				% equations above, but that is numerically not as stable (only sdot is constrained, and
				% v1 occurs only quadratically so is not driven hard enough to zero).
				
				% s = 0
				f(nxg+1) = x(is);
				dfdx(nxg+1,is) = 1;
				
				% v1 = 0;
				f(nxg+2) = x(iv1);
				dfdx(nxg+2,iv1) = 1;
			end
		
			% actuator equation: RA*phidot - v1 - v2 + L*M = 0
			f(nxg+3) = model.RA * x(ikneevel) - x(iv1) - x(iv2) + model.L*x(iM) ;
			dfdx(nxg+3, ikneevel) = model.RA;
			dfdx(nxg+3,iv1) = -1;
			dfdx(nxg+3,iv2) = -1;
			dfdx(nxg+3, iM) = model.L;
			
			% valve equation for valve 2: u2^2 * C2max^2 * (M + RA*B2*v2) + RA * |v2| * v2 = 0;
			f(nxg+4) = scale * (u(iu2)^2 *  model.C2max^2 * ( x(iM) + model.RA*model.B2*x(iv2) ) ...
				+ model.RA * abs(x(iv2)) * x(iv2));
			dfdx(nxg+4,iM) = scale * u(iu2)^2 * model.C2max^2;			
			dfdx(nxg+4,iv2) = scale * (u(iu2)^2 * model.C2max^2 * model.RA * model.B2 + 2 * model.RA * x(iv2) * sign(x(iv2)) );
			dfdu(nxg+4,iu2) = scale * 2*u(iu2) * model.C2max^2 * (x(iM) + model.RA*model.B2*x(iv2));
		else
			error('optim:dyn error: unknown model type: %s', model.type);
		end
	end

end
