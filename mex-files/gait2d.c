/*=================================================================
 *
 * gait2d.c
 *
 * Implicit differential equation for 2D musculoskeletal model : f(x,dx/dt,u) = 0
 * This is the source code for the MEX function gait2d.mexw32
 * The musculoskeletal model is documented in the file gait2d.doc
 * The function documentation is in gait2d.m
 *
 * Copyright 2009-2011 Orchard Kinetics LLC
 *
 *=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "gait2d.h"

// Compilation settings
#define VERBOSE 0					// set this to 1 for debugging

// size of the model (some other size constants are in gait2d.h)
#define NMUS 16						// number of muscles 
#define NSTATES 2*NDOF+2*NMUS		// number of system state variables, 2*NDOF + 2*NMUS

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// define struct that holds muscle properties
typedef struct {
	double Lceopt; 			// Optimal length of CE (m)
	double Width;			// Width of CE force-length relationship (m)
	double Fmax;			// Maximal isometric force of CE (N)
	double Vmax;			// Max. contraction velocity in Lceopt/s
	double Tact, Tdeact;	// activation and deactivation time constants
	double gmax;			// maximum eccentric force
	double SEEslack;		// Slack length of the SEE, in meters
	double PEEslack;		// Slack length of the PEE, relative to Lceopt
	double umax;			// strain of SEE at Fmax load
	double kPEE;			// Stiffness parameter of PEE, in Fmax/Lceopt^2
	double L0;				// Muscle+tendon length when all DOFs are zero
	double MA[NMOM];		// Moment arms (positive when muscle causes ccw rotation of distal joint, model facing +X)
	// the following parameters are derived from other parameters during initialization
	double c1,c2;			// Activation-deactivation rate constants
	double c3;				// Continuity parameter for eccentric force-velocity relationship
	double kSEE;			// Stiffness parameter of SEE, in Fmax/m^2
} muscleprop;

// global variables
static int dynamics_evaluations = 0;
static int initialized = 0;
static double zeros[NSTATES];			// contains zeros
static muscleprop muscles[NMUS];		// contains all muscle properties
static param_struct param;				// contains other model parameters
static int nonzeromomentarms;			// number of nonzero moment arms (nonzeros in dL/dq), for memory allocation

// ===================================================================================
// ExtractParameters: extract model parameters from the parameter matrix
// ===================================================================================
void ExtractParameters(double par[]) {
	// The parameter matrix is obtained by: readxls('gait2d_par.xls').
	// par(1,1) is the first numerical element in the XLS file.
	// par contains NaN for the text cells in the XLS file	
	int i,j,s;
	#define MAXSECTIONS 10
	int sectionrows[MAXSECTIONS];
	
	// Locate sections of parameter matrix, and store them in proper place here
	int nrows = 0;
	int label;
	while ((label = (int) par[nrows++]) != 999) {
		if (label >= 0 && label < MAXSECTIONS) sectionrows[label] = nrows-1;
	};
	// Read general parameters
	s = sectionrows[1] + 1 + nrows;
	param.ContactStiff	= par[s++];
	param.ContactDamp	= par[s++];
	param.ContactY0		= par[s++];
	param.ContactFric	= par[s++];
	param.ContactV0		= par[s++];
	param.HillA			= par[s++];
	param.Gmax			= par[s++];
	
	// Read multibody model parameters
	s = sectionrows[2] + 2 + nrows;
	param.TrunkMass 	= par[s++];
	param.ThighMass 	= par[s++];
	param.ShankMass 	= par[s++];
	param.FootMass 		= par[s++];
	s = sectionrows[2] + 2 + 2*nrows;
	param.TrunkInertia 	= par[s++];
	param.ThighInertia 	= par[s++];
	param.ShankInertia 	= par[s++];
	param.FootInertia 	= par[s++];
	s = sectionrows[2] + 5 + 3*nrows;
	param.FootCMx		= par[s++];
	s = sectionrows[2] + 2 + 4*nrows;
	param.TrunkCMy 		= par[s++];
	param.ThighCMy 		= par[s++];
	param.ShankCMy 		= par[s++];
	param.FootCMy 		= par[s++];
	s = sectionrows[2] + 3 + 5*nrows;
	param.ThighLen 		= par[s++];
	param.ShankLen 		= par[s++];
	
	// Read muscle parameters
	for (i=0; i<NMUS; i++) {
		s = sectionrows[3] + 2 + i + nrows;
		muscles[i].Fmax 	= par[s]; s = s + nrows;
		muscles[i].Lceopt 	= par[s]; s = s + nrows;
		muscles[i].Width 	= par[s]; s = s + nrows;
		muscles[i].PEEslack	= par[s]; s = s + nrows;
		muscles[i].SEEslack	= par[s]; s = s + nrows;
		muscles[i].kPEE		= par[s]; s = s + nrows;
		muscles[i].umax		= par[s]; s = s + nrows;
		muscles[i].Vmax		= par[s]; s = s + nrows;
		muscles[i].Tact		= par[s]; s = s + nrows;
		muscles[i].Tdeact	= par[s]; s = s + nrows;
		muscles[i].L0		= par[s]; s = s + nrows;
		for (j=0; j<NMOM; j++) {
			muscles[i].MA[j] = par[s]; s = s + nrows;	
		}
	}

	// Determine the number of nonzero moment arms, for memory allocation
	nonzeromomentarms = 0;
	for (i=0; i<NMUS; i++)
		for (j=0; j<NMOM; j++)
			if (muscles[i].MA[j] != 0) nonzeromomentarms++;
			
	// Read contact model parameters
	s = sectionrows[4] + 2 + nrows;
	param.ContactHeelX	= par[s++];
	param.ContactToeX	= par[s++];
	s = sectionrows[4] + 2 + 2*nrows;
	param.ContactY		= par[s++];
	
	// Read joint range of motion parameters
	for (i=0; i<NMOM; i++) {
		s = sectionrows[5] + 2 + i+ nrows;
		param.MinAngle[i] 	= par[s];
		param.MaxAngle[i] 	= par[s+nrows];
		param.JointK[i]		= par[s+2*nrows];
		param.JointD0[i]	= par[s+3*nrows]*M_PI/180;
		param.JointB[i]		= par[s+4*nrows];
	}
	
#if VERBOSE
	printf("\nGeneral parameters:\n");
	printf("	HillA			%8.4f\n", param.HillA);          
	printf("	Gmax			%8.4f\n", param.Gmax);          

	printf("\nBody segment parameters:\n");
	printf("		Mass			Inertia			CMx			CMy			Length\n");
	printf("Trunk	%7.3f			%7.3f			zero		%7.3f		n/a\n", param.TrunkMass,param.TrunkInertia,param.TrunkCMy);
	printf("Thigh	%7.3f			%7.3f			zero		%7.3f		%7.3f\n", param.ThighMass,param.ThighInertia,param.ThighCMy,param.ThighLen);
	printf("Shank	%7.3f			%7.3f			zero		%7.3f		%7.3f\n", param.ShankMass,param.ShankInertia,param.ShankCMy,param.ShankLen);
	printf("Foot	%7.3f			%7.3f			%7.3f		%7.3f		n/a\n", param.FootMass,param.FootInertia,param.FootCMx,param.FootCMy);

	printf("\nMuscle parameters:\n");
	printf("muscle      Fmax   Lceopt  Width   PEEslack  SEEslack  kPEE   umax   Vmax   Tact  Tdeact\n");
	printf("----------------------------------------------------------------------------------------\n");
	for (i=0; i<NMUS; i++) {
		printf("%5d    %7.2f %7.4f %7.3f %7.3f %7.4f %7.3f %7.3f %7.2f %7.3f %7.3f\n", i+1,
			muscles[i].Fmax, muscles[i].Lceopt, muscles[i].Width, muscles[i].PEEslack, 
			muscles[i].SEEslack, muscles[i].kPEE, muscles[i].umax, muscles[i].Vmax, 
			muscles[i].Tact, muscles[i].Tdeact);
	}
	printf("----------------------------------------------------------------------------------------\n");
	printf("\nMuscle path models:\n");
	printf("muscle      L0	   dRhip	dRknee	dRankle	dLhip	dLknee	dLankle\n");
	printf("------------------------------------------------------------------\n");
	for (i=0; i<NMUS; i++) {
		printf("%5d    %7.3f", i+1, muscles[i].L0);
		for (j=0; j<NMOM; j++)
			if (muscles[i].MA[j] == 0)
				printf("        ");
			else
				printf(" %7.3f", muscles[i].MA[j]);
		printf("\n");
	}
	printf("\n");
	printf("------------------------------------------------------------------\n");
	
	printf("\nGround contact parameters:\n");
	printf("	Stiffness		%8.1f N/m\n", 	param.ContactStiff); 
	printf("	Damping			%8.3f\n", 		param.ContactDamp); 
	printf("	ContactY		%8.3f m\n", 	param.ContactY);
	printf("	ContactY0		%8.3f m\n", 	param.ContactY0);
	printf("	ContactFric		%8.3f\n", 		param.ContactFric);
	printf("	ContactV0		%8.3f m/s\n", 	param.ContactV0); 
	printf("	ContactX heel	%8.3f m\n", 	param.ContactHeelX); 
	printf("	ContactX toe	%8.3f m\n", 	param.ContactToeX); 
	
	printf("\nJoint range of motion parameters:\n");
	printf("Joint    Min.angle   Max.angle   JointK   JointD0   JointB\n");
	printf("--------------------------------\n");
	for (i=0; i<NMOM; i++) {
		printf("%5d    %8.2f  %8.2f\n", i,param.MinAngle[i], param.MaxAngle[i], param.JointK[i], param.JointD0[i], param.JointB[i]);
	}
	
	printf("--------------------------------\n");
#endif

// Preprocessing and error checking
	for (i=0; i<NMOM; i++) {
		param.MinAngle[i] = param.MinAngle[i]*M_PI/180.0;
		param.MaxAngle[i] = param.MaxAngle[i]*M_PI/180.0;
		if (param.MinAngle[i] > param.MaxAngle[i]) {
			printf("Error in joint %d\n", i);
			mexErrMsgTxt("Max angle must be greater than min angle.");
		}
	}
	for (i=0; i<NMUS; i++) {
		muscles[i].kSEE = 1.0/(muscles[i].umax*muscles[i].umax*muscles[i].SEEslack*muscles[i].SEEslack);
		muscles[i].c2 = 1.0/muscles[i].Tdeact;
		muscles[i].c1 = 1.0/muscles[i].Tact - muscles[i].c2;
		if (muscles[i].c2 < 0 || muscles[i].c1 < 0) {
			printf("Error in muscle %d\n", i);
			mexErrMsgTxt("Muscle time constants must be positive, and deactivation > activation.");
		}
		muscles[i].gmax = param.Gmax;		// for now, gmax is same for all muscles
		muscles[i].c3 = muscles[i].Vmax * param.HillA * (muscles[i].gmax - 1.0) / (param.HillA + 1.0);
	}
}

// ================================================================================
//  Softpos: soft positive function y = Softpos(x) and its derivative
// ================================================================================
double Softpos(double x, double x0, double *dydx) {
	double s,y;
	s = sqrt(x*x + x0*x0);
	y = 0.5*(s + x);
	*dydx = 0.5*(x/s + 1.0); 
	return y;
}

// ===================================================================================
// MusclePath: calculate muscle-tendon length and its derivatives w.r.t. joint angles
// ===================================================================================
double MusclePath(muscleprop *m, double ang[NMOM], double dLm_dang[NMOM]) {

	// this is the 2D model with constant moment arms (which are stored in MA)
	// in the 3D model we will have polynomial models for Lm(ang) and
	
	int i;
	double Lm;
	
	Lm = m->L0;
	for (i=0; i<NMOM; i++) {
		Lm = Lm - m->MA[i]*ang[i];
		dLm_dang[i] = -m->MA[i];
	}
	return Lm;
}


// ===================================================================================
// MuscleDynamics: the implicit muscle dynamics, returns f(a,Lce,Lcedot,Lm) and its derivatives
// ===================================================================================
double MuscleDynamics(muscleprop *m,										// muscle parameters (input)
	double a, double Lce, double Lcedot, double Lm,							// the input variables
	double *df_da, double *df_dLce, double *df_dLcedot, double *df_dLm, 	// the gradients (output)													
	double *force, double *dforce_dLm, double *dforce_dLce,					// muscle force (output) and derivatives
	double *powerCE) {														// power generated by CE (only valid if Lcedot satisfies dynamics)

	double F1, dF1_dLce;
	double f,x,k1;
	double F2, dF2_dLcedot;
	double Fpee, dFpee_dLce;
	double Fsee, dFsee_dLce, dFsee_dLm;
	double Fdamp, dFdamp_dLcedot;
	double Fce;
	
	// F1 is the isometric force-length relationship at max activation, normalized to Fmax
	x = (Lce - 1.0)/m->Width;
	F1 = exp(-x*x);
	dF1_dLce = -2.0*x*F1 / m->Width;
	
	// F2 is the normalized force-velocity relationship
	if (Lcedot < 0) {
		// concentric contraction
		x = (m->Vmax - Lcedot/param.HillA);
		F2 = (m->Vmax + Lcedot)/x;
		dF2_dLcedot = (1.0 + F2/param.HillA)/x;
	}
	else {
		// eccentric contraction
		x = Lcedot + m->c3;
		F2 = (m->gmax*Lcedot + m->c3) / x;
		dF2_dLcedot = (m->gmax - F2)/x;
	}
	
	// Fdamp is a small viscous damping of CE (0.001 of Fmax at 1 Lceopt/s) to ensure df/dLcedot is never zero
	Fdamp = 0.001*Lcedot;
	dFdamp_dLcedot = 0.001;
	
	// Fpee is the PEE force-length relationship, normalized to Fmax
	k1 = 0.01 * m->Lceopt;			// stiffness of the linear term is 0.01 Fmax/meter
	x = (Lce - m->PEEslack);		// elongation of PEE, in Lceopt units
	Fpee = k1*x;					// low stiffness linear term
	dFpee_dLce = k1;
	if (x>0) {						// add quadratic term for positive deformation
		Fpee = Fpee + m->kPEE*x*x;
		dFpee_dLce = dFpee_dLce + 2*m->kPEE*x;
	}
	
	//  Fsee is the SEE force-length relationship, normalized to Fmax
	k1 = 0.01;										// stiffness of the linear term is 0.01 Fmax/meter	
	x = Lm - Lce*m->Lceopt - m->SEEslack;			// elongation of SEE, in meters
	Fsee = k1*x;									// low stiffness linear term
	dFsee_dLce = -k1*m->Lceopt;
	dFsee_dLm = k1;
	if (x>0) {										// add quadratic term for positive deformation
		Fsee = Fsee + m->kSEE*x*x;
		dFsee_dLce = dFsee_dLce - 2*m->kSEE*m->Lceopt*x;
		dFsee_dLm = dFsee_dLm + 2*m->kSEE*x;
	}

	// Compute f, the force imbalance in the muscle contraction, and its derivatives
	Fce = a*F1*F2 + Fdamp;
	f = Fsee - Fce - Fpee;
	*df_da = -F1*F2;
	*df_dLce = dFsee_dLce - a*dF1_dLce*F2 - dFpee_dLce;
	*df_dLcedot = -a*F1*dF2_dLcedot - dFdamp_dLcedot;
	*df_dLm = dFsee_dLm;
	
	// Muscle force is the force in SEE
	*force = m->Fmax*Fsee;
	*dforce_dLm = m->Fmax * dFsee_dLm;
	*dforce_dLce = m->Fmax * dFsee_dLce;
	
	// power (W) generated by CE (positive when shortening)
	// Fce was in Fmax units, and Lce was in Lceopt units, so some scaling needed
	*powerCE = -m->Fmax*Fce*Lcedot*Lce*m->Lceopt;
	
	// Return the imbalance
	return f;
}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// working variables
	int i, j, k, nrows, ncols;
	char *command, *parametername;
	double *par, ang, angvel, dmax, dmin, ddmax, ddmin;
	double massfactor, lengthfactor, inertiafactor;
	int derivatives;
	
	// muscle variables
	double Lm[NMUS];				// muscle+tendon length, based on skeleton kinematic state
	double dLm_dang[NMUS][NMOM];	// derivatives of Lm with respect to joint angles
	double force[NMUS];				// muscle forces
	double CEpower[NMUS];			// power generated by CE of each muscle
	double g[NMUS];					// muscle force imbalance
	double dg_da[NMUS], dg_dLce[NMUS], dg_dLcedot[NMUS], dg_dLm[NMUS]; 	// derivatives of muscle imbalance
	double h[NMUS];					// activation dynamics imbalance
	double dh_da[NMUS], dh_dadot[NMUS], dh_du[NMUS];
	double dforce_dLm[NMUS], dforce_dLce[NMUS];

	// multibody dynamics variables
	double *q, *qd, *qdd;
	double mom[NMOM];
	double zero[NDOF];
	double dz_dq[NDOF][NDOF];
	double dz_dqd[NDOF][NDOF];
	double dz_dqdd[NDOF][NDOF];
	double dz_dmom[NDOF][NMOM];
	double dmom_dang[NMOM][NMOM];
	double dmom_dangvel[NMOM];		// assumed to be a diagonal matrix (no coupling between joints)
	double dmom_dLce[NMOM][NMUS];
	double GRF[4];
	double dGRF_dq[4][NDOF];
	double dGRF_dqd[4][NDOF];
	double Stick[NSTICK][2];
	double tmp[16];

	// accelerometer model variables
	double acc[42];
	double dacc_dq[42][NDOF];
	double dacc_dqd[42][NDOF];
	double dacc_dqdd[42][NDOF];
	
	// MEX function pointers to inputs and outputs
	// for inputs
	double *x, *xdot, *u, *M;
	// for dynamic model
	double *f, *df_dx, *df_dxdot, *df_du, *df_dM;
	mwIndex *df_dx_irs, *df_dx_jcs;
	mwIndex *df_dxdot_irs, *df_dxdot_jcs;
	mwIndex *df_du_irs, *df_du_jcs;
	mwIndex *df_dM_irs, *df_dM_jcs;
	mwIndex Ndf_dx, Ndf_dxdot, Ndf_du, Ndf_dM;
	// for muscle forces and moments
	double *mforces, *moments;
	// for grf model
	double *grf, *dgrf_dx, *stick;
	mwIndex Ndgrf_dx;
	mwIndex *dgrf_dx_irs, *dgrf_dx_jcs;
	// for accelerometer model
	double *fa, *dfa_dq, *dfa_dqd, *dfa_dqdd;
	mwIndex Ndfa_dq;
	mwIndex *dfa_dq_irs, *dfa_dq_jcs;
	mwIndex Ndfa_dqd;
	mwIndex *dfa_dqd_irs, *dfa_dqd_jcs;
	mwIndex Ndfa_dqdd;
	mwIndex *dfa_dqdd_irs, *dfa_dqdd_jcs;
	// the number of nonzeros of the sparse accelerometer Jacobians is already known:
	#define NNZ_DA_DQ 	38
	#define NNZ_DA_DQD 	58
	#define NNZ_DA_DQDD 86
	
	// command flags
	int CMDinitialize, CMDdynamics, CMDjointmoments, CMDstick, CMDmuscleforces;
	int CMDgrf, CMDget, CMDset, CMDmuscleCEpower, CMDaccelerometer, CMDscale;
		
	// get the first argument, see if it is a command string
	if (nrhs < 1)
		mexErrMsgTxt("At least one input argument is needed.");
    if (!mxIsChar(prhs[0])) 
		mexErrMsgTxt("First argument must be a string.");
 	if (mxGetM(prhs[0])!=1)
		mexErrMsgTxt("First argument must be a row vector.");
	command = mxArrayToString(prhs[0]);
	if(command == NULL) 
	  mexErrMsgTxt("First argument was not a command string.");
	  
	// decode the command string
	CMDinitialize = 0;
	CMDdynamics = 0;
	CMDjointmoments = 0;
	CMDstick = 0;
	CMDmuscleforces = 0;
	CMDmuscleCEpower = 0;
	CMDgrf = 0;
	CMDget = 0;
	CMDset = 0;
	CMDaccelerometer = 0;
	CMDscale = 0;

	if 		(strcmp(command, "Dynamics") == 0) 		CMDdynamics = 1;		// the most common use first
	else if (strcmp(command, "Initialize") == 0) 	CMDinitialize = 1;
	else if	(strcmp(command, "Jointmoments") == 0) 	CMDjointmoments = 1;
	else if (strcmp(command, "Stick") == 0) 		CMDstick = 1;
	else if	(strcmp(command, "Muscleforces") == 0) 	CMDmuscleforces = 1;
	else if	(strcmp(command, "MuscleCEpower") == 0) CMDmuscleCEpower = 1;
	else if (strcmp(command, "GRF") == 0) 			CMDgrf = 1;
	else if (strcmp(command, "Get") == 0) 			CMDget = 1;
	else if (strcmp(command, "Set") == 0) 			CMDset = 1;
	else if (strcmp(command, "Accelerometer") == 0) CMDaccelerometer = 1;
	else if (strcmp(command, "Scale") == 0) 		CMDscale = 1;
	else mexErrMsgTxt("gait2d: Command not recognized.");
	
	if (CMDinitialize) {
		if (nrhs != 2)
			mexErrMsgTxt("gait2d:Initialize: Two arguments required.");
		if (nlhs != 0)
			mexErrMsgTxt("gait2d:Initialize: No output argument allowed.");
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
			mexErrMsgTxt("gait2d: Incorrect type for par, must be double.");
		par = mxGetPr(prhs[1]);
		printf("***************************************************\n");
		printf("*               MEX Function gait2d               *\n");  
		printf("*       (c) 2010-2011 Orchard Kinetics LLC        *\n");
		printf("***************************************************\n");
		printf("Initializing model...\n");
		ExtractParameters(par);
		initialized = 1959;
		dynamics_evaluations = 0;
		for (i=0; i<NSTATES; i++) zeros[i] = 0.0;
		return;
	}
	
	// Not initializing, so see if Initialize has been done already
	if (initialized != 1959)
		mexErrMsgTxt("gait2d: model was not initialized.");
		
	// Is it a "Scale" command?
	if (CMDscale) {
		if (nrhs != 3) {
			mexErrMsgTxt("gait2d: Scale: exactly three input arguments needed.");				
		}
		
		// make sure we have initialized with the original parameters from Winter, and model was not scaled yet:
		if ((fabs(param.TrunkMass - 50.85)>1e-6) || (fabs(param.TrunkCMy - 0.315504)>1e-6)) {
			printf("trunk mass: %f\n",param.TrunkMass);
			printf("trunk CM: %f\n",param.TrunkCMy);
			mexErrMsgTxt("gait2d: Scale: model was already scaled.");				
		}
		
		// check the input arguments
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect type for Body Mass, must be double.");
		}
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if ((nrows != 1) || (ncols != 1) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect size for Body Mass, must be a scalar.");
		}
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect type for Body Mass, must be double.");
		}
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if ((nrows != 1) || (ncols != 1) ) {
			mexErrMsgTxt("gait2d: Scale: Incorrect size for Body Mass, must be a scalar.");
		}
		
		// get the body mass and height and calculate scaling factors
		massfactor = *mxGetPr(prhs[1]) / 75.0;		// divide subject mass by Winter standard body mass
		lengthfactor = *mxGetPr(prhs[2]) / 1.8;		// divide subject height by Winter standard body height
		inertiafactor = massfactor * lengthfactor * lengthfactor;

		// apply the scaling factors
		param.TrunkMass *= massfactor;
		param.ThighMass *= massfactor;
		param.ShankMass *= massfactor;
		param.FootMass *= massfactor;

		param.TrunkInertia *= inertiafactor;
		param.ThighInertia *= inertiafactor;
		param.ShankInertia *= inertiafactor;
		param.FootInertia *= inertiafactor;

		param.ThighLen *= lengthfactor;
		param.ShankLen *= lengthfactor;

		param.TrunkCMy *= lengthfactor;
		param.ThighCMy *= lengthfactor;
		param.ShankCMy *= lengthfactor;
		param.FootCMx *= lengthfactor;
		param.FootCMy *= lengthfactor;

		param.ContactY0 *= lengthfactor;
		param.ContactHeelX *= lengthfactor;
		param.ContactToeX  *= lengthfactor;
		
		return;
	}
		
	// Is it a "Get" command?
	if (CMDget) {
		if (nrhs != 2) {
			mexErrMsgTxt("gait2d: Get: exactly two input arguments needed.");	
		}
		if (!mxIsChar(prhs[1])) 
			mexErrMsgTxt("gait2d: Get: second argument must be a string.");
		if (mxGetM(prhs[1])!=1)
		mexErrMsgTxt("gait2d: Get second argument must be a string.");
		parametername= mxArrayToString(prhs[1]);
		if (parametername == NULL) 
			mexErrMsgTxt("gait2d: Get: second argument was not a command string.");
			
		// Decode the parameter name
		if (strcmp(parametername, "Lceopt") == 0) {
			plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Lceopt;
		}
		else if (strcmp(parametername, "Fmax") == 0) {
			plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			for (i=0; i<NMUS; i++) f[i] = muscles[i].Fmax;
		}
		else if (strcmp(parametername, "Evaluations") == 0) {
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			*f = dynamics_evaluations;
		}
		else if (strcmp(parametername, "Total Mass") == 0) {
			plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			f = mxGetPr(plhs[0]);
			*f = param.TrunkMass + 2*(param.ThighMass + param.ShankMass + param.FootMass);
		}
		else {
			printf("Parameter name: %s\n", parametername);
			mexErrMsgTxt("gait2d: getparameter: parameter name not recognized.");		
		}
		return;
	}

	// Is it a "Set" command?
	if (CMDset) {
		if (nrhs != 3 && nrhs != 4) {
			mexErrMsgTxt("gait2d: Set: three or four input arguments needed.");	
		}
		if (!mxIsChar(prhs[1])) 
			mexErrMsgTxt("gait2d: Set: second argument must be a string.");
		if (mxGetM(prhs[1])!=1)
		mexErrMsgTxt("gait2d: Set: second argument must be a string.");
		parametername= mxArrayToString(prhs[1]);
		if (parametername == NULL) 
			mexErrMsgTxt("gait2d: Set: second argument was not a command string.");
			
		// Decode the parameter name
		if (strcmp(parametername, "Extra Mass") == 0) {
			if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgTxt("gait2d: Incorrect type for Extra Mass, must be double.");
			nrows = mxGetM(prhs[2]);
			ncols = mxGetN(prhs[2]);
			if ((nrows != 1) || (ncols != 4) )
				mexErrMsgTxt("gait2d: Incorrect size for Extra Mass, must be a 1x4 matrix.");
			f = mxGetPr(prhs[2]);
			if (f[0] != 3 && f[2] != 0.0) {
				mexErrMsgTxt("gait2d: On Trunk, Thigh, Shank, extra mass must be at X=0.");
			}
			if (f[0] == 0) {
				printf("Adding %8.3f kg to Trunk at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.TrunkCMy = (param.TrunkMass * param.TrunkCMy + f[1]*f[3])/(param.TrunkMass + f[1]);
				// new mass is sum
				param.TrunkMass = param.TrunkMass + f[1];
				// new inertia according to parallel axes theorem
				param.TrunkInertia = param.TrunkInertia + f[1]*(f[3]-param.TrunkCMy)*(f[3]-param.TrunkCMy);
			}
			else if (f[0] == 1) {
				printf("Adding %8.3f kg to Thigh segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.ThighCMy = (param.ThighMass * param.ThighCMy + f[1]*f[3])/(param.ThighMass + f[1]);
				// new mass is sum
				param.ThighMass = param.ThighMass + f[1];
				// new inertia according to parallel axes theorem
				param.ThighInertia = param.ThighInertia + f[1]*(f[3]-param.ThighCMy)*(f[3]-param.ThighCMy);
			}
			else if (f[0] == 2) {
				printf("Adding %8.3f kg to Shank segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.ShankCMy = (param.ShankMass * param.ShankCMy + f[1]*f[3])/(param.ShankMass + f[1]);
				// new mass is sum
				param.ShankMass = param.ShankMass + f[1];
				// new inertia according to parallel axes theorem
				param.ShankInertia = param.ShankInertia + f[1]*(f[3]-param.ShankCMy)*(f[3]-param.ShankCMy);
			}
			else if (f[0] == 3) {
				printf("Adding %8.3f kg to Foot segments at %8.3f %8.3f\n", f[1], f[2], f[3]);
				// new CM is weighted average
				param.FootCMx = (param.FootMass * param.FootCMx + f[1]*f[2])/(param.FootMass + f[1]);
				param.FootCMy = (param.FootMass * param.FootCMy + f[1]*f[3])/(param.FootMass + f[1]);
				// new mass is sum
				param.FootMass = param.FootMass + f[1];
				// new inertia according to parallel axes theorem
				param.FootInertia = param.FootInertia 
					+ f[1]*(f[2]-param.FootCMx)*(f[2]-param.FootCMx)
					+ f[1]*(f[3]-param.FootCMy)*(f[3]-param.FootCMy);
			}
			else {
				mexErrMsgTxt("gait2d: Incorrect segment number in Set Extra Mass.  Must be 1 or 2 or 3.");
			}
				
		}
		else if (strcmp(parametername, "Right BKA") == 0) {
			// Get the ankle stiffness, specified by user as 3rd argument of the MEX function
			if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgTxt("gait2d: Incorrect type for ankle stiffness, must be double.");
			nrows = mxGetM(prhs[2]);
			ncols = mxGetN(prhs[2]);
			if ((nrows != 1) || (ncols != 1) )
				mexErrMsgTxt("gait2d: Incorrect size for ankle stiffness, must be a scalar.");
            
            // Get the ankle damping, specified by user as 4th argument of the MEX function
			if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
				mexErrMsgTxt("gait2d: Incorrect type for ankle damping, must be double.");
			nrows = mxGetM(prhs[3]);
			ncols = mxGetN(prhs[3]);
			if ((nrows != 1) || (ncols != 1) )
				mexErrMsgTxt("gait2d: Incorrect size for ankle damping, must be a scalar.");

			// Remove the ankle muscles in the right leg
			for (i=5; i<8; i++) muscles[i].Fmax = 0.0;

			// Passive elastic ankle on right leg
			param.MinAngle[2] = 0.0;
			param.MaxAngle[2] = 0.0;
			param.JointK[2] = *mxGetPr(prhs[2]);
            param.JointB[2] = *mxGetPr(prhs[3]);
		}
		else if (strcmp(parametername, "Right AKA") == 0) {
			// Get the ankle stiffness, specified by user as 3rd argument of the MEX function
			if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) )
				mexErrMsgTxt("gait2d: Incorrect type for ankle stiffness, must be double.");
			nrows = mxGetM(prhs[2]);
			ncols = mxGetN(prhs[2]);
			if ((nrows != 1) || (ncols != 1) )
				mexErrMsgTxt("gait2d: Incorrect size for ankle stiffness, must be a scalar.");
			f = mxGetPr(prhs[2]);

			// Remove the ankle muscles in the right leg
			for (i=2; i<8; i++) muscles[i].Fmax = 0.0;

			// Passive elastic ankle on right leg
			param.MinAngle[2] = 0.0;
			param.MaxAngle[2] = 0.0;
			param.JointK[2] = *mxGetPr(prhs[2]);		
		}
		else {
			printf("Parameter name: %s\n", parametername);
			mexErrMsgTxt("gait2d: Set: parameter name not recognized.");		
		}
		return;
	}
	
	// If command is "Accelerometer", we need q,qdot,qdotdot then we can finish up
	if (CMDaccelerometer) {
	
		// get the inputs from Matlab
		if (nrhs < 4) mexErrMsgTxt("gait2d: Accelerometer command needs q,qd,qdd.");
		nrows = mxGetM(prhs[1]);
		ncols = mxGetN(prhs[1]);
		if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) mexErrMsgTxt("gait2d: Incorrect type for q, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for q, must be a 9 x 1 column vector.");
		q = mxGetPr(prhs[1]);
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) mexErrMsgTxt("gait2d: Incorrect type for qd, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for qd, must be a 9 x 1 column vector.");
		qd = mxGetPr(prhs[2]);
		nrows = mxGetM(prhs[3]);
		ncols = mxGetN(prhs[3]);
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) mexErrMsgTxt("gait2d: Incorrect type for qd, must be double.");
		if ((nrows != NDOF) || (ncols != 1)) mexErrMsgTxt("gait2d: Incorrect size for qd, must be a 9 x 1 column vector.");
		qdd = mxGetPr(prhs[3]);
		
		// call the Autolev C code
		gait2d_acc(&param, q, qd, qdd, acc, dacc_dq, dacc_dqd, dacc_dqdd);

		// Create the return argument fa
		if (nlhs < 1) mexErrMsgTxt("gait2d: Accelerometer command needs at least one output.");
		plhs[0] = mxCreateDoubleMatrix(42, 1, mxREAL);
		fa = mxGetPr(plhs[0]);
		for (i=0; i<42; i++) {
			fa[i] = acc[i];
		}
		
		// Create the sparse Jacobian dfa_dq
		if (nlhs < 2) return;
		plhs[1] = mxCreateSparse(42, NDOF, NNZ_DA_DQ, 0);
		dfa_dq = mxGetPr(plhs[1]);
		dfa_dq_irs = mxGetIr(plhs[1]);
		dfa_dq_jcs = mxGetJc(plhs[1]);
		Ndfa_dq = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dq_jcs[i] = Ndfa_dq;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				dfa_dq_irs[Ndfa_dq] = j;			// store row number of this matrix element
				if (dacc_dq[j][i] != 0){				// if this matrix element is not zero
					dfa_dq[Ndfa_dq++] = dacc_dq[j][i];	// store its value
				}
			}
		}
		dfa_dq_jcs[NDOF] = Ndfa_dq;		// store final element number
		if (Ndfa_dq > NNZ_DA_DQ) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dq,NNZ_DA_DQ);
			mexErrMsgTxt("gait2d: dfa_dq has incorrect number of nonzeros.");
		}
		
		// Create the sparse Jacobian dfa_dqd
		if (nlhs < 3) return;
		plhs[2] = mxCreateSparse(42, NDOF, NNZ_DA_DQD, 0);
		dfa_dqd = mxGetPr(plhs[2]);
		dfa_dqd_irs = mxGetIr(plhs[2]);
		dfa_dqd_jcs = mxGetJc(plhs[2]);
		Ndfa_dqd = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dqd_jcs[i] = Ndfa_dqd;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				dfa_dqd_irs[Ndfa_dqd] = j;			// store row number of this matrix element
				if (dacc_dqd[j][i] != 0){					// if this matrix element is not zero
					dfa_dqd[Ndfa_dqd++] = dacc_dqd[j][i];	// store its value
				}
			}
		}
		dfa_dqd_jcs[NDOF] = Ndfa_dqd;		// store final element number
		if (Ndfa_dqd > NNZ_DA_DQD) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dqd,NNZ_DA_DQD);
			mexErrMsgTxt("gait2d: dfa_dqd has incorrect number of nonzeros.");
		}
		
		// Create the sparse Jacobian dfa_dqdd
		if (nlhs < 4) return;
		plhs[3] = mxCreateSparse(42, NDOF, NNZ_DA_DQDD, 0);
		dfa_dqdd = mxGetPr(plhs[3]);
		dfa_dqdd_irs = mxGetIr(plhs[3]);
		dfa_dqdd_jcs = mxGetJc(plhs[3]);
		Ndfa_dqdd = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of fa w.r.t. q, columns 1..NDOF
			dfa_dqdd_jcs[i] = Ndfa_dqdd;		// store element number where this column starts			
			for (j=0; j<42; j++) {			// go through all elements of this row
				
				dfa_dqdd_irs[Ndfa_dqdd] = j;			// store row number of this matrix element
				if (dacc_dqdd[j][i] != 0){						// if this matrix element is not zero
					dfa_dqdd[Ndfa_dqdd++] = dacc_dqdd[j][i];	// store its value
				}
			}
		}
		dfa_dqdd_jcs[NDOF] = Ndfa_dqdd;		// store final element number
		if (Ndfa_dqdd > NNZ_DA_DQDD) {
			printf("Number of nonzeros was %d and should be %d\n",Ndfa_dqdd,NNZ_DA_DQDD);
			mexErrMsgTxt("gait2d: dfa_dqdd has incorrect number of nonzeros.");
		}
		
		return;

	}

	// State of model is the second argument, for all other uses of this MEX function
	if (nrhs < 2)
		mexErrMsgTxt("gait2d: state vector x is needed as second argument.");
	nrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) )
		mexErrMsgTxt("gait2d: Incorrect type for state vector x, must be double.");
	if ((nrows != NSTATES) || (ncols != 1) )
		mexErrMsgTxt("gait2d: Incorrect size for state vector x, must be a 50 x 1 column vector.");
	x = mxGetPr(prhs[1]);
	
	// Use zeros for xdot and u, if they are not important
	xdot = zeros;		// xdot now points to zeros
	u = zeros;
	
	// If command is "Dynamics" or "MuscleCEpower", we also need xdot
	if (CMDdynamics || CMDmuscleCEpower || CMDaccelerometer) {
		if (nrhs < 3) mexErrMsgTxt("gait2d: xdot (argument 3) required.");
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) mexErrMsgTxt("gait2d: Incorrect type for xdot, must be double.");
		if ((nrows != NSTATES) || (ncols != 1) ) mexErrMsgTxt("gait2d: Incorrect size for xdot, must be a 50 x 1 column vector.");
		xdot = mxGetPr(prhs[2]);
	}

	if (CMDdynamics) {
		if (nrhs < 4 || nrhs > 5) {
			mexErrMsgTxt("gait2d(Dynamics): 4 or 5 inputs required.");
		}
		nrows = mxGetM(prhs[3]);
		ncols = mxGetN(prhs[3]);
		if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) )
			mexErrMsgTxt("gait2d: Incorrect type for u, must be double.");
		if ((nrows != NMUS) || (ncols != 1) )
			mexErrMsgTxt("gait2d(Dynamics): Incorrect size for u, must be a 16 x 1 column vector.");
		u = mxGetPr(prhs[3]);
		
		if (nrhs == 5) {
			// additional joint moments were given as inputs
			nrows = mxGetM(prhs[4]);
			ncols = mxGetN(prhs[4]);
			if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[3]) )
				mexErrMsgTxt("gait2d: Incorrect type for M, must be double.");
			if ((nrows != NMOM) || (ncols != 1) )
				mexErrMsgTxt("gait2d(Dynamics): Incorrect size for M, must be a 6 x 1 column vector.");
			M = mxGetPr(prhs[4]);
		}

		// Determine whether derivatives are requested
		if (nlhs == 1) 
			derivatives = 0;
		else {
			if (nlhs > 5 || nlhs < 1)
				mexErrMsgTxt("gait2d(Dynamics): must have at least 1 and no more than 5 outputs.");
			if (nlhs == 5 && nrhs != 5)
				mexErrMsgTxt("gait2d(Dynamics): 5th output df_dM requires 5th input M");
			derivatives = 1;
		}
	}
	
	// Compute the muscle dynamics, and get muscle forces
	for(i=0; i<NMUS; i++) {
	
		// Calculate dynamics residual of muscle activation model
		double rate = muscles[i].c1*u[i] + muscles[i].c2;
		h[i] = xdot[2*NDOF+NMUS+i] - rate*(u[i] - x[2*NDOF+NMUS+i]);
		if (derivatives) {
			dh_da[i] = rate;
			dh_du[i] = -rate - muscles[i].c1*(u[i] - x[2*NDOF+NMUS+i]);
			dh_dadot[i] = 1.0;
		}
	
		// Calculate muscle length Lm and derivatives dLm/dq from generalized coordinates in x
		Lm[i] = MusclePath(&muscles[i], &x[3], dLm_dang[i]);
	
		// Calculate muscle force imbalance, normalized to Fmax
		g[i] = MuscleDynamics(&muscles[i],
			x[2*NDOF+NMUS+i],		// active state of muscle i
			x[2*NDOF+i],			// Lce of muscle i
			xdot[2*NDOF+i],			// Lcedot
			Lm[i],					// muscle length
			&dg_da[i], &dg_dLce[i], &dg_dLcedot[i], &dg_dLm[i],
			&force[i], &dforce_dLm[i], &dforce_dLce[i],
			&CEpower[i]);				
	}
	
	// Compute the joint moments
	
	for (i=0; i<NMOM; i++) {
		// initialize derivatives to zero
		if (derivatives) {
			for (j=0; j<NMOM; j++) dmom_dang[i][j] = 0.0;
			for (j=0; j<NMUS; j++) dmom_dLce[i][j] = 0.0;
			dmom_dangvel[i] = 0.0;							// this is a diagonal matrix so we only store diagonal
		}
		
		// start with passive joint moment
		ang = x[i+3];											// joint angle is one of the state variables
		angvel = x[NDOF+i+3];									// the corresponding angular velocity
		// calculate deformations (positive when outside ROM)
		dmax = Softpos(ang - param.MaxAngle[i], param.JointD0[i], &ddmax);
		dmin = Softpos(param.MinAngle[i] - ang, param.JointD0[i], &ddmin);

		// calculate moments due to both endpoints and add them, also add damping
		mom[i] = -param.JointK[i] * (dmax - dmin) - param.JointB[i] * angvel;

		if (derivatives) {
			dmom_dang[i][i] = -param.JointK[i] * (ddmax + ddmin);
			dmom_dangvel[i] = -param.JointB[i];
		}
		
		// add the muscle moments
		for (j=0; j<NMUS; j++) {
			mom[i] = mom[i] - dLm_dang[j][i]*force[j];		// moment arm is -dLm/dang
			if (derivatives) {
				for (k=0; k<NMOM; k++) {
					dmom_dang[i][k] = dmom_dang[i][k] - dLm_dang[j][i]*dforce_dLm[j]*dLm_dang[j][k];
				}
				dmom_dLce[i][j] = dmom_dLce[i][j] - dLm_dang[j][i]*dforce_dLce[j];
			}
		}
		
		// add the extra moments given as inputs
		if (CMDdynamics && nrhs == 5) {
			mom[i] = mom[i] + M[i];
		}
	}
	
	// Assemble the MEX function outputs for the "Jointmoments" command
	if (CMDjointmoments) {
		plhs[0] = mxCreateDoubleMatrix(NMOM, 1, mxREAL);
		moments = mxGetPr(plhs[0]);
		// fill the column vector
		for (j=0; j<NMOM; j++)
			*moments++ = mom[j];
		return;
	}
	
	// Assemble the MEX function outputs for the "Muscleforces" command
	if (CMDmuscleforces) {
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);	
		// fill the column vector with muscle forces
		for (j=0; j<NMUS; j++)
				*mforces++ = force[j];
		return;
	}
	if (CMDmuscleCEpower) {
		plhs[0] = mxCreateDoubleMatrix(NMUS, 1, mxREAL);
		mforces = mxGetPr(plhs[0]);	
		// fill the column vector with CE power
		for (j=0; j<NMUS; j++)
				*mforces++ = CEpower[j];
				// verify that contraction dynamics equation was satisfied
				if (g[i] > 1e-4) {
					mexErrMsgTxt("gait2d(MuscleCEpower): Contraction dynamics not satisfied.");			
				}
		return;
	}

	// Call the C function that was generated by Autolev
	q = &x[0];
	qd = &x[NDOF];	
	qdd = &xdot[NDOF];
#if VERBOSE
	for (i=0; i<NDOF; i++) 
		printf("dof: %2d q=%10.5f  qd=%10.5f  qdd=%10.5f\n",i,q[i],qd[i],qdd[i]);
	for (i=0; i<NMOM; i++)
		printf("mom[%2d] = %10.5f\n",i,mom[i]);
#endif
	gait2d_dyn(&param, q, qd, qdd, mom, zero, dz_dq, dz_dqd, dz_dqdd, dz_dmom, GRF, dGRF_dq, dGRF_dqd, Stick, tmp);
#if VERBOSE
	for (i=0; i<NDOF; i++) 
		printf("zero[%2d] =  %10.5f\n",i,zero[i]);
	printf("Right foot:   GRFx=%10.5f   GRFy=%10.5f\n",GRF[0],GRF[1]);
	printf("Left  foot:   GRFx=%10.5f   GRFy=%10.5f\n",GRF[2],GRF[3]);
#endif

	// Assemble the MEX function outputs for the "Dynamics" command
	if (CMDdynamics) {
	
		// Create matrix for the return argument f
		plhs[0] = mxCreateDoubleMatrix(NSTATES, 1, mxREAL);
		f = mxGetPr(plhs[0]);
			
		// the first NDOF rows of the IDE (implicit differential equation) are: qdot-dq/dt = 0
		for (i=0; i<NDOF; i++) f[i] = x[NDOF+i] - xdot[i];
		// the next NDOF rows of the IDE are the equations of motion from Autolev
		for (i=0; i<NDOF; i++) f[NDOF+i] = zero[i];
		// the next NMUS rows of the IDE are the muscle contraction dynamics
		for (i=0; i<NMUS; i++) f[2*NDOF+i] = g[i];
		// the final NMUS rows of the IDE are the muscle activation dynamics: da/dt - (u-a)(c1*u + c2) = 0
		for (i=0; i<NMUS; i++) f[2*NDOF+NMUS+i] =  h[i];
		
		if (derivatives) {
			// The sparse Jacobians have to be filled in column order, using Matlab sparse data structure
			// --------Jacobian df/dx
			plhs[1] = mxCreateSparse(NSTATES,NSTATES, 
				NDOF*(1 + 2*NDOF + NMUS) + 3*NMUS + nonzeromomentarms, 0);
			df_dx = mxGetPr(plhs[1]);
			df_dx_irs = mxGetIr(plhs[1]);
			df_dx_jcs = mxGetJc(plhs[1]);
			Ndf_dx = 0;

			// derivatives with respect to q, columns 1..NDOF of x
			for (i=0; i<NDOF; i++) {			
				df_dx_jcs[i] = Ndf_dx;				// store element number where this column starts
				
				// derivatives of ZERO with respect to q are in rows NDOF+1 to 2*NDOF
				for (j=0; j<NDOF; j++) {
					df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
					df_dx[Ndf_dx] = dz_dq[j][i];	// store the value of this matrix element
					// add the contributions dz/dmom * dmom/dq, these are nonzero when q is an angle
					if (i>2) for (k=0; k<NMOM; k++) {
						df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][k]*dmom_dang[k][i-3];
					}
					Ndf_dx++;
				}
				
				// derivatives of muscle imbalance with respect to q are in rows 2*NDOF+1 to 2*NDOF+NMUS
				for (j=0; j<NMUS; j++) {
					// element only exists if i is an angle
					if (i>2) if (dLm_dang[j][i-3] != 0) {
						df_dx_irs[Ndf_dx] = 2*NDOF+j;		// store row number of this matrix element
						df_dx[Ndf_dx] = dg_dLm[j]*dLm_dang[j][i-3];	// store the value of this matrix element
						Ndf_dx++;
					}
				}
			}
			
			// derivatives with respect to qdot, columns NDOF+1 to 2*NDOF of x
			for (i=0; i<NDOF; i++) {
				df_dx_jcs[NDOF+i] = Ndf_dx;		// store element number where this column starts
				
				// derivatives of (qdot-dq/dt) with respect to qdot are diagonal, in rows 1 to NDOF
				df_dx_irs[Ndf_dx] = i;		// store row number of this matrix element
				df_dx[Ndf_dx] = 1.0;		// store the value of this matrix element
				Ndf_dx++;
				
				// derivatives of ZERO with respect to qdot are in rows NDOF+1 to 2*NDOF
				for (j=0; j<NDOF; j++) {
					// vertical acceleration zero[1] is not dependent on horizontal velocity qdot[0]
					// "if" is needed because Autolev does not consistently produce hard zeros for those
					if (i != 0 || j != 1) {		
						df_dx_irs[Ndf_dx] = NDOF+j;			// store row number of this matrix element
						df_dx[Ndf_dx] = dz_dqd[j][i];		// store the value of this matrix element
						// add the contributions dz/dmom * dmom/dqdot, these are nonzero when qdot is an ang. vel
						if (i >= 3) {
							df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][i-3]*dmom_dangvel[i-3];
						}
						Ndf_dx++;
					}
				}		
			}		
			
			// derivatives with respect to Lce, columns 2*NDOF+1 to 2*NDOF+NMUS
			for (i=0; i<NMUS; i++) {			
				df_dx_jcs[2*NDOF+i] = Ndf_dx;		// store element number where this column starts
			
				// derivatives of ZERO with respect to Lce are in rows NDOF+1 to 2*NDOF
				for (j=0; j<NDOF; j++) {
					df_dx_irs[Ndf_dx] = NDOF+j;		// store row number of this matrix element
					df_dx[Ndf_dx] = 0.0;			// store the value of this matrix element
					for (k=0; k<NMOM; k++) {
						// next line is a bit tricky, neg sign for moment arm, neg sign for dF/dLce cancel
						df_dx[Ndf_dx] = df_dx[Ndf_dx] + dz_dmom[j][k]*dmom_dLce[k][i];
					}
					Ndf_dx++;
				}
				
				// derivatives of muscle force balance with respect to Lce are diagonal, rows 2*NDOF to 2*NDOF+NMUS
				df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
				df_dx[Ndf_dx] = dg_dLce[i];		// store the value of this matrix element
				Ndf_dx++;				
			}
			
			// derivatives with respect to Act, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
			for (i=0; i<NMUS; i++) {			
				df_dx_jcs[2*NDOF+NMUS+i] = Ndf_dx;		// store element number where this column starts
			
				// derivatives of muscle force balance with respect to Act are diagonal, rows 2*NDOF+1 to 2*NDOF+NMUS
				df_dx_irs[Ndf_dx] = 2*NDOF+i;		// store row number of this matrix element
				df_dx[Ndf_dx] = dg_da[i];		// store the value of this matrix element
				Ndf_dx++;				
			
				// derivatives of activation dynamics with respect to Act are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
				df_dx_irs[Ndf_dx] = 2*NDOF+NMUS+i;		// store row number of this matrix element
				df_dx[Ndf_dx] = dh_da[i];		// store the value of this matrix element
				Ndf_dx++;				
			}
			df_dx_jcs[NSTATES] = Ndf_dx;		// store final element number

			// --------Jacobian df/dxdot
			if (nlhs > 2) {
				plhs[2] = mxCreateSparse(NSTATES, NSTATES, NDOF*(NDOF+1) + 2*NMUS, 0);
				df_dxdot = mxGetPr(plhs[2]);
				df_dxdot_irs = mxGetIr(plhs[2]);
				df_dxdot_jcs = mxGetJc(plhs[2]);
				Ndf_dxdot = 0;

				// derivatives with respect to dq/dt, columns 1..NDOF of xdot
				for (i=0; i<NDOF; i++) {			
					df_dxdot_jcs[i] = Ndf_dxdot;			// store element number where this column starts

					// derivatives of (qdot-dq/dt) with respect to dq/dt are diagonal, in rows 1 to NDOF
					df_dxdot_irs[Ndf_dxdot] = i;		// store row number of this matrix element
					df_dxdot[Ndf_dxdot] = -1.0;			// store the value of this matrix element
					Ndf_dxdot++;
				}
				
				// derivatives with respect to dqdot/dt, columns NDOF+1 to 2*NDOF of xdot
				for (i=0; i<NDOF; i++) {
					df_dxdot_jcs[NDOF+i] = Ndf_dxdot;		// store element number where this column starts
					
					// derivatives of ZERO with respect to qdd are in rows NDOF to 2*NDOF
					for (j=0; j<NDOF; j++) {
						df_dxdot_irs[Ndf_dxdot] = NDOF+j;		// store row number of this matrix element
						df_dxdot[Ndf_dxdot] = dz_dqdd[j][i];	// store the value of this matrix element
						Ndf_dxdot++;
					}
				}
				
				// derivatives with respect to Lcedot, columns 2*NDOF+1 to 2*NDOF+NMUS
				for (i=0; i<NMUS; i++) {			
					df_dxdot_jcs[2*NDOF+i] = Ndf_dxdot;		// store element number where this column starts
								
					// derivatives of muscle force balance with respect to Lcedot are diagonal, rows 2*NDOF to 2*NDOF+NMUS
					df_dxdot_irs[Ndf_dxdot] = 2*NDOF+i;		// store row number of this matrix element
					df_dxdot[Ndf_dxdot] = dg_dLcedot[i];		// store the value of this matrix element
					Ndf_dxdot++;
				}
				
				// derivatives with respect to Actdot, columns 2*NDOF+NMUS+1 to 2*NDOF+2*NMUS
				for (i=0; i<NMUS; i++) {			
					df_dxdot_jcs[2*NDOF+NMUS+i] = Ndf_dxdot;		// store element number where this column starts
					
					// derivatives of activation dynamics with respect to Actdot are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
					df_dxdot_irs[Ndf_dxdot] = 2*NDOF+NMUS+i;		// store row number of this matrix element
					df_dxdot[Ndf_dxdot] = dh_dadot[i];					// store the value of this matrix element
					Ndf_dxdot++;				
				}
				df_dxdot_jcs[NSTATES] = Ndf_dxdot;		// store final element number
			}

			// --------Jacobian df/du
			if (nlhs > 3) {
				plhs[3] = mxCreateSparse(NSTATES, NMUS, NMUS, 0);
				df_du = mxGetPr(plhs[3]);
				df_du_irs = mxGetIr(plhs[3]);
				df_du_jcs = mxGetJc(plhs[3]);
				Ndf_du = 0;
				
				// derivatives with respect to u of each muscle
				for (i=0; i<NMUS; i++) {			
					df_du_jcs[i] = Ndf_du;		// store element number where this column starts
				
					// derivatives of activation dynamics with respect to u are diagonal, rows 2*NDOF+NMUS to 2*NDOF+2*NMUS
					df_du_irs[Ndf_du] = 2*NDOF+NMUS+i;		// store row number of this matrix element
					df_du[Ndf_du] = dh_du[i];				// store the value of this matrix element
					Ndf_du++;				
				}
				df_du_jcs[NMUS] = Ndf_du;		// store final element number
			}
			
			//-----------Jacobian df/dM (if requested) we already have made sure that M inputs were given
			if (nlhs > 4) {
				plhs[4] = mxCreateSparse(NSTATES, NMOM, NDOF*NMOM, 0);  // NOTE: actual number of nonzeros will be a lot smaller!
				df_dM = mxGetPr(plhs[4]);
				df_dM_irs = mxGetIr(plhs[4]);
				df_dM_jcs = mxGetJc(plhs[4]);
				Ndf_dM = 0;
				
				// derivatives with respect to M at each joint
				for (i=0; i<NMOM; i++) {			
					df_dM_jcs[i] = Ndf_dM;		// store element number where this column starts
					
					for (j=0; j<NDOF; j++) {	// we only need to go into the NDOF equation of motion rows
						// derivatives of Zero with respect to each moment, rows NDOF+1 to 2*NDOF
						if (dz_dmom[j][i] != 0) {
							df_dM_irs[Ndf_dM] = NDOF+j;			// store row number of this matrix element
							df_dM[Ndf_dM] = dz_dmom[j][i];		// store the value of this matrix element
							Ndf_dM++;			
						}
					}
				}
				df_dM_jcs[NMOM] = Ndf_dM;		// store final element number
			}
		}
		dynamics_evaluations++;
		return;
	}

	// Assemble the MEX function outputs for the "GRF" command
	if (CMDgrf) {
		plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
		grf = mxGetPr(plhs[0]);	
		// fill the columns
		for (i=0; i<4; i++) *grf++ = GRF[i];
		if (nlhs < 2) return;
		
		// also return derivatives dGRF/dx, a 4 x 50 matrix with at most 4 x 18 nonzeros
		plhs[1] = mxCreateSparse(4, NSTATES, 4*18, 0);
		dgrf_dx = mxGetPr(plhs[1]);
		dgrf_dx_irs = mxGetIr(plhs[1]);
		dgrf_dx_jcs = mxGetJc(plhs[1]);
		Ndgrf_dx = 0;
		for (i=0; i<NDOF; i++) {			// derivatives of GRF w.r.t. q, columns 1..NDOF of x		
			dgrf_dx_jcs[i] = Ndgrf_dx;		// store element number where this column starts			
			for (j=0; j<4; j++) {
				dgrf_dx_irs[Ndgrf_dx] = j;			// store row number of this matrix element
				dgrf_dx[Ndgrf_dx] = dGRF_dq[j][i];	// store the value of this matrix element
				Ndgrf_dx++;	
			}
		}
		for (i=0; i<NDOF; i++) {				// derivatives of GRF w.r.t. qdot, columns NDOF+1..2*NDOF of x
			dgrf_dx_jcs[NDOF+i] = Ndgrf_dx;		// store element number where this column starts			
			for (j=0; j<4; j++) {
				dgrf_dx_irs[Ndgrf_dx] = j;				// store row number of this matrix element
				dgrf_dx[Ndgrf_dx] = dGRF_dqd[j][i];	// store the value of this matrix element
				Ndgrf_dx++;	
			}
		}
		for (i=2*NDOF; i<NSTATES; i++) {	// derivatives of GRF w.r.t. other state vars are zero
			dgrf_dx_jcs[i] = Ndgrf_dx;		// store element number where this column starts						
		}
		dgrf_dx_jcs[NSTATES] = Ndgrf_dx;		// store final element number
		if (nlhs < 3) return;
		
		// also return the tmp array
		plhs[2] = mxCreateDoubleMatrix(16, 1, mxREAL);
		f = mxGetPr(plhs[2]);
		for (i=0; i<16; i++) f[i] = tmp[i];
		if (nlhs < 4) return;

		// now it is an error
		mexErrMsgTxt("gait2d(GRF): GRF function has 1 or 2 or 3 outputs only.");			
		return;	
	}
	
	// Assemble the MEX function outputs for the "Stick" command
	if (CMDstick) {
		plhs[0] = mxCreateDoubleMatrix(NSTICK, 2, mxREAL);
		stick = mxGetPr(plhs[0]);	
		// fill the columns
		for (i=0; i<2; i++)
			for (j=0; j<NSTICK; j++)
				*stick++ = Stick[j][i];
		return;
	}
			
    mexErrMsgTxt("gait2d: Programming error.");

}
