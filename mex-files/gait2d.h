// walk2d.h
// This file defines the data structure that holds model parameters.
// It is needed to pass model parameters to the Autolev generated C code
// in walk2d_al.c

#define NDOF 9		/* number of kinematic degrees of freedom */
#define NMOM 6		/* number of joint moments */
#define NSTICK 10	/* number of stick figure points */

typedef struct {

// Body segment parameters
	double TrunkMass, TrunkInertia, TrunkCMy;
	double ThighMass, ThighInertia, ThighCMy, ThighLen;
	double ShankMass, ShankInertia, ShankCMy, ShankLen;
	double FootMass, FootInertia, FootCMx, FootCMy;

// Parameters for the ground contact model
	double ContactY;
	double ContactHeelX, ContactToeX;
	double ContactStiff, ContactDamp, ContactY0, ContactV0, ContactFric;

// Other parameters	
	double MinAngle[NMOM], MaxAngle[NMOM];				// passive joint range of motion
	double JointK[NMOM], JointD0[NMOM], JointB[NMOM];	// passive joint stiffness and damping parameters
	double HillA;							// Normalized Hill parameter a/Fmax for F-v relationship (usually 0.25)
	double Gmax;							// Maximum eccentric muscle force, relative to Fmax
	
} param_struct;

// prototype for the Autolev C function for multibody dynamics
void gait2d_dyn(param_struct* par, 
	double q[NDOF], 
	double qd[NDOF], 
	double qdd[NDOF],
	double mom[NMOM],
	double Zero[NDOF],
	double dz_dq[NDOF][NDOF], 
	double dz_dqd[NDOF][NDOF],
	double dz_dqdd[NDOF][NDOF],
	double dz_dmom[NDOF][NMOM],
	double GRF[4],
	double dGRF_dq[4][NDOF],
	double dGRF_dqd[4][NDOF],
	double Stick[NSTICK][2], double tmp[16]);
// prototype for the Autolev C function for accelerometer model
void gait2d_acc(param_struct* par,
	double q[NDOF],
	double qd[NDOF],
	double qdd[NDOF],
	double acc[42],
	double dacc_dq[42][NDOF],
	double dacc_dqd[42][NDOF],
	double dacc_dqdd[42][NDOF]);
