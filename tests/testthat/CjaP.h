/***
   NAME
     CjaP
   PURPOSE
     Module computes the internal equilibrium of the stage-structured biomass model
     with a single resource and pulsed reproduction and a single predator

    Last modification: AMdR - Jan 26, 2018
 ***/

/*
 *====================================================================================================================================
 *  PROGRAM SETTINGS
 *====================================================================================================================================
 */

#define RESFLUC                   0                                                 // Set resource dynamics to Semi-chemostat ( = 0) or fluctuating dynamics ( = 1)

/*
 *====================================================================================================================================
 *  DEFINITION OF PROBLEM DIMENSIONS AND NUMERICAL SETTINGS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define EQUATIONS_DIM             4
#define EXTRAOUTPUT_DIM           6
#define PARAMETER_NR              24

// Numerical settings: Optional (default values adopted otherwise)
#define ODESOLVE_MAX_STEP         0.1                                               // Largest step size in odesolver
#define ODESOLVE_REL_ERR          1.0E-10

#define DYTOL                     1.0E-8                                            // Variable tolerance
#define RHSTOL                    1.0E-9                                            // Function tolerance
#define ALLOWNEGATIVE             0                                                 // Negative solution components allowed?


/*
 *====================================================================================================================================
 *  DEFINITION OF ALIASES
 *====================================================================================================================================
 */
// Define aliases for the parameters
#define DELTAR                    parameter[0]
#define RMAX                      parameter[1]

// Parameters for consumer general
#define MC                        parameter[2]
#define TC                        parameter[3]
#define SIGMAC                    parameter[4]
#define ZC                        parameter[5]                                      // Body size maturation ratio
#define MUC                       parameter[6]                                      // Background mortality
#define QC                        parameter[7]
#define HC                        parameter[8]

#define MUCPLUS                   parameter[9]
#define MUCJPLUS                  parameter[10]
#define MUCAPLUS                  parameter[11]

// Parameters for consumer related to discrete reproduction
#define EPS                       parameter[12]                                     // Conversion factor maintenance rate gonads
#define INTERVAL                  parameter[13]                                     // Duration of period between spawning events
#define PSI                       parameter[14]                                     // Rep. enrg. all. 2 cont. rep.

// Parameters for predator
#define MP                        parameter[15]
#define TP                        parameter[16]
#define SIGMAP                    parameter[17]
#define MUP                       parameter[18]
#define PHI                       parameter[19]                                     // Preference for Juvenile stage of consumer
#define HP                        parameter[20]

#define MUPLUSP                   parameter[21]

#define FUNCRESP                  parameter[22]

#define ARESF                     parameter[23]


/*
 *====================================================================================================================================
 *  DEFINITION OF NAMES AND DEFAULT VALUES OF THE PARAMETERS
 *====================================================================================================================================
 */
// At least two parameters should be specified in this array
char *parameternames[PARAMETER_NR] = {"Delta", "Rmax",     "Mc",  "Tc", "Sigmac", "Zc",     "Muc", "Qc",  "Hc", "MucPlus", "MucjPlus", "MucaPlus",
                                      "EPS",   "Interval", "Psi", "Mp", "Tp",     "Sigmap", "Mup", "Phi", "Hp", "Muplusp", "FuncResp", "AresF"};

// These are the default parameters values
double parameter[PARAMETER_NR] = {0.1, 30.0, 1.0, 0.1,  0.5,   0.1, 0.015, 1.5, 3.0, 0.0, 0.0, 0.0,
                                  0.0, 70.0, 0.0, 0.32, 0.032, 0.5, 0.005, 0.0, 3.0, 0.0, 1.0, 0.0};

/*
 *====================================================================================================================================
 *  DEFINITION OF THE SYSTEM OF EQUATIONS TO SOLVE
 *====================================================================================================================================
 */

#undef MAX_EXP
#define MAX_EXP                   50.0

double Maturation(double z, double nuj, double muj)

{
  double  logz, tmp, tres, matrate = 0.0;

  logz = log(z);
  tres = muj/(1.0 - MAX_EXP/logz);
  if (nuj < tres)
    matrate = 0.0;
  else
    {
      tmp = 1.0 - muj/nuj;
      if (fabs(tmp) < 1.0E-6)
        matrate = tmp/2 - 1/logz;
      else
        matrate = tmp/(1.0 - exp(tmp*logz));
    }
  matrate *= nuj;

  return matrate;
}


/*==================================================================================================================================*/

// The ODE system defining the change in state variables during the growth period

#define ODEDIM 11
static double                     StoredVals[ODEDIM];
static int                        OdeDim = ODEDIM;

#define R                         argument[0]                                       // Resource
#define CJ                        argument[1]                                       // Juvenile consumers
#define CA                        argument[2]                                       // Adult consumers
#define P                         argument[3]                                       // Predator biomass
#define B                         argument[4]                                       // energy storage

#define DRDT                      derivative[0]                                     // Resource
#define DCJDT                     derivative[1]                                     // Juvenile consumers
#define DCADT                     derivative[2]                                     // Adult consumers
#define DPDT                      derivative[3]                                     // Predator
#define DBDT                      derivative[4]                                     // Energy storage

void WithinSeason(double t, double *argument, double *derivative)
{
  static double ingest_R = 0.0;
  static double nu_J, nu_A;
  static double mort_J, mort_A, maturation, encP;
#if (RESFLUC == 1)
  static double rfluc;
#endif

  ingest_R = MC*((2 - QC)*CJ + QC*CA)*R/((1 - FUNCRESP) + FUNCRESP*(R + HC));
  nu_J     = SIGMAC*(2 - QC)*MC*R/((1 - FUNCRESP) + FUNCRESP*(R + HC)) - TC;
  nu_A     = SIGMAC*QC*MC*R/((1 - FUNCRESP) + FUNCRESP*(R + HC)) - TC*(1 + EPS*B/CA);

  encP = PHI*CJ + (1 - PHI)*(CA + B);

  // mort_J   = MUC - (nu_J - max(nu_J, 0.0));
  // Starvation mortality (max(nu_J, 0.0)-nu_J) included here !
  mort_J = MUC + MUCJPLUS + MUCPLUS + max(-nu_J, 0.0) + MP*PHI*P/((1 - FUNCRESP) + FUNCRESP*(encP + HP));
  mort_A = MUC + MUCAPLUS + MUCPLUS + max(-nu_A, 0.0) + MP*(1 - PHI)*P/((1 - FUNCRESP) + FUNCRESP*(encP + HP));

  maturation = Maturation(ZC, nu_J, mort_J)*CJ;

#if (RESFLUC == 0)
  DRDT = DELTAR*(RMAX - R) - ingest_R;
#elif (RESFLUC == 1)
  rfluc = RMAX*(1 + ARESF*sin(2*M_PI*t/INTERVAL));
  DRDT  = DELTAR*(rfluc - R) - ingest_R;
#endif

  DCJDT = max(nu_J, 0.0)*CJ - maturation - mort_J*CJ + PSI*max(nu_A, 0.0)*CA;
  DCADT = maturation - mort_A*CA;
  DPDT  = (SIGMAP*MP*encP/((1 - FUNCRESP) + FUNCRESP*(encP + HP)) - TP - MUP - MUPLUSP)*P;
  DBDT  = (1 - PSI)*max(nu_A, 0.0)*CA - mort_A*B;

  // Integrate the following ODEs only for output purposes
  if (OdeDim == ODEDIM)
    {
      derivative[5]  = R;
      derivative[6]  = CJ;
      derivative[7]  = CA + B;
      derivative[8]  = CJ + CA + B;
      derivative[9]  = P;
      derivative[10] = (SIGMAP*MP*encP/((1 - FUNCRESP) + FUNCRESP*(encP + HP)) - TP - MUP - MUPLUSP);
    }

  return;
}


/*==================================================================================================================================*/

// Routine specifying the system of equalities from which to solve for
// R, J and A at equilibrium
#define PERIOD 1

int Equations(double *argument, double *result)

{
  int     period;
  double  tval, tend, x[ODEDIM];

  //================================================================================
  // Set the initial point for the ODEs

  memset(x, 0, ODEDIM*sizeof(double));
  x[0] = R;
  x[1] = CJ;
  x[2] = CA;
  x[3] = P;

  if (result)
    OdeDim = 5;
  else
    OdeDim = ODEDIM;
  tval     = 0.0;
  tend     = INTERVAL;

  for (period = 0; period < PERIOD; period++)
    {
      tend = (period + 1)*INTERVAL;
      // Integrate up to end of the growing phase
      if (odesolve(x, OdeDim, &tval, tend, WithinSeason, NULL) == FAILURE)
        {
          ErrorMsg(__FILE__, __LINE__, "Integration failed!");
          return FAILURE;
        }
      if (!result) memcpy(StoredVals, x, ODEDIM*sizeof(double));
      x[1] += x[4];                                                                 // Add reproductive mass to juveniles
      x[4] = 0.0;                                                                   // Reset reproductive mass
    }

  //================================================================================
  // Compute the final values of the fixed point equation F(y)=0,

  if (result)
    {
      result[0] = (R - x[0]);
      result[1] = (CJ - x[1]);
      result[2] = (CA - x[2]);
      result[3] = (P - x[3]);
    }

  return SUCCES;
}


/*==================================================================================================================================*/

// Define all variables to be written to the output file (column-organized ASCII file)

int DefineExtraOutput(double *argument, double *ExtraOutput)

{
  // Invoke the routine that sets the right-hand side for setting the output variables
  if (Equations(argument, NULL) == FAILURE) return FAILURE;

  ExtraOutput[0] = StoredVals[5]/INTERVAL;                                          // Average resource
  ExtraOutput[1] = StoredVals[6]/INTERVAL;                                          // Average Cj
  ExtraOutput[2] = StoredVals[7]/INTERVAL;                                          // Avergae Ca+Cb
  ExtraOutput[3] = StoredVals[8]/INTERVAL;                                          // Average Cj+Ca+Cb
  ExtraOutput[4] = StoredVals[9]/INTERVAL;                                          // Average P
  ExtraOutput[5] = StoredVals[10]/INTERVAL;                                         // P.C. growth rate P

  return SUCCES;
}


/*==================================================================================================================================*/
