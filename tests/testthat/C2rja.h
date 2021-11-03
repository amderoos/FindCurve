/*
  NAME
    C2rja
  PURPOSE
    Module computes the internal equilibrium of the stage-structured biomass model
    with two resources, a possible niche shift at maturation and pulsed reproduction

  Last modification: AMdR - Jan 26, 2018
 */

/*
 *====================================================================================================================================
 *  DEFINITION OF PROBLEM DIMENSIONS AND NUMERICAL SETTINGS
 *====================================================================================================================================
 */
// Dimension settings: Required
#define EQUATIONS_DIM             4
#define EXTRAOUTPUT_DIM           6
#define PARAMETER_NR              14

// Numerical settings: Optional (default values adopted otherwise)
#define DYTOL                     1.0E-8                                            // Variable tolerance
#define RHSTOL                    1.0E-9                                            // Function tolerance
#define ALLOWNEGATIVE             0                                                 // Negative solution components allowed?


/*
 *====================================================================================================================================
 *  DEFINITION OF ALIASES
 *====================================================================================================================================
 */
// Define aliases for the parameters
#define DELTA                     parameter[0]
#define R1MAX                     parameter[1]
#define R2MAX                     parameter[2]

// Parameters for consumer
#define SIGMA                     parameter[3]

#define ZC                        parameter[4]

#define IMAX                      parameter[5]
#define Q                         parameter[6]
#define E                         parameter[7]
#define MUC                       parameter[8]

#define ETA                       parameter[9]

#define MUCPLUS                   parameter[10]
#define MUCJPLUS                  parameter[11]
#define MUCAPLUS                  parameter[12]

// Parameters for consumer related to discrete reproduction
#define INTERVAL                  parameter[13]

/*
 *====================================================================================================================================
 *  DEFINITION OF NAMES AND DEFAULT VALUES OF THE PARAMETERS
 *====================================================================================================================================
 */
// At least two parameters should be specified in this array
char *parameternames[PARAMETER_NR] = {"Delta", "Rmax1", "Rmax2", "Sigma",   "Zc",       "Imax",     "Q",
                                      "E",     "Muc",   "Eta",   "MucPlus", "MucjPlus", "MucaPlus", "Interval"};

// These are the default parameters values
double parameter[PARAMETER_NR] = {10.0, 2.0, 2.0, 0.5, 0.1, 100.0, 10.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

/*
 *====================================================================================================================================
 *  DEFINITION OF THE SYSTEM OF EQUATIONS TO SOLVE
 *====================================================================================================================================
 */

#undef MAX_EXP
#define MAX_EXP                   50.0

double Maturation(double z, double nuj, double muj)

{
  double logz, tmp, tres, matrate = 0.0;

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

#define ODEDIM                    11
static double                     StoredVals[ODEDIM];
static int                        OdeDim = ODEDIM;

#define R1                        argument[0]                                       // Resource #1
#define R2                        argument[1]                                       // Resource #2
#define CJ                        argument[2]                                       // Juvenile consumers
#define CA                        argument[3]                                       // Adult consumers
#define B                         argument[4]                                       // energy storage

#define DR1DT                     derivative[0]                                     // Resource #1
#define DR2DT                     derivative[1]                                     // Resource #2
#define DCJDT                     derivative[2]                                     // Juvenile consumers
#define DCADT                     derivative[3]                                     // Adult consumers
#define DBDT                      derivative[4]                                     // Energy storage

void WithinSeason(double t, double *argument, double *derivative)
{
  static double ingest_R1, ingest_R2;
  static double nu_J, nu_A;
  static double mort_J, mort_A, maturation;

  nu_J = SIGMA*IMAX*R1 - Q;
  nu_A = SIGMA*IMAX*((ETA*R1 + (1.0 - ETA)*R2)) - Q*(1 + E*B/CA);

  ingest_R1 = IMAX*(R1*CJ + ETA*R1*CA);
  ingest_R2 = IMAX*(1.0 - ETA)*R2*CA;

  mort_J = MUC + MUCJPLUS + MUCPLUS;                                                // Starvation mortality (max(nu_J, 0.0)-nu_J) not included here !
  mort_A = MUC + MUCAPLUS + MUCPLUS + max(nu_A, 0.0) - nu_A;

  maturation = Maturation(ZC, nu_J, mort_J)*CJ;

  DR1DT = DELTA*(R1MAX - R1) - ingest_R1;
  DR2DT = DELTA*(R2MAX - R2) - ingest_R2;

  DCJDT = nu_J*CJ - maturation - mort_J*CJ;
  DCADT = maturation - mort_A*CA;
  DBDT  = max(nu_A, 0.0)*CA - mort_A*B;

  // Integrate the following ODEs only for output purposes
  if (OdeDim == ODEDIM)
    {
      derivative[5]  = R1;
      derivative[6]  = R2;
      derivative[7]  = CJ;
      derivative[8]  = CA;
      derivative[9]  = CA + B;
      derivative[10] = CJ + CA + B;
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
  x[0] = R1;
  x[1] = R2;
  x[2] = CJ;
  x[3] = CA;

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
      x[2] += x[4];                                                                 // Add reproductive mass to juveniles
      x[4] = 0.0;                                                                   // Reset reproductive mass
    }

  //================================================================================
  // Compute the final values of the fixed point equation F(y)=0,

  if (result)
    {
      result[0] = (R1 - x[0]);
      result[1] = (R2 - x[1]);
      result[2] = (CJ - x[2]);
      result[3] = (CA - x[3]);
    }

  return SUCCES;
}


/*==================================================================================================================================*/

// Define all variables to be written to the output file (column-organized ASCII file)

int DefineExtraOutput(double *argument, double *ExtraOutput)

{
  // Invoke the routine that sets the right-hand side for setting the output variables
  if (Equations(argument, NULL) == FAILURE) return FAILURE;

  ExtraOutput[0] = StoredVals[5]/INTERVAL;                                          // Average resource 1
  ExtraOutput[1] = StoredVals[6]/INTERVAL;                                          // Average resource 2
  ExtraOutput[2] = StoredVals[7]/INTERVAL;                                          // Avergae Cj
  ExtraOutput[3] = StoredVals[8]/INTERVAL;                                          // Average Ca
  ExtraOutput[4] = StoredVals[9]/INTERVAL;                                          // Average Ca+Cb
  ExtraOutput[5] = StoredVals[10]/INTERVAL;                                         // Average Cj+Ca+Cb

  return SUCCES;
}


/*==================================================================================================================================*/
