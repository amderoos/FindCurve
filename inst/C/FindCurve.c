/*
  NAME
    FindCurve

  PURPOSE
    Generic, problem-independent specification for curve continuation problems
    specified by a system of non-linear equations determining a fixed point of
    an arbitrary number of unknowns to be determined. All problem-specific
    fixed-point conditions are specified in an include file

    Copyright (C) 2015, Andre M. de Roos, University of Amsterdam

    This file is part of the FindCurve software package.

    FindCurve is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    FindCurve is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FindCurve. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Oct 15, 2023
*/

#include "globals.h"

/*
 *====================================================================================================================================
 *  Import the population and environment dimension settings
 *====================================================================================================================================
 */
static int                        DefaultBifparone;
static int                        DefaultBifpartwo;

// Global flags to tailor execution
static int                        TestRun  = 0, ReportLevel = 1, pntdim;
static int                        BPdetection = 1;
static int                        EXTdetection = 1;
static int                        LPdetection = 1;
static int                        AllowNegative = 0;
static int                        DoSingle = 0;

static char                       ContinuationString[MAX_STR_LEN];
static char                       runname[MAXPATHLEN];
static char                       progname[MAXPATHLEN];

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
static char                       pntstring[MAX_STR_LEN];
static char                       curvestring[MAX_STR_LEN];
static char                       parstring[MAX_STR_LEN];
static char                       optstring[MAX_STR_LEN];
#endif

#define FINDCURVE                 1

#if defined(PROBLEMHEADER)
#define HEADERNAME <PROBLEMHEADER>
#include HEADERNAME                                                                 // Include header file
#else
#error No header file defined!
#endif

#include "defaults.h"

#if !defined(EQUATIONS_DIM) || (EQUATIONS_DIM < 1)
#error EQUATIONS_DIM should be larger than 0
#endif

#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
#undef EXTRAOUTPUT_DIM
#define EXTRAOUTPUT_DIM           0
#endif

#if !defined(PARAMETER_NR) || (PARAMETER_NR < 3)
#error PARAMETER_NR should be defined larger than 2
#endif

#undef MAX_PNTDIM
#define MAX_PNTDIM                (2*EQUATIONS_DIM + 3)                             // N variables, 2 parameters, 1 eigenvalue & N eigenvector elements

/*
 *====================================================================================================================================
 *  Definition of global variables and parameters
 *====================================================================================================================================
 */

// These are the variables to solve for
static double                     Sol[EQUATIONS_DIM];


// Global variables for other purposes

static double                     point[MAX_PNTDIM];
static double                     point_scales[MAX_PNTDIM];
static double                     pntmin[MAX_PNTDIM];
static double                     pntmax[MAX_PNTDIM];

/*
 *====================================================================================================================================
 *  Implementation of problem specification routines
 *====================================================================================================================================
 */

int EQsystem(double *argument, double *result)

{
  int index = 0, i;

  //******************************************************************************
  // Map current estimate of solution to global variables
  parameter[Bifparone] = argument[index]*pnt_scale[index];
  index++;

  for (i = 0; i < EQUATIONS_DIM; i++)
    {
      Sol[i] = argument[index]*pnt_scale[index];
      index++;
    }

  if (CurveType != EQ) parameter[Bifpartwo] = argument[index]*pnt_scale[index];
  index++;

  // Compute the basic system of equations
  Equations(Sol, result);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  if (checkInterrupt()) return FAILURE;
#endif

  return SUCCES;
}


/*==================================================================================================================================*/

int DerivEquations(double *argument, double *result)

{
  int     index = 0, i;
  double  storedval[EQUATIONS_DIM];
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  double  extra[EXTRAOUTPUT_DIM], storedpar;
#endif

  for (i = 0; i < EQUATIONS_DIM; i++)
    {
      storedval[i] = Sol[i];
      Sol[i]       = argument[index];
      index++;
    }

#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  storedpar         = parameter[EXTpar];
  parameter[EXTpar] = argument[index]*storedpar;
#endif

  // Compute the basic system of equations
  Equations(Sol, result);

#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  // Also add the additional output variables
  DefineExtraOutput(Sol, extra);

  result[EQUATIONS_DIM] = extra[EXTfun];
#endif

  // Reset the global variables
  for (i = 0; i < EQUATIONS_DIM; i++) Sol[i] = storedval[i];

#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  parameter[EXTpar]                          = storedpar;
#endif

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  if (checkInterrupt()) return FAILURE;
#endif

  return SUCCES;
}


/*==================================================================================================================================*/
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)

int EXTcondition(double *argument, double *dFdp)

{
  const int cfrpntdim = EQUATIONS_DIM + 1;                                          // All variables plus derivative parameter
  int       index     = 0, i, retval;
  double    cfrpnt[EQUATIONS_DIM + 1];
  double    dydp[EQUATIONS_DIM];

  memset(cfrpnt, 0, cfrpntdim*sizeof(double));                                      // Initialize point to 0
  memset(dydp, 0, (cfrpntdim - 1)*sizeof(double));                                  // Initialize dy/dp to 0

  //================================================================================
  // Map current estimate of solution to global variables
  parameter[Bifparone] = argument[index]*pnt_scale[index];
  index++;

  for (i = 0; i < EQUATIONS_DIM; i++)
    {
      Sol[i] = argument[index]*pnt_scale[index];
      index++;
    }

  if (CurveType != EQ) parameter[Bifpartwo] = argument[index]*pnt_scale[index];
  index++;

  for (i = 0; i < EQUATIONS_DIM; i++) cfrpnt[i] = Sol[i];
  cfrpnt[EQUATIONS_DIM] = 1.0;                                                      // Assign relative diff. parameter to last element

  retval = CurveFuncDeriv(cfrpntdim, cfrpnt, dydp, 1, dFdp, DerivEquations, CENTRAL, DYTOL);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  if (checkInterrupt()) return FAILURE;
#endif

  return retval;
}

#endif

/*==================================================================================================================================*/

int AllEquations(double *argument, double *result)

{
  int i, j, retval = SUCCES, savedTestRun;

  // Compute the basic system of equations
  EQsystem(argument, result);

  if (TestRun)
    {
      STDOUT("\n\nParameter #1:\t\t\t\t%18.6G", parameter[Bifparone]);
      for (i = 0; i < EQUATIONS_DIM; i++) STDOUT("\nValue of variable #%d:\t\t\t%18.6G", i, Sol[i]);
      if (CurveType != EQ) STDOUT("\nParameter #2:\t\t\t\t%18.6G", parameter[Bifpartwo]);
      if (CurveType == BP)
        {
          j = i;
          STDOUT("\nValue of eigenvalue:\t\t\t%18.6G", argument[j++]);
          for (i = 0; i < EQUATIONS_DIM; i++) STDOUT("\nValue of eigenvector component #%d:\t%18.6G", i, argument[j++]);
        }
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  if (checkInterrupt()) return FAILURE;
#endif

  //================================================================================
  // Add the final value in case of BP, LP or EXT continuation

  savedTestRun = TestRun;
  TestRun      = 0;

  if (CurveType == LP)
    retval = LPcondition(pntdim, argument, EQsystem, CENTRAL, 1, result + EQUATIONS_DIM, DYTOL);
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  else if (CurveType == EXT)
    retval = EXTcondition(argument, result + EQUATIONS_DIM);
#endif
  else if (CurveType == BP)
    retval = BPconditions(pntdim, argument, EQsystem, CENTRAL, 1, result, RHSTOL);
  else if (CurveType == EQ)
    {
      if (LocalizeType == LP)
        retval = LPcondition(pntdim, argument, EQsystem, CENTRAL, 0, result + EQUATIONS_DIM, DYTOL);
      else if (LocalizeType == BP)
        retval = BPconditions(pntdim, argument, EQsystem, CENTRAL, 0, result, RHSTOL);
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
      else if (LocalizeType == EXT)
        retval = EXTcondition(argument, result + EQUATIONS_DIM);
#endif
    }
  TestRun = savedTestRun;

  if (TestRun)
    {
      STDOUT("\n\n");
      if (CurveType == BP)
        {
          for (i = 0; i < EQUATIONS_DIM; i++) STDOUT("\nValue of equation %2d:\t\t\t%18.6G", i, result[i]);
          j      = i;
          STDOUT("\nInproduct with parameter derivative:\t%18.6G", result[j++]);
          for (i = 0; i < EQUATIONS_DIM; i++) STDOUT("\nInproduct with Jacobian component #%d:\t%18.6G", i, result[j++]);
          STDOUT("\nUnit norm of eigenvector:\t\t%18.6G", result[j++]);
        }
      else
        for (i = 0; i < pntdim - 1; i++) STDOUT("\nValue of equation %2d:\t\t\t%18.6G", i, result[i]);
      STDOUT("\n");
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  if (checkInterrupt()) return FAILURE;
#endif

  return retval;
}


/*==================================================================================================================================*/

int DefineOutput(double *x, double *output)

{
  int     outnr = 0, i;
  double  result[MAX_PNTDIM];
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  double  extra[EXTRAOUTPUT_DIM];
#endif

  AllEquations(x, result);

  // There are maximally (EQUATIONS_DIM+2) values in the point vector
  // to solve for, which occurs when continuing an LP for this system.
  // In all cases the (EQUATIONS_DIM+2) values are written to the output file
  output[outnr++] = parameter[Bifparone];
  for (i = 0; i < EQUATIONS_DIM; i++) output[outnr++] = Sol[i];
  output[outnr++] = parameter[Bifpartwo];

#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  // Also add the additional output variables
  DefineExtraOutput(Sol, extra);

  for (i = 0; i < EXTRAOUTPUT_DIM; i++) output[outnr++] = extra[i];
#endif

  // Should we indeed add how accurate the solution was?
  output[outnr++] = NRM2(pntdim - 1, result, 1);

  return outnr;
}


/*
 *====================================================================================================================================
 *  Implementation of the routine that computes the entire curve
 *====================================================================================================================================
 */

void ComputeCurve(const int argc, char **argv)
{
  register int  i, colnr;
  int           pntnr = 0, outmax;
  int           cycles, last = 1, retval, unknowns;
  int           SkipOutput = 0;
  int           Look4BP = 0, Look4LP = 0;
  double        oldpoint[MAX_PNTDIM], tmpVec[MAX_PNTDIM];
  double        tanvec[MAX_PNTDIM], oldvec[MAX_PNTDIM];
  double        detJ;
  char          bifname[MAX_STR_LEN], errname[MAX_STR_LEN], outname[MAX_STR_LEN];
  struct stat   buffer;
#if (DETECTSPECIALPOINTS == 1)
  double        oldBPval = 0, oldLPval = 0;
  double        lastBPpoint[MAX_PNTDIM], lastLPpoint[MAX_PNTDIM];
#endif
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
  int           Look4EXT  = 0;
  double        fncderiv, oldEXTval = 0, lastEXTpoint[MAX_PNTDIM];
#endif

  unknowns = (CurveType == BP) ? (pntdim + 1 + EQUATIONS_DIM) : pntdim;
  for (i = 0; i < MAX_PNTDIM; i++) point_scales[i] = 1.0;
  pnt_scale = point_scales;
  (void)SetScales(point, pntdim);
  CurveEnd = 0;

  // Initialize additional variables in BP continuation
  if (CurveType == BP) BPconditions(pntdim, point, EQsystem, CENTRAL, 1, NULL, RHSTOL);

  if (TestRun)
    {
      double rhs[MAX_PNTDIM], rhsnorm;

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
      STDOUT("Executing : ");
      for (i = 0; i < argc; i++) STDOUT("%s ", argv[i]);
#else
      STDOUT("\n\nExecuting : ");
      STDOUT("FindCurve('%s', '%s', %s, %.6G, %s, %s, %s)", progname, ContinuationString, pntstring, curvestep, curvestring, parstring, optstring);
#endif
      STDOUT("\n\n");

      STDOUT("Parameter values  : \n");
      for (i = 0; i < PARAMETER_NR; i++)
        {
          if (!(i % 3)) STDOUT("\n");
          STDOUT("\t%-10s:", parameternames[i]);
          STDOUT("  %-13G", parameter[i]);
        }
      STDOUT("\n");
      fflush(NULL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif

      AllEquations(point, rhs);

      // Compute rhsnorm
      rhsnorm = NRM2(unknowns - 1, rhs, 1);
      rhsnorm = rhsnorm/(1.0 + rhsnorm);

      STDOUT("\n\nNorm of fixed point conditions:\t\t%18.6G", rhsnorm);

      outmax = DefineOutput(point, Output);

      STDOUT("\n\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      STDOUT("#      %2d: Par.1", colnr++);
      for (i = 0; i < EQUATIONS_DIM; i++) STDOUT("       %2d: x[%2d]", colnr++, i);
      STDOUT("       %2d: Par.2", colnr++);
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
      for (i = 0; i < EXTRAOUTPUT_DIM; i++) STDOUT("       %2d: O[%2d]", colnr++, i);
#endif
      STDOUT("    %2d: RHS norm\n", colnr++);
      STDOUT("#\n");
      if (outmax) PrettyPrintArray(NULL, outmax, Output);

#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif

      return;
    }

  if (strlen(runname))
    {
      sprintf(bifname, "%s.bif", runname);
      sprintf(errname, "%s.err", runname);
      sprintf(outname, "%s.out", runname);
    }
  else
    {
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)
      strcpy(progname, argv[0]);
#endif
      i = 0;
      while (1)
        {
          sprintf(bifname, "%s-%s-%04d.bif", progname, ContinuationString, i);
          sprintf(errname, "%s-%s-%04d.err", progname, ContinuationString, i);
          sprintf(outname, "%s-%s-%04d.out", progname, ContinuationString, i);
          if (stat(bifname, &buffer) && stat(errname, &buffer) && stat(outname, &buffer)) break;
          i++;
        }
      sprintf(runname, "%s-%s-%04d", progname, ContinuationString, i);
    }

  if (CurveType == EQ) biffile = fopen(bifname, "w");
  errfile                      = fopen(errname, "w");
  outfile                      = fopen(outname, "w");

  if (outfile)
    {
      fprintf(outfile, "#\n# Executing : ");

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      fprintf(outfile, "FindCurve('%s', '%s', %s, %.6G, %s, %s, %s)", progname, ContinuationString, pntstring, curvestep, curvestring, parstring, optstring);
#elif defined(R_PACKAGE)
      fprintf(outfile, "FindCurve(\"%s\", \"%s\", %s, %.6G, %s, %s, %s)", progname, ContinuationString, pntstring, curvestep, curvestring, parstring, optstring);
#else
      for (i = 0; i < argc; i++) fprintf(outfile, "%s ", argv[i]);
#endif
      fprintf(outfile, "\n#\n");

      fprintf(outfile, "# Parameter values  : \n#");
      for (i = 0; i < PARAMETER_NR; i++)
        {
          if (!(i % 3)) fprintf(outfile, "\n# ");
          fprintf(outfile, "\t%-10s:", parameternames[i]);
          fprintf(outfile, "  %-13G", parameter[i]);
        }
      fprintf(outfile, "\n#\n");
      fprintf(outfile, "# Index of bifurcation parameter #1               : %d\n", Bifparone);
      if (CurveType != EQ) fprintf(outfile, "# Index of bifurcation parameter #2               : %d\n", Bifpartwo);
      fprintf(outfile, "#\n");
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)   // In command-line model follow C convention of 0 start index
      colnr = 0;
#else
      colnr = 1;
#endif
      fprintf(outfile, "#      %2d: Par.1", colnr++);
      for (i = 0; i < EQUATIONS_DIM; i++) fprintf(outfile, "       %2d: x[%2d]", colnr++, i);
      fprintf(outfile, "       %2d: Par.2", colnr++);
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
      for (i = 0; i < EXTRAOUTPUT_DIM; i++) fprintf(outfile, "       %2d: O[%2d]", colnr++, i);
#endif
      fprintf(outfile, "    %2d: RHS norm\n", colnr++);
      fprintf(outfile, "#\n");
      fflush(outfile);
    }

  memset((void *)tanvec, 0, unknowns*sizeof(double));
  tanvec[0] = 1.0;

  // Continue the curve
  while (1)
    {
      // Compute fixed point with new varied parameter
      cycles = 0;
      retval = SUCCES;
      while (Stepreduce <= MAX_STEPREDUCE)
        {
          retval = FindPoint(unknowns, point, NULL, tanvec, DYTOL, RHSTOL, MAXITER, AllEquations);

          if (retval == SUCCES) break;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
#endif

          NumProcError(__FILE__, __LINE__, retval);

          // If unsuccesfull while repeating point exit
          if (!Stepchange)
            {
              ErrorMsg(__FILE__, __LINE__, "Failed to locate a solution point after scaling\n");
#if (defined(R_PACKAGE))
              R_FlushConsole();
              R_ProcessEvents();
#endif
              return;
            }

          // Generate prediction of solution point with smaller step size
          cycles++;
          Stepreduce *= 4;
          COPY(pntdim, oldpoint, 1, point, 1);

          //   AXPY(pntdim, (curvestep/(Stepreduce*pnt_scale[0])), tanvec, 1, point, 1);
          AXPY(pntdim, (curvestep/Stepreduce), tanvec, 1, point, 1);

          // Initialize additional variables in BP continuation
          if (CurveType == BP) BPconditions(pntdim, point, EQsystem, CENTRAL, 1, NULL, RHSTOL);

          ReportMsg("\n\nPrediction :\t");
          for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", point[i]*pnt_scale[i]);
          ReportMsg("\n");
        }

      // If unsuccesfull exit
      if (retval != SUCCES)
        {
          ErrorMsg(__FILE__, __LINE__, "Failed to locate a solution point\n");
          break;
        }

      // New solution point
      ReportMsg("\nNew point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", point[i]*pnt_scale[i]);
      ReportMsg("\n");

      if (!DoSingle)
        {
          // Increase step when both this and previous point were located with
          // the current step size
          if ((last && !cycles) && (Stepreduce > 1)) Stepreduce /= 2;
          last = (!cycles);

          // Save the old tangent vector and compute the new tangent vector
          // Ignore the additional variables in a BP continuation
          retval = TangentVec(unknowns, point, NULL, tanvec, AllEquations, &detJ, DYTOL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
          if (checkInterrupt()) return;
#endif

          // Signal curve stop if one of the components has become negative or parameter is out of bounds
          // and this is not the curve beginning
          if (pntnr > 20)
            {
              for (i = 0; i < pntdim; i++)
                {
                  CurveEnd = CurveEnd || ((point[i]*pnt_scale[i]) < pntmin[i] - DYTOL) || ((point[i]*pnt_scale[i]) > pntmax[i] + DYTOL);
                  if (!AllowNegative) CurveEnd = CurveEnd || (point[i] <= -DYTOL);
                }
            }

          // Sanitize the point:
          for (i = 0; i < pntdim; i++)
            if (fabs(point[i]*pnt_scale[i]) < 1.0E-13) point[i] = 0.0;

#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
          retval = EXTcondition(point, &fncderiv);
#endif

          /*
           * ===========================================================================
           * The following code section is the only section that is specific to the
           * problem of continuation of structured population curves. All other parts of
           * the main routine are generic. The following lines detect bifurcation points
           */
#if (DETECTSPECIALPOINTS == 1)
          if (CurveType == EQ)
            {
              if (BPdetection && Look4BP && (detJ*oldBPval < 0))
                {
                  if (fabs(detJ) < fabs(oldBPval))
                    COPY(pntdim, point, 1, tmpVec, 1);
                  else
                    COPY(pntdim, lastBPpoint, 1, tmpVec, 1);
                  retval = LocateBP(pntdim, tmpVec, AllEquations, DYTOL, RHSTOL);
                  if (retval == SUCCES)
                    {
                      Look4BP    = 0;
                      SkipOutput = 1;
                    }
                }
              else if (LPdetection && Look4LP && (tanvec[0]*oldLPval < 0))
                {
                  if (fabs(tanvec[0]) < fabs(oldLPval))
                    COPY(pntdim, point, 1, tmpVec, 1);
                  else
                    COPY(pntdim, lastLPpoint, 1, tmpVec, 1);
                  retval = LocateLP(pntdim, tmpVec, AllEquations, DYTOL, RHSTOL);
                  if (retval == SUCCES)
                    {
                      Look4LP    = 0;
                      SkipOutput = 1;
                    }
                }
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
              else if (EXTdetection && Look4EXT && (retval == SUCCES) && (fncderiv*oldEXTval < 0))
                {
                  if (fabs(fncderiv) < fabs(oldEXTval))
                    COPY(pntdim, point, 1, tmpVec, 1);
                  else
                    COPY(pntdim, lastEXTpoint, 1, tmpVec, 1);
                  retval = LocateEXT(pntdim, tmpVec, AllEquations, DYTOL, RHSTOL);
                  if (retval == SUCCES)
                    {
                      Look4EXT   = 0;
                      SkipOutput = 1;
                    }
                }
#endif
              if (BPdetection && (fabs(detJ) > RHSTOL))
                {
                  Look4BP  = 1;
                  oldBPval = detJ;
                  COPY(pntdim, point, 1, lastBPpoint, 1);
                }
              if (LPdetection && (fabs(tanvec[0]) > RHSTOL))
                {
                  Look4LP  = 1;
                  oldLPval = tanvec[0];
                  COPY(pntdim, point, 1, lastLPpoint, 1);
                }
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
              if (EXTdetection && (fabs(fncderiv) > RHSTOL))
                {
                  Look4EXT  = 1;
                  oldEXTval = fncderiv;
                  COPY(pntdim, point, 1, lastEXTpoint, 1);
                }
#endif
            }
#endif
        }

      if ((!SkipOutput) && ((pntnr % ReportLevel) == 0))
        {
          // Generate output: invoked after setting CurveEnd to allow for output of last solution point on branch
          for (i = 0; i < pntdim; i++)
            {
              if (!i) STDOUT("%16.8E", point[i]*pnt_scale[i]);
              else STDOUT(", %16.8E", point[i]*pnt_scale[i]);
            }
          STDOUT("\n");
#if (defined(R_PACKAGE))
          R_FlushConsole();
          R_ProcessEvents();
#endif

          DoOutput = 1;
          outmax = DefineOutput(point, Output);
          DoOutput = 0;
          if (outmax) PrettyPrintArray(outfile, outmax, Output);
        }
      SkipOutput = 0;

      // Exit at end of curve after generation of last output point
      if (CurveEnd | DoSingle) break;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      // check for a Ctrl-C event
      if (checkInterrupt()) return;
#endif

      pntnr++;

      // Scale the point vector anew if necessary and redo the current point
      if ((retval = SetScales(point, pntdim)))
        {
          ReportMsg("\n\nVariable %d rescaled!\n", retval);

          // Point is rescaled: Save the previously computed tangent vector in the current solution point
          // to preserve the direction after rescaling
          // Ignore the additional variables in a BP continuation
          COPY(pntdim, tanvec, 1, tmpVec, 1);

          // Compute the new tangent vector. If the inner product with the previous tangent vector is
          // negative, the direction has been reversed and is hence not preserved. Therefore, scale the
          // tangent vector with a factor -1 to flip the direction.
          // Ignore the additional variables in a BP continuation
          retval = TangentVec(unknowns, point, NULL, tanvec, AllEquations, &detJ, DYTOL);
          if (DOT(pntdim, tmpVec, 1, tanvec, 1) < 0) SCAL(unknowns, -1, tanvec, 1);
          Look4BP    = 0;
          Look4LP    = 0;
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
          Look4EXT   = 0;
#endif
          Stepchange = 0;
          SkipOutput = 1;
        }
      // Otherwise generate output and predict new point on the curve
      else
        {
          // Store info on current solution point
          COPY(unknowns, point, 1, oldpoint, 1);
          COPY(pntdim, tanvec, 1, oldvec, 1);
          Stepchange = 1;

          ReportMsg("Determinant of Jacobian         : %.8G\n", detJ);
#if defined(EXTRAOUTPUT_DIM) && (EXTRAOUTPUT_DIM >= 1)
          ReportMsg("Derivative of objective function: %.8G\n", fncderiv);
#endif
          // Determine appropriate stepsize
          ReportMsg("Tangent vector in component    0: %.8G\n", tanvec[0]);
          ReportMsg("Targeted  step in component    0: %.8G\n", curvestep*tanvec[0]/Stepreduce);

          if ((Stepreduce == 1) && (fabs(curvestep) > fabs(Maxcurvestep))) curvestep = sign(curvestep)*fabs(Maxcurvestep);
          AXPY(pntdim, (curvestep/Stepreduce), tanvec, 1, point, 1);

          // Initialize additional variables in BP continuation
          if (CurveType == BP) BPconditions(pntdim, point, EQsystem, CENTRAL, 1, NULL, RHSTOL);

          ReportMsg("Realized  step in component    0: %.8G\n", curvestep*tanvec[0]/Stepreduce);

          // Prediction
          ReportMsg("\n\nPrediction :\t");
          for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", point[i]*pnt_scale[i]);
          ReportMsg("\n");
        }
      fflush(NULL);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
      mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  STDOUT("\n");
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif
  return;
}


/*==================================================================================================================================*/

void InitialiseVars(void)

{
  int                             i;

#if (defined(_MSC_VER) && (_MSC_VER < 1500)) || (defined(R_PACKAGE) && defined(_WIN32))
  (void)_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  // Initialize some variables
  errfile    = NULL;
  outfile    = NULL;
  Stepchange = 0;
  Stepreduce = 1;
  strcpy(runname, "");

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
  CtrlCPressed = 0;
#endif

  CurveType = UNDEFINED;
  ConsoleSilent = 0;
  DefaultBifparone = Bifparone;
  DefaultBifpartwo = Bifpartwo;
  DoOutput = 0;
  for (i = 0; i < MAX_PNTDIM; i++) pntmin[i] = -DBL_MAX;
  for (i = 0; i < MAX_PNTDIM; i++) pntmax[i] = DBL_MAX;

  // Get the machine precisions
  epsMach = dlamch("Epsilon" FCONE);
  
#if defined(R_PACKAGE)
  STDOUT("\n");
#else
  fprintf(stderr, "\n");
#endif

  Jacobian_Min_Step   = JACOBIAN_MIN_STEP;
  Jacobian_Step       = JACOBIAN_STEP;
  Jacobian_Updates    = JACOBIAN_UPDATES;
  odesolveMem         = NULL;

  // Set the integration options
  Odesolve_Init_Step  = ODESOLVE_INIT_STEP;
  Odesolve_Fixed_Step = ODESOLVE_FIXED_STEP;
  Odesolve_Min_Step   = ODESOLVE_MIN_STEP;
  Odesolve_Max_Step   = ODESOLVE_MAX_STEP;
  Odesolve_Abs_Err    = ODESOLVE_ABS_ERR;
  Odesolve_Rel_Err    = ODESOLVE_REL_ERR;
  Odesolve_Func_Tol   = ODESOLVE_FUNC_TOL;

  return;
}


/*
 *====================================================================================================================================
 *  UNIX shell interface function main() and supporting functions.
 *====================================================================================================================================
 */

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE) && !defined(R_PACKAGE)

static void Usage(char *progname)

{
  int   i;
  char  tmpstr[MAX_STR_LEN];
  char  bifstr[MAX_STR_LEN];
  char  desc[MAX_STR_LEN];
  char  varstr[MAX_STR_LEN];

  switch (CurveType)
    {
      case BP:
        strcpy(bifstr, "BP");
        strcpy(desc,
               "Aim:\tContinuation of a transcritical bifurcation in a curve determined by a system of equations\n\tas a function of two parameters");
        break;
      case EQ:
        strcpy(bifstr, "EQ");
        strcpy(desc, "Aim:\tContinuation of a curve determined by a system of equations as a function of a single parameter");
        break;
      case LP:
        strcpy(bifstr, "LP");
        strcpy(desc,
               "Aim:\tContinuation of a saddle-node bifurcation in a curve determined by a system of equations\n\tas a function of two parameters");
        break;
      case EXT:
        strcpy(bifstr, "EXT");
        strcpy(desc, "Aim:\tContinuation of a maximum or minimum function value along a curve determined by a system of equations\n\tas a function "
                     "of two parameters");
        break;
      default:
        strcpy(bifstr, "<BP, EQ, LP or EXT>");
        strcpy(desc, "Aim:\tContinuation of a curve determined by a system of equations, as well as transcritical and saddle-node bifurcations\n\t");
        strcat(desc, "or a maximum or minimum function value along this curve as a function of one or two parameters");
    }
  if (CurveType == UNDEFINED)
    strcpy(varstr, "<Initial values of all variables>");
  else
    {
      strcpy(varstr, "Par.1");
      for (i = 0; i < EQUATIONS_DIM; i++)
        {
          sprintf(tmpstr, " x[%d]", i);
          strcat(varstr, tmpstr);
        }
      if (!(CurveType == EQ)) strcat(varstr, " Par.2");
    }

  fprintf(stderr, "\nUsage:\t%s [<options>] %s %s %s", progname, bifstr, varstr, "<step> <min. par.1> <max. par.1>");
  if (!(CurveType == EQ)) fprintf(stderr, " <min. par.2> <max. par.2>");
  fprintf(stderr, "\n\n%s\n\n", desc);
  fprintf(stderr, "Possible options are:\n\n");
  fprintf(stderr, "\t-par1   <index>: Index of first bifurcation parameter\n");
  fprintf(stderr, "\t-par2   <index>: Index of second bifurcation parameter\n");
  fprintf(stderr, "\t-EXTfun <index>: Index of the element in the ExtraOutput[] vector, for which to test for maximum and minimum values along the curve\n");
  fprintf(stderr, "\t-EXTpar <index>: Index of the parameter, with respect to which to test for maximum and minimum values along the curve\n");
  fprintf(stderr, "\t-report <value>: Interval of reporting computed output to console. Minimum value of 1 implies output of every point.\n");
  fprintf(stderr, "\t-negative      : Allow negative values in the solution point\n");
  fprintf(stderr, "\t-noBP          : Do not check for branching points while computing equilibrium curves\n");
  fprintf(stderr, "\t-noEXT         : Do not check for extremum points while computing equilibrium curves\n");
  fprintf(stderr, "\t-noLP          : Do not check for limit points while computing equilibrium curves\n");
  fprintf(stderr, "\t-silent        : Do not report error messages to the console\n");
  fprintf(stderr, "\t-single        : Only compute the first point of the solution curve, do not continue the curve\n");
  fprintf(stderr, "\t-test          : Perform only a single integration over the life history, reporting dynamics of survival, R0, i-state and "
                  "interaction variables\n");
  fprintf(stderr, "\nThe values for -par1, -par2, -EXTfun and -EXTpar default to 0\n");
  fprintf(stderr, "\n");
  exit(1);

  return;
}


/*
 *====================================================================================================================================
 *  Implementation of the generic main routine.
 *====================================================================================================================================
 */

int main(int argc, char **argv)
{
  register int  i, j;
  int           my_argc, tmpint;
  char          **argpnt1 = NULL, **argpnt2 = NULL, *my_argv[argc];

  if (argc < 2) Usage(argv[0]);

  InitialiseVars();

  argpnt1 = argv;
  argpnt2 = my_argv;
  my_argc = 0;

  // Store the program name
  *argpnt2 =*argpnt1;
  my_argc++;
  argpnt1++;
  argpnt2++;

  while (*argpnt1)
    {
      if (!strcmp(*argpnt1, "-?") || !strcmp(*argpnt1, "--help"))
        {
          Usage(argv[0]);
        }
      else if (!strcmp(*argpnt1, "-par1"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "No index of first bifurcation parameter specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            {
              fprintf(stderr, "Index of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
              Usage(argv[0]);
            }
          Bifparone = tmpint;
        }
      else if (!strcmp(*argpnt1, "-par2"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "No index of second bifurcation parameter specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            {
              fprintf(stderr, "Index of second parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
              Usage(argv[0]);
            }
          Bifpartwo = tmpint;
        }
      else if (!strcmp(*argpnt1, "-EXTfun"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          fprintf(stderr, "EXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
          Usage(argv[0]);
#endif
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "Index of the element in the ExtraOutput[] vector, for which to test for maximum and minimum values along the curve, "
                              "not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= EXTRAOUTPUT_DIM))
            {
              fprintf(stderr, "Index of the element in the ExtraOutput[] vector for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                      tmpint, EXTRAOUTPUT_DIM);
              Usage(argv[0]);
            }
          EXTfun = tmpint;
        }
      else if (!strcmp(*argpnt1, "-EXTpar"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          fprintf(stderr, "EXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
          Usage(argv[0]);
#endif
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "Index of the parameter, with respect to which to test for maximum and minimum values along the curve, not specified!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            {
              fprintf(stderr, "Index of the parameter for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
              Usage(argv[0]);
            }
          EXTpar = tmpint;
        }
      else if (!strcmp(*argpnt1, "-report"))
        {
          argpnt1++;
          if (!*argpnt1)
            {
              fprintf(stderr, "\nNo index of i-state variable specified for argument -report!\n");
              Usage(argv[0]);
            }
          tmpint = atoi(*argpnt1);
          ReportLevel = max(tmpint, 1);
        }
      else if ((!strcmp(*argpnt1, "-negative")))
        {
          AllowNegative = 1;
        }
      else if ((!strcmp(*argpnt1, "-noBP")))
        {
          BPdetection = 0;
        }
      else if ((!strcmp(*argpnt1, "-noEXT")))
        {
          EXTdetection = 0;
        }
      else if ((!strcmp(*argpnt1, "-noLP")))
        {
          LPdetection = 0;
        }
      else if ((!strcmp(*argpnt1, "-silent")))
        {
          ConsoleSilent = 1;
        }
      else if ((!strcmp(*argpnt1, "-single")))
        {
          DoSingle = 1;
        }
      else if ((!strcmp(*argpnt1, "-test")))
        {
          TestRun = 1;
        }
      else if ((!strncmp(*argpnt1, "--", 2)))
        {
          fprintf(stderr, "Unknown command line option: %s\n", *argpnt1);
          Usage(argv[0]);
        }
      else if ((!strncmp(*argpnt1, "-", 1)) && isalpha(*(*argpnt1 + 1)))
        {
          fprintf(stderr, "Unknown command line option: %s\n", *argpnt1);
          Usage(argv[0]);
        }
      else
        {
          *argpnt2 =*argpnt1;
          my_argc++;
          argpnt2++;
        }
      argpnt1++;
    }

  if (!strcmp(*(my_argv + 1), "BP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = BP;
      strcpy(ContinuationString, *(my_argv + 1));
    }
  else if (!strcmp(*(my_argv + 1), "EQ"))
    {
      pntdim    = EQUATIONS_DIM + 1;
      CurveType = EQ;
      strcpy(ContinuationString, *(my_argv + 1));
    }
  else if (!strcmp(*(my_argv + 1), "LP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = LP;
      strcpy(ContinuationString, *(my_argv + 1));
    }
  else if (!strcmp(*(my_argv + 1), "EXT"))
    {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
      fprintf(stderr, "EXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
      Usage(argv[0]);
#endif
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = EXT;
      strcpy(ContinuationString, *(my_argv + 1));
    }
  else
    {
      fprintf(stderr, "Curve type undefined! Define curve as EQ, BP, LP or EXT!\n");
      Usage(my_argv[0]);
    }

  if (Bifpartwo == Bifparone)
    {
      if (((CurveType == BP) || (CurveType == LP) || (CurveType == EXT)))
        {
          fprintf(stderr, "Index of first and second bifurcation parameter the same!\n");
          fprintf(stderr, "Two parameter continuation not possible!\n");
          Usage(argv[0]);
        }
      // Here we only end up in case of EQ continuation
      // Bifparone reset on command-line to Bifpartwo
      if ((Bifparone != DefaultBifparone) && (Bifpartwo != DefaultBifparone))
        Bifpartwo = DefaultBifparone;
      else if ((Bifparone != DefaultBifpartwo) && (Bifpartwo != DefaultBifpartwo))
        Bifpartwo = DefaultBifpartwo;
      else
        {
          for (i = 0; i < PARAMETER_NR; i++)
            if (i != Bifparone) break;
          Bifpartwo = i;
        }
    }

  if ((CurveType == EQ) && (my_argc != pntdim + 5))
    {
      fprintf(stderr, "Wrong number of command-line values (%d) for curve type EQ, should be %d.\n", my_argc, pntdim + 5);
      Usage(argv[0]);
    }
  else if ((CurveType != EQ) && (my_argc != pntdim + 7))
    {
      fprintf(stderr, "Wrong number of command-line values (%d) for curve type BP, LP or EXT, should be %d.\n", my_argc, pntdim + 7);
      Usage(argv[0]);
    }

  // Map all command-line variables into the argument vector
  memset((void *)point, 0, MAX_PNTDIM*sizeof(double));
  for (i = 0, j = 2; i < pntdim; i++, j++) point[i] = atof(my_argv[j]);
  parameter[Bifparone]                      = point[0];
  if (CurveType != EQ) parameter[Bifpartwo] = point[pntdim - 1];

  Maxcurvestep = curvestep = atof(my_argv[j++]);

  pntmin[0] = atof(my_argv[j++]);
  pntmax[0] = atof(my_argv[j++]);

  if (CurveType != EQ)
    {
      pntmin[pntdim - 1] = atof(my_argv[j++]);
      pntmax[pntdim - 1] = atof(my_argv[j++]);
    }

  ComputeCurve(argc, argv);

  return 0;
}


/*==================================================================================================================================*/
#elif defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

static void CloseStreams(void)
{
  fflush(NULL);
  if (biffile) fclose(biffile);
  if (errfile) fclose(errfile);
  if (outfile) fclose(outfile);

  ResetCurve();
  if (odesolveMem) free(odesolveMem);
  odesolveMem = NULL;

  return;
}

// The gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t          nrows, ncols;
  double          curveVals[2*MAX_PNTDIM], tmpdouble;
  int             tmpint, irhs;
  const mxArray   *cell_element_ptr;
  char            optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstring[MAX_STR_LEN];
  mwIndex         i, j;
  size_t          total_num_of_cells, buflen;
  int             status;

  InitialiseVars();

  // check for proper number of arguments
  if (nrhs != 6)
    mexErrMsgIdAndTxt("MATLAB:FindCurve:nrhs",
                      "\nIncorrect number of command-line arguments.\n\nUse: %s(%s)\n\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s\n%12s: %s",
                      mexFunctionName(), "<type>, <point>, <step>, <bounds>, <parameters>, <options>", "<type>",
                      "Type of bifurcation to perform (BP, EQ, LP or EXT)", "<point>", "Initial point of the bifurcation", "<step>",
                      "Step size in bifurcation", "<bounds>", "Minimum and maximum parameter value for the bifurcation (2 or 4 values)",
                      "<parameters>", "Array of parameter values to use (empty array or PARAMETER_NR in length)", "<options>",
                      "Possible bifurcation options: par1, par2, EXTfun, EXTpar, negative, noBP, noEXT, noLP, single, silent or test");

  // check for proper number of output variables
  if (nlhs != 1) mexErrMsgIdAndTxt("MATLAB:FindCurve:nlhs", "A single output argument is required.");

  //============================== Process the options argument ======================================================================
  // Extract the contents of MATLAB cell into the C array
  irhs = 5;
  if (!mxIsCell(prhs[irhs])) mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nOptions should be specified as a cell array!\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(optstring, "{");
  for (i = 0; i < total_num_of_cells; i++)
    {
      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optname, buflen);
      if (!((!strcmp(optname, "par1")) || (!strcmp(optname, "par2")) || (!strcmp(optname, "EXTfun")) || (!strcmp(optname, "EXTpar")) ||
          (!strcmp(optname, "negative")) || (!strcmp(optname, "noBP")) || (!strcmp(optname, "noEXT")) || (!strcmp(optname, "noLP")) ||
          (!strcmp(optname, "single")) || (!strcmp(optname, "silent")) || (!strcmp(optname, "test")) || (!strcmp(optname, "report"))))
        mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nIllegal option %s!\n", optname);

      if (!strcmp(optname, "negative") || !strcmp(optname, "noBP") || !strcmp(optname, "noEXT") || !strcmp(optname, "noLP") ||
          !strcmp(optname, "single") || !strcmp(optname, "silent") || !strcmp(optname, "test"))
        {
        if (!strcmp(optname, "negative")) AllowNegative = 1;
        if (!strcmp(optname, "noBP")) BPdetection = 0;
        if (!strcmp(optname, "noEXT")) EXTdetection = 0;
        if (!strcmp(optname, "noLP")) LPdetection = 0;
        if (!strcmp(optname, "silent")) ConsoleSilent  = 1;
        if (!strcmp(optname, "single")) DoSingle  = 1;
        if (!strcmp(optname, "test")) TestRun     = 1;

          // optstring still equal to "{"
          if (strlen(optstring) == 1)
            strcat(optstring, "'");
          else
            strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "'");
          continue;
        }

      if (!(++i < total_num_of_cells)) mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nNo value specified for option %s!\n", optname);

      cell_element_ptr = mxGetCell(prhs[irhs], i);
      buflen           = mxGetN(cell_element_ptr)*sizeof(mxChar) + 1;
      status           = mxGetString(cell_element_ptr, optval, buflen);
      if (status) mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nError in retrieving value for option %s!\n", optname);

      // optstring still equal to "{"
      if (strlen(optstring) == 1)
        strcat(optstring, "'");
      else
        strcat(optstring, ", '");
      strcat(optstring, optname);
      strcat(optstring, "', '");
      strcat(optstring, optval);
      strcat(optstring, "'");

      if (!strcmp(optname, "par1"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nIndex of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          Bifparone = tmpint;
        }
      else if (!strcmp(optname, "par2"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            mexErrMsgIdAndTxt("MATLAB:FindCurve:options", "\nIndex of second parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          Bifpartwo = tmpint;
        }
      else if (!strcmp(optname, "EXTfun"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          mexErrMsgIdAndTxt(
                "MATLAB:FindCurve:options",
                "\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EXTRAOUTPUT_DIM))
            mexErrMsgIdAndTxt(
                "MATLAB:FindCurve:options",
                "\nIndex of the element in the ExtraOutput[] vector for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, EXTRAOUTPUT_DIM);
          EXTfun = tmpint;
        }
      else if (!strcmp(optname, "EXTpar"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          mexErrMsgIdAndTxt(
                "MATLAB:FindCurve:options",
                "\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            mexErrMsgIdAndTxt("MATLAB:FindCurve:options",
                              "\nIndex of the parameter for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          EXTpar = tmpint;
        }
      else if (!strcmp(optname, "report"))
        {
          tmpint = atoi(optval);
          ReportLevel = max(tmpint, 1);
        }
    }
  strcat(optstring, "}");

  //============================== Process the bifurcation type argument =============================================================
  // Input must be a string
  irhs = 0;
  if ((mxIsChar(prhs[irhs]) != 1) || (mxGetM(prhs[irhs]) != 1))
    mexErrMsgIdAndTxt("MATLAB:FindCurve:inputNotString", "Bifurcation type must be a string (BP, EQ, LP or EXT).");

  // Copy the string data from prhs[irhs] into a C string input_ buf.
  strcpy(ContinuationString, mxArrayToString(prhs[irhs]));
  if (!strcmp(ContinuationString, "BP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = BP;
    }
  else if (!strcmp(ContinuationString, "EQ"))
    {
      pntdim    = EQUATIONS_DIM + 1;
      CurveType = EQ;
    }
  else if (!strcmp(ContinuationString, "LP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = LP;
    }
  else if (!strcmp(ContinuationString, "EXT"))
    {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
      mexErrMsgIdAndTxt(
            "MATLAB:FindCurve:curveType",
            "\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = EXT;
    }

  if (Bifpartwo == Bifparone)
    {
      if (((CurveType == BP) || (CurveType == LP) || (CurveType == EXT)))
        mexErrMsgIdAndTxt("MATLAB:FindCurve:bifpars",
                          "\nIndex of first and second bifurcation parameter the same!\nTwo parameter continuation not possible!\n");

      // Here we only end up in case of EQ continuation
      // Bifparone reset on command-line to Bifpartwo
      if ((Bifparone != DefaultBifparone) && (Bifpartwo != DefaultBifparone))
        Bifpartwo = DefaultBifparone;
      else if ((Bifparone != DefaultBifpartwo) && (Bifpartwo != DefaultBifpartwo))
        Bifpartwo = DefaultBifpartwo;
      else
        {
          for (i = 0; i < PARAMETER_NR; i++)
            if (i != Bifparone) break;
          Bifpartwo = (int)i;
        }
      mexWarnMsgIdAndTxt("MATLAB:FindCurve:bifpars", "\nBifurcation parameter #2 reset to parameter #%d for locating branching points\n", Bifpartwo);
    }

  //================================ Process the parameters argument =================================================================

  irhs  = 4;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((ncols == PARAMETER_NR) && (nrows == 1))
    memcpy(parameter, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  else if (ncols)
    mexWarnMsgIdAndTxt("MATLAB:FindCurve:parameters", "\nParameter argument ignored as it is not a row vector of length PARAMETER_NR\n");

  total_num_of_cells = mxGetNumberOfElements(prhs[irhs]);
  strcpy(parstring, "[");
  for (i = 0; i < total_num_of_cells; i++)
    {
      if (i) strcat(parstring, " ");
      memcpy(&tmpdouble, mxGetPr(prhs[irhs]) + i, mxGetElementSize(prhs[irhs]));
      sprintf(tmpstring, "%.6G", tmpdouble);
      strcat(parstring, tmpstring);
    }
  strcat(parstring, "]");

  //============================== Process the curve parameters argument =============================================================

  irhs  = 3;
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if (CurveType == EQ)
    {
      if ((ncols == 2) && (nrows == 1))
        {
          memcpy(curveVals, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
          pntmin[0] = curveVals[0];
          pntmax[0] = curveVals[1];
        }
      else if ((ncols == 2*pntdim) && (nrows == 1))
        {
          memcpy(curveVals, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
          for (i = 0, j = 0; i < ncols; i++, j++)
            {
              pntmin[j] = curveVals[i];
              i++;
              pntmax[j] = curveVals[i];
            }
        }
      else
        mexErrMsgIdAndTxt("MATLAB:FindCurve:bounds", "\nFor %s continuation bounds argument must be a row vector of length 2 or length %d.\n",
                          ContinuationString, 2*pntdim);
    }
  else
    {
      if ((ncols == 4) && (nrows == 1))
        {
          memcpy(curveVals, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
          pntmin[0] = curveVals[0];
          pntmax[0] = curveVals[1];

          pntmin[pntdim - 1] = curveVals[2];
          pntmax[pntdim - 1] = curveVals[3];
        }
      else if ((ncols == 2*pntdim) && (nrows == 1))
        {
          memcpy(curveVals, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
          for (i = 0, j = 0; i < ncols; i++, j++)
            {
              pntmin[j] = curveVals[i];
              i++;
              pntmax[j] = curveVals[i];
            }
        }
      else
        mexErrMsgIdAndTxt("MATLAB:FindCurve:bounds", "\nFor %s continuation bounds argument must be a row vector of length 4 or length %d.\n",
                          ContinuationString, 2*pntdim);
    }

  strcpy(curvestring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(curvestring, " ");
      sprintf(tmpstring, "%.6G", curveVals[i]);
      strcat(curvestring, tmpstring);
    }
  strcat(curvestring, "]");

  //================================ Process the step size argument ==================================================================
  // Get the step size, make sure it is scalar
  irhs = 2;
  if (!mxIsDouble(prhs[irhs]) || mxIsComplex(prhs[irhs]) || mxGetNumberOfElements(prhs[irhs]) != 1)
    mexErrMsgIdAndTxt("MATLAB:FindCurve:stepsize", "\nStep size must be a scalar value.\n");
  Maxcurvestep = mxGetScalar(prhs[irhs]);
  curvestep    = Maxcurvestep;

  //============================== Process the initial point argument ================================================================

  irhs = 1;
  memset((void *)point, 0, MAX_PNTDIM*sizeof(double));
  nrows = mxGetM(prhs[irhs]);
  ncols = mxGetN(prhs[irhs]);
  if ((nrows != 1) || (ncols != pntdim))
    mexErrMsgIdAndTxt("MATLAB:FindCurve:point", "\nInitial point for bifurcation must be a row vector of length %d\n", pntdim);
  memcpy(point, mxGetPr(prhs[irhs]), ncols*mxGetElementSize(prhs[irhs]));
  parameter[Bifparone]                      = point[0];
  if (CurveType != EQ) parameter[Bifpartwo] = point[pntdim - 1];

  strcpy(pntstring, "[");
  for (i = 0; i < (nrows*ncols); i++)
    {
      if (i) strcat(pntstring, " ");
      sprintf(tmpstring, "%.6G", point[i]);
      strcat(pntstring, tmpstring);
    }
  strcat(pntstring, "]");

  //============================= Get the program name ===============================================================================
  // Get the name of the mex file
  strcpy(progname, mexFunctionName());

  mexAtExit(CloseStreams);

  // call the computational routine
  ComputeCurve(0, NULL);

  if (TestRun)
    plhs[0] = mxCreateString("");
  else
    plhs[0] = mxCreateString(runname);

  return;
}


/*
 *====================================================================================================================================
 * R interface function
 *====================================================================================================================================
 */

#elif defined(R_PACKAGE)

SEXP FindCurve(SEXP moduleName, SEXP bifType, SEXP initVals, SEXP stepsize, SEXP curveVals, SEXP parVals, SEXP optVals)

{
  int   i, j, ncols, tmpint;
  char  optname[MAX_STR_LEN], optval[MAX_STR_LEN], tmpstr[MAX_STR_LEN];
  SEXP  resfil;

  InitialiseVars();

  //============================== Process the options argument ======================================================================

  if (length(optVals) && (!isString(optVals)))
    error("\nOptions should either be specified as NULL or as a vector of strings\n\n");
  strcpy(optstring, "");

  ncols =length(optVals);

  strcpy(optstring, "");
  for (i = 0; i < ncols; i++)
    {
      strcpy(optname, CHAR(STRING_ELT(optVals, i)));
      if (!((!strcmp(optname, "par1")) || (!strcmp(optname, "par2")) || (!strcmp(optname, "EXTfun")) || (!strcmp(optname, "EXTpar")) ||
            (!strcmp(optname, "negative")) || (!strcmp(optname, "noBP")) || (!strcmp(optname, "noEXT")) || (!strcmp(optname, "noLP")) ||
            (!strcmp(optname, "single")) || (!strcmp(optname, "silent")) || (!strcmp(optname, "test")) || (!strcmp(optname, "report"))))
        error("\nIllegal option %s!\n\n", optname);

      if (!strcmp(optname, "negative") || !strcmp(optname, "noBP") || !strcmp(optname, "noEXT") || !strcmp(optname, "noLP") ||
          !strcmp(optname, "single") || !strcmp(optname, "silent") || !strcmp(optname, "test"))
        {
          if (!strcmp(optname, "negative")) AllowNegative = 1;
          if (!strcmp(optname, "noBP")) BPdetection = 0;
          if (!strcmp(optname, "noEXT")) EXTdetection = 0;
          if (!strcmp(optname, "noLP")) LPdetection = 0;
          if (!strcmp(optname, "single")) DoSingle  = 1;
          if (!strcmp(optname, "silent")) ConsoleSilent  = 1;
          if (!strcmp(optname, "test")) TestRun     = 1;

          // optstring still equal to ""
          if (!strlen(optstring))
            strcat(optstring, "c('");
          else
            strcat(optstring, ", '");
          strcat(optstring, optname);
          strcat(optstring, "'");
          continue;
        }

      if (!(++i < ncols)) error("\nNo value specified for option %s!\n\n", optname);

      strcpy(optval, CHAR(STRING_ELT(optVals, i)));

      // optstring still equal to ""
      if (!strlen(optstring))
        strcat(optstring, "c('");
      else
        strcat(optstring, ", '");
      strcat(optstring, optname);
      strcat(optstring, "', '");
      strcat(optstring, optval);
      strcat(optstring, "'");

      if (!strcmp(optname, "par1"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            error("\nIndex of first parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          Bifparone = tmpint;
        }
      else if (!strcmp(optname, "par2"))
        {
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            error("\nIndex of second parameter (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          Bifpartwo = tmpint;
        }
      else if (!strcmp(optname, "EXTfun"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          error("\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= EXTRAOUTPUT_DIM))
            error("\nIndex of the element in the ExtraOutput[] vector for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n",
                  tmpint, EXTRAOUTPUT_DIM);
          EXTfun = tmpint;
        }
      else if (!strcmp(optname, "EXTpar"))
        {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
          error("\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
          tmpint = atoi(optval);
          if ((tmpint < 0) || (tmpint >= PARAMETER_NR))
            error("\nIndex of the parameter for EXT continuation (%d) not in the appropriate range (0 <= i < %d)!\n", tmpint, PARAMETER_NR);
          EXTpar = tmpint;
        }
      else if (!strcmp(optname, "report"))
        {
          tmpint = atoi(optval);
          ReportLevel = max(tmpint, 1);
        }
    }

  if (strlen(optstring))
    strcat(optstring, ")");
  else
    strcpy(optstring, "NULL");

  //============================== Process the bifurcation type argument =============================================================

  if ((!isString(bifType)) || (length(moduleName) != 1))
    error("\nBifurcation type must be a single string (BP, EQ, LP or EXT).\n\n");

  strcpy(ContinuationString, CHAR(STRING_ELT(bifType, 0)));
  if (!((!strcmp(ContinuationString, "BP")) || (!strcmp(ContinuationString, "EQ")) || (!strcmp(ContinuationString, "LP")) ||
        (!strcmp(ContinuationString, "EXT"))))
    error("\nBifurcation type must be a single string (BP, EQ, LP or EXT).\n\n");

  if (!strcmp(ContinuationString, "BP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = BP;
    }
  else if (!strcmp(ContinuationString, "EQ"))
    {
      pntdim    = EQUATIONS_DIM + 1;
      CurveType = EQ;
    }
  else if (!strcmp(ContinuationString, "LP"))
    {
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = LP;
    }
  else if (!strcmp(ContinuationString, "EXT"))
    {
#if !defined(EXTRAOUTPUT_DIM) || (EXTRAOUTPUT_DIM < 1)
      error("\nEXT continuation is only possible if EXTRAOUTPUT_DIM is assigned a positive value (now %d)!\n", EXTRAOUTPUT_DIM);
#endif
      pntdim    = EQUATIONS_DIM + 2;
      CurveType = EXT;
    }

  if (Bifpartwo == Bifparone)
    {
      if (((CurveType == BP) || (CurveType == LP) || (CurveType == EXT)))
        error("\nIndex of first and second bifurcation parameter the same!\nTwo parameter continuation not possible!\n");

      // Here we only end up in case of EQ continuation
      // Bifparone reset on command-line to Bifpartwo
      if ((Bifparone != DefaultBifparone) && (Bifpartwo != DefaultBifparone))
        Bifpartwo = DefaultBifparone;
      else if ((Bifparone != DefaultBifpartwo) && (Bifpartwo != DefaultBifpartwo))
        Bifpartwo = DefaultBifpartwo;
      else
        {
          for (i = 0; i < PARAMETER_NR; i++)
            if (i != Bifparone) break;
          Bifpartwo = (int)i;
        }
      warning("\nBifurcation parameter #2 reset to parameter #%d for locating branching points\n", Bifpartwo);
    }

  //================================ Process the parameters argument =================================================================

  ncols = length(parVals);
  if (isReal(parVals) && (ncols == PARAMETER_NR))
    memcpy(parameter, REAL(parVals), ncols*sizeof(double));
  else if (ncols)
    warning("\nParameter argument ignored as it is not a row vector of length PARAMETER_NR\n\n");

  if (ncols)
    {
      strcpy(parstring, "c(");
      for (i = 0; i < ncols; i++)
        {
          if (i) strcat(parstring, ", ");
          sprintf(tmpstr, "%.6G", REAL(parVals)[i]);
          strcat(parstring, tmpstr);
        }
      strcat(parstring, ")");
    }
  else
    strcpy(parstring, "NULL");

  //============================== Process the curve parameters argument =============================================================

  ncols = length(curveVals);
  if (!isReal(curveVals) || (!length(curveVals)))
    error("\nThe curve settings argument should be a vector of real values\n\n");

  if (CurveType == EQ)
    {
      if (ncols == 2)
        {
          pntmin[0] = REAL(curveVals)[0];
          pntmax[0] = REAL(curveVals)[1];
        }
      else if (ncols == 2*pntdim)
        {
          for (i = 0, j = 0; i < ncols; i++, j++)
            {
              pntmin[j] = REAL(curveVals)[i];
              i++;
              pntmax[j] = REAL(curveVals)[i];
            }
        }
      else
        error("\nFor %s continuation bounds argument must be a row vector of length 2 or length %d.\n", ContinuationString, 2*pntdim);
    }
  else
    {
      if (ncols == 4)
        {
          pntmin[0] = REAL(curveVals)[0];
          pntmax[0] = REAL(curveVals)[1];

          pntmin[pntdim - 1] = REAL(curveVals)[2];
          pntmax[pntdim - 1] = REAL(curveVals)[3];
        }
      else if (ncols == 2*pntdim)
        {
          for (i = 0, j = 0; i < ncols; i++, j++)
            {
              pntmin[j] = REAL(curveVals)[i];
              i++;
              pntmax[j] = REAL(curveVals)[i];
            }
        }
      else
        error("\nFor %s continuation bounds argument must be a row vector of length 4 or length %d.\n", ContinuationString, 2*pntdim);
    }

  strcpy(curvestring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(curvestring, ", ");
      sprintf(tmpstr, "%.6G", REAL(curveVals)[i]);
      strcat(curvestring, tmpstr);
    }
  strcat(curvestring, ")");

  //================================ Process the step size argument ==================================================================

  if (isReal(stepsize) && (length(stepsize) == 1))
    {
      Maxcurvestep = REAL(stepsize)[0];
      curvestep    = Maxcurvestep;
    }
  else
    error("\nStep size argument should be a single, real value\n\n");

  //============================== Process the initial point argument ================================================================

  memset((void *)point, 0, MAX_PNTDIM*sizeof(double));
  ncols = length(initVals);
  if (!isReal(initVals) || (ncols != pntdim))
    error("\nInitial point for bifurcation must be a row vector of length %d\n", pntdim);

  memcpy(point, REAL(initVals), ncols*sizeof(double));
  parameter[Bifparone]                      = point[0];
  if (CurveType != EQ) parameter[Bifpartwo] = point[pntdim - 1];

  strcpy(pntstring, "c(");
  for (i = 0; i < ncols; i++)
    {
      if (i) strcat(pntstring, ", ");
      sprintf(tmpstr, "%.6G", point[i]);
      strcat(pntstring, tmpstr);
    }
  strcat(pntstring, ")");

  //============================= Get the program name ===============================================================================
  // Get the name of the module file
  if (isString(moduleName) && (length(moduleName) == 1))
    strcpy(progname, CHAR(STRING_ELT(moduleName, 0)));
  else
    error("\nModel name argument must be a single string\n\n");

  // call the computational routine
  ComputeCurve(0, NULL);

  fflush(NULL);
  if (biffile) fclose(biffile);
  if (errfile) fclose(errfile);
  if (outfile) fclose(outfile);

  ResetCurve();
  if (odesolveMem) free(odesolveMem);
  odesolveMem = NULL;

  if (TestRun)
    PROTECT(resfil = mkString(""));
  else
    PROTECT(resfil = mkString(runname));

  UNPROTECT(1);

  return resfil;
}

/*==================================================================================================================================*/
#endif
