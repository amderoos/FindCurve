/***
   NAME
     curve
   DESCRIPTION
     This module implements routines that are specifically used in locating
     points on the equilibrium branch of a structured population model.

    Copyright (C) 2015, Andre M. de Roos, University of Amsterdam

    This file is part of the FindCurve software package.

    This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this software. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Jan 24, 2018
***/
#ifndef CURVE
#define CURVE
#endif
#include "globals.h"


/*
 *====================================================================================================================================
 *  Some numerical settings
 *====================================================================================================================================
 */

#ifndef RHSMAX
#define RHSMAX                    0.99
#endif
#ifndef CENTRALDIFF
#define CENTRALDIFF               1                                                 // Compute derivative by central difference
#endif
#ifndef FUNCTOL
#define FUNCTOL                   1.0E-8
#endif

#define STEP_DOUBLE               4
#define STEP_HALF                 (MAXITER/2)

#define FEMTO                     1.0E-15

#define SAFETY                    0.8
#define MIN_SCALE                 1.0E-6
#define MAX_SCALE                 1.0E6

#ifdef _MSC_VER
#include <float.h>
#define issane(a)                 ((_fpclass(a) == _FPCLASS_NN) || (_fpclass(a) == _FPCLASS_NZ) || (_fpclass(a) == _FPCLASS_PZ) || (_fpclass(a) == _FPCLASS_PN))
#else
#define issane(a)                 ((fpclassify(a) == FP_ZERO) || (fpclassify(a) == FP_NORMAL))
#endif

static int                        *oldscale;
static int                        FirstTangent = 1;
static int                        fast_iters = 0, slow_iters = 0;
static int                        LPImmediateReturn  = 0;
static int                        BPImmediateReturn  = 0;
static int                        EXTImmediateReturn = 0;
const double                      upper = ((2.0 - SAFETY)*10.0), lower = SAFETY;

static double                     *dDETBaseMem    = NULL;
static LAPACK_SIZE_T              *iDETBaseMem    = NULL;
static int                        dDETBaseMemDim = 0, iDETBaseMemDim = 0;

static double                     *dEVBaseMem    = NULL;
static LAPACK_SIZE_T              *iEVBaseMem    = NULL;
static int                        dEVBaseMemDim = 0, iEVBaseMemDim = 0;

static double                     *dSLSBaseMem    = NULL;
static LAPACK_SIZE_T              *iSLSBaseMem    = NULL;
static int                        dSLSBaseMemDim = 0, iSLSBaseMemDim = 0;


/*==================================================================================================================================*/

double anorm(int rows, int cols, double *a)
{
  register int  i;
  double        tmp, maxval = 0.0;

  for (i = 0; i < rows; i++)
    {
      tmp    = ASUM(cols, a + i*cols, 1);
      maxval = max(maxval, tmp);
    }

  return maxval;
}


/*==================================================================================================================================*/

int SetScales(double *point, int pntdim)

/*
 * This routine scales the vector of problem variables to within reasonable
 * bounds. This makes them more comparable, which is advantageous for the
 * computations.
 */

{
  register int  i, scaleset = 0, newscale;
  double        tmp, val;

  if (!oldscale)
    {
      oldscale = calloc(pntdim, sizeof(int));
      for (i = 0; i < pntdim; i++) oldscale[i] = 0;
    }
  if (!point) return scaleset;

  for (i = 0; i < pntdim; i++)
    {
      val = fabs(point[i]);
      if ((val < lower) || (val > upper))
        {
          tmp      = max(val*pnt_scale[i], MIN_SCALE);
          tmp      = min(tmp, MAX_SCALE);
          tmp      = floor(log10(tmp) + FUNCTOL);
          newscale = (int)tmp;

          if (newscale != oldscale[i])
            {
              tmp = pow(10.0, tmp);
              point[i] *= pnt_scale[i]/tmp;
              pnt_scale[i] = tmp;
              oldscale[i]  = newscale;
              scaleset     = i + 1;
            }
        }
    }

  return scaleset;
}


/*==================================================================================================================================*/

int FindPoint(const int pntdim, double *guess, double *JacImport, double *tanvec, double ytol, double rhstol, const int max_iter,
              int (*fnc)(double *, double *))

/*
 * FindPoint -  Routine locates a point on a curve determined by a
 *              system of non-linear, algebraic equations.
 *              The iteration adjusts the vector-elements following a simple
 *              Newton-Chord method with Broyden update (see Kuznetsov pg. 418).
 *              Pseudo-arclength continuation is used to continue past curve folds.
 *
 * Arguments -  pntdim    : The dimension of the solution point on the curve.
 *                          The dimension of the system of equations is
 *                          assumed to be exactly 1 less.
 *              guess     : Pointer to an array containing the initial point
 *                          to start the iteration from. The first element of
 *                          the vector is assumed to be non-adjustable parameter.
 *              ytol      : Tolerance determining when change in y equals zero.
 *              rhstol    : Tolerance determining when RHS equals zero.
 *              max_iter  : Maximum stepnumber allowed in iteration.
 *              fnc       : Pointer to function specifying the system of
 *                          equations. The function must have a (double)
 *                          pointer as first argument, containing the point
 *                          in which to evaluate the system and a (double)
 *                          pointer as second argument, containing the
 *                          results after evaluation of the equations.
 */

{
  register int  iter, i, j;
  int           rhsdim;
  int           pntdim2 = pntdim*pntdim;
  int           retcode = NO_CONVERGENCE;
  double        ynorm, dynorm, rhsnorm;
  double        *basemem, *y, *tv, *dy, *rhs;
  double        *Jac, *JacCopy;

  y = basemem = calloc((4*pntdim + 2*pntdim2), sizeof(double));
  if (!basemem) return ReportMemError("FindPoint");

  tv      = y + pntdim;
  dy      = tv + pntdim;
  rhs     = dy + pntdim;
  Jac     = rhs + pntdim;
  JacCopy = Jac + pntdim2;

  // If tangent vector is given, we are doing pseudo-arc length continuation
  rhsdim = pntdim - (tanvec != NULL);

  ReportMsg("\nLocating new solution point\n");
  COPY(pntdim, guess, 1, y, 1);
  memset((void *)dy, 0, pntdim*sizeof(double));

  // The iteration loop
  for (iter = 0; iter < max_iter; iter++)
    {
      // Compute norm of Y and of dY
      ynorm  = NRM2(pntdim, y, 1);
      dynorm = NRM2(pntdim, dy, 1);
      if (!issane(ynorm) || !issane(dynorm))
        {
          ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }

      // Compute RHS and its norm
      memset((void *)rhs, 0, pntdim*sizeof(double));
      if ((*fnc)(y, rhs) == FAILURE)
        {
          ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
          retcode = FAILED_EVALUATION;
          break;
        }
      // The dimension of rhs is pntdim, for example in case of BP localisation
      rhsnorm = NRM2(pntdim, rhs, 1);
      if (!iter) ReportMsg("Start:");
      for (i = 0; i < pntdim; i++) ReportMsg("\t%12.5E", y[i]*pnt_scale[i]);
      ReportMsg("\tdY  norm: %12.5E\tRHS norm: %12.5E\n", dynorm, rhsnorm);

      // Return if converged or diverged
      if ((!issane(rhsnorm)) || (rhsnorm/(1.0 + rhsnorm) > RHSMAX))
        {
          ErrorMsg(__FILE__, __LINE__, "Norm overflow in FindPoint");
          retcode = NORM_OVERFLOW;
          break;
        }
      // The dimension of rhs is pntdim, for example in case of BP localisation
      else if ((rhsnorm < pntdim*rhstol) && (dynorm < pntdim*ytol))
        {
          COPY(pntdim, y, 1, guess, 1);
          if ((Stepchange) && (Stepreduce == 1))
            {
              if (iter < STEP_DOUBLE)
                {
                  fast_iters++;
                  slow_iters = 0;
                  if (fast_iters == 2)
                    {
                      curvestep *= 1.5;
                      fast_iters = 0;
                    }
                }
              else if (iter > STEP_HALF)
                {
                  slow_iters++;
                  fast_iters = 0;
                  if (slow_iters == 2)
                    {
                      curvestep *= 0.5;
                      slow_iters = 0;
                    }
                }
              else
                {
                  fast_iters = 0;
                  slow_iters = 0;
                }
            }

          retcode = SUCCES;
          break;
        }

      // Compute Jacobian every Jacobian_Updates steps, otherwise the Jacobian is updated
      // via a Broyden update (see below)
      if ((iter == 0) && JacImport)
        memcpy(Jac, JacImport, pntdim*rhsdim*sizeof(double));
      else if (!(iter % Jacobian_Updates))
        {
          ReportMsg("%-s", "Computing jacobian");
          memset((void *)Jac, 0, (pntdim*pntdim)*sizeof(double));
          /*
           * Notice that the Jacobian is stored as
           *
           *          |dF1/dy1 ... dFn/dy1|
           *          |dF1/dy2 ... dFn/dy2|
           *      J = |   .           .   |
           *          |   .           .   |
           *          |dF1/dyn ... dFn/dyn|
           *
           * The matrix is hence stored in column-wise (fortran) style.
           * From a C perspective this means that all coefficients pertaining to yi are to be found
           * in ROW i (as opposed to column i).
           *
           * Solving J.dy = -F(y) with dy = (dy1 ... dyn) and F(y) = (F1(y) ... Fn(y)) requires the variable
           * trans[1] to be defined as {"N"} (see programs/various/testlapack.c for details).
           */
          Jacobian(pntdim, y, rhsdim, Jac, fnc, FORWARD);                           // Compute J = F_x(X^k)
          ReportMsg(".....Ok!\n");
        }
      else                                                                          // Broyden update of Jacobian
        {
          dynorm = DOT(pntdim, dy, 1, dy, 1);
          // See 10.7 on pg. 419 in Kuznetsov. Notice though that Jac is the transposed jacobian
          for (i = 0; i < pntdim; i++)
            for (j = 0; j < rhsdim; j++) Jac[j + i*rhsdim] += rhs[j]*dy[i]/dynorm;
        }

      // Extend the Jacobian matrix to include an additional row for the tangent
      // vector if we are doing pseudo-archlength continuation
      // If tangent is not given rhsdim == pntdim and we follow simple Newton
      memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
      for (i = 0; i < pntdim; i++)
        COPY(rhsdim, Jac + i*rhsdim, 1, JacCopy + i*pntdim, 1);                     // Extract dF/dx
      memset((void *)dy, 0, pntdim*sizeof(double));
      AXPY(rhsdim, -1.0, rhs, 1, dy, 1);

      if (tanvec != NULL)
        {
          // When tangent is present, find new point via pseudo-arclength continuation
          for (i = 0; i < pntdim; i++) *(JacCopy + i*pntdim + rhsdim) = tanvec[i];
          COPY(pntdim, y, 1, tv, 1);
          AXPY(pntdim, -1.0, guess, 1, tv, 1);
          dy[rhsdim] = DOT(pntdim, tv, 1, tanvec, 1);
        }

      // Solve the linear system
      retcode = SolveLinearSystem(pntdim, JacCopy, dy, ytol);
      if (retcode != SUCCES) break;

      // Adjust point
      AXPY(pntdim, 1.0, dy, 1, y, 1);
      retcode = NO_CONVERGENCE;
    }

  free(basemem);

  if (retcode == SUCCES)
    ReportMsg("New solution point found\n\n");
  else
    ReportMsg("Locating new solution point failed\n\n");

  return retcode;
}


/*==================================================================================================================================*/

int TangentVec(const int pntdim, double *sol, double *JacExport, double *tanvec, int (*fnc)(double *, double *), double *det, const double ytol)

/*
 * TangentVec - Routine determines the direction of the curve defined by the
 *              system of equations
 *
 *                  F(y) = 0
 *
 *              The point y is considered to have a dimension of exactly 1
 *              larger than the number of equations (i.e. the dimension of
 *              F(y)).
 *
 * Arguments -  pntdim  : The dimension of the solution point on the curve.
 *              y       : Pointer to an array containing the fixed point
 *              tanvec  : Pointer to return tangent vector
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 */

{
  register int  j;
  int           rhsdim = pntdim - 1, pntdim2 = pntdim*pntdim, retcode;
  double        norm;
  double        *basemem, *y, *Jac, *JacCopy;

  y = basemem = calloc((pntdim + 2*pntdim2), sizeof(double));
  if (!basemem) return ReportMemError("TangentVec");

  Jac     = y + pntdim;
  JacCopy = Jac + pntdim2;

  // Initialize
  COPY(pntdim, sol, 1, y, 1);
  norm = NRM2(pntdim, y, 1);
  if (!issane(norm))
    {
      ErrorMsg(__FILE__, __LINE__, "Norm overflow in curvedir");
      free(basemem);
      return NORM_OVERFLOW;
    }

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence).
  ReportMsg("\nComputing curve direction     ");
  Jacobian(pntdim, y, rhsdim, JacCopy, fnc, CENTRALDIFF);
  if (JacExport) memcpy(JacExport, JacCopy, pntdim*rhsdim*sizeof(double));
  ReportMsg(".....Ok!\n");

  // Append the current tangent vector as the last row to the jacobian to
  // preserve direction. See the matcont manual at
  // http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  // Notice, however, it is here added as the last COLUMN because of the
  // Fortran column-wise storage!

  for (j = 0; j < pntdim; j++)
    {
      COPY(rhsdim, JacCopy + j*rhsdim, 1, Jac + j*pntdim, 1);                       // Extract dF/dy
      *(Jac + j*pntdim + rhsdim) = tanvec[j];
    }

  memset((void *)JacCopy, 0, (pntdim*pntdim)*sizeof(double));
  COPY(pntdim2, Jac, 1, JacCopy, 1);
  memset((void *)tanvec, 0, pntdim*sizeof(double));
  tanvec[rhsdim] = 1.0;
  Stepchange     = 0;

  // Solve the linear system
  retcode = SolveLinearSystem(pntdim, JacCopy, tanvec, ytol);
  if (retcode != SUCCES)
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to solve for tangent vector in TangentVec()");
      memset((void *)tanvec, 0, pntdim*sizeof(double));
      tanvec[0] = 1.0;
      free(basemem);
      return retcode;
    }

  if (det)
    {
      // Replace the last row of the (saved) Jacobian with the newly computed tangent vector
      // to compute the determinant for BP detection
      for (j = 0; j < pntdim; j++)
        {
          COPY(rhsdim, Jac + j*pntdim, 1, JacCopy + j*pntdim, 1);
          *(JacCopy + j*pntdim + rhsdim) = tanvec[j];
        }
      Determinant(pntdim, JacCopy, det, NULL);
    }
  norm = NRM2(pntdim, tanvec, 1); /* Normalize and store      */
  SCAL(pntdim, 1.0/norm, tanvec, 1);

  if (FirstTangent && (tanvec[0] < 0.0)) SCAL(pntdim, -1.0, tanvec, 1);
  FirstTangent = 0;

  free(basemem);

  return SUCCES;
}


/*==================================================================================================================================*/

int Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac, int (*fnc)(double *, double *), int method)
/*
 * Routine determines the Jacobian of the n-dimensional function F(y) w.r.t. the m-dimensional
 * variable y at the current point given by 'pnt'. The routine hence returns in 'jac' the
 * following matrix of partial derivatives:
 *
 *            |dF1/dy1 ... dFn/dy1|
 *            |   .           .   |
 *      Df =  |   .           .   |
 *            |   .           .   |
 *            |dF1/dym ... dFn/dym|
 *
 * Notice that all coefficients pertaining to yi are to be found in ROW i (as opposed to column i).
 * The matrix is hence stored in column-wise (fortran) style.
 */

{
  register int  j;
  double        *basemem, *y, *rhs;
  double        ydif, yfac, old;

  y = basemem = calloc(2*pntdim, sizeof(double));
  if (!basemem) return ReportMemError("Jacobian");

  rhs = y + pntdim;

  // Initialize
  COPY(pntdim, pnt, 1, y, 1);
  if (method == FORWARD)
    {
      if ((*fnc)(y, rhs) == FAILURE)
        {
          ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
          free(basemem);
          return FAILURE;
        }
    }

  memset((void *)jac, 0, (pntdim*fncdim)*sizeof(double));
  for (j = 0; j < pntdim; j++)
    {
      old  = y[j];
      ydif = max(fabs(Jacobian_Step*old), Jacobian_Min_Step/pnt_scale[j]);
      y[j] = old + ydif;
      // Trick to reduce precision errors. See Num. Recipes 9.7, pg. 388
      ydif = y[j] - old;
      yfac = 1.0/((1 + (method == CENTRALDIFF))*ydif);

      if ((*fnc)(y, jac + j*fncdim) == FAILURE)
        {
          ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
          free(basemem);
          return FAILURE;
        }

      if (method == CENTRALDIFF)
        {
          y[j] = old - ydif;
          memset((void *)rhs, 0, fncdim*sizeof(double));
          if ((*fnc)(y, rhs) == FAILURE)
            {
              ErrorMsg(__FILE__, __LINE__, "Right-hand side computation failed");
              free(basemem);
              return FAILURE;
            }
        }
      AXPY(fncdim, -1.0, rhs, 1, jac + j*fncdim, 1);
      SCAL(fncdim, yfac, jac + j*fncdim, 1);

      y[j] = old;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt())
        {
          free(basemem);
          return FAILURE;
        }
#endif
    }

  free(basemem);
  return SUCCES;
}


/*==================================================================================================================================*/

int LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), const int method, const int lpcurve, double *curvedir, const double ytol)

/*
 * LPcondition -  Routine computes the factor determining the location of a limit point, i.e.
 *                the parameter component of the tangent vector. This component always has
 *                index 0 in the vector of the solution point.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the parameters) plus the dimension of the vector of
 *                        state variables in case of lpcurve = 1, otherwise the dimension
 *                        of y equals 1 plus the vector of state variables.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables y. The last element is assumed to be the
 *                        second parameter in the LP continuation
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              lpcurve : Routine called for detection of LP in EQ curve (0) or
 *                        during computation of LP curve (1)
 */
{
  register int  j;
  const int     lppntdim = pntdim - lpcurve;
  int           rhsdim = lppntdim - 1, lppntdim2 = lppntdim*lppntdim, retcode;
  int           maxind;
  double        *basemem, *rhs, *Jac, *JacCopy, *tanvec;
  double        norm, maxcond, cond;

  // Prevent recurrence
  *curvedir = 0.0;
  if (LPImmediateReturn) return SUCCES;

  rhs = basemem = calloc((pntdim + pntdim*pntdim + lppntdim*lppntdim + lppntdim), sizeof(double));
  if (!basemem) return ReportMemError("LPcondition");

  LPImmediateReturn = 1;

  Jac     = rhs + pntdim;
  JacCopy = Jac + pntdim*pntdim;
  tanvec  = JacCopy + lppntdim*lppntdim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  // Notice that when continuing a LP curve (lpcurve = 1) we have to call Jacobian() with the
  // full dimension (pntdim) of y to pass the entire argument vector to Equation().
  // As a result, the Jacobian matrix will have one row more than we really need,
  // representing the derivatives w.r.t. to the 2nd parameter of the LP continuation.
  // This additional row, with index lppntdim = pntdim-1 will be ignored.
  // Also notice that the last column will be unassigned, as the systems of equations returned
  // only has size rhsdim = pntdim-lpcurve-1. This last column is hence explicitly set to 0
  Jacobian(pntdim, y, lppntdim, Jac, fnc, method);
  for (j = 0; j < lppntdim; j++) Jac[j*lppntdim + rhsdim] = 0.0;

  /*
    * When lpcurve = 1 the Jacobian equals the following (n+2)x(n+1) matrix of partial
    * derivatives:
    *
    *           |dF1/dp1 ... dFn/dp1  0|
    *           |dF1/dy1 ... dFn/dy1  0|
    *           |   .           .     0|
    *      Df = |   .           .     0|
    *           |   .           .     0|
    *           |dF1/dyn ... dFn/dyn  0|
    *           |dF1/dp2 ... dFn/dp2  0|
    *
    * Otherwise, when lpcurve = 0 the (n+1)x(n+1) matrix
    *
    *           |dF1/dp1 ... dFn/dp1  0|
    *           |dF1/dy1 ... dFn/dy1  0|
    *           |   .           .     0|
    *      Df = |   .           .     0|
    *           |   .           .     0|
    *           |dF1/dyn ... dFn/dyn  0|
    *
    * In which n = pntdim-lpcurve-1 (i.e. equal to the number of state variables).
    * Notice that all coefficients pertaining to yi are to be found in ROW i (as 
    * opposed to column i). The matrix is hence stored in column-wise (fortran) style.
    */

  // Additional call to reset global variables
  (*fnc)(y, rhs);

  // Find the most non-singular matrix (largest inverse condition)
  for (j = 0, maxcond = 0.0, maxind = -1; j < lppntdim; j++)
    {
      // LU decompose the matrix and compute determinant
      COPY(lppntdim2, Jac, 1, JacCopy, 1);
      JacCopy[j*lppntdim + rhsdim] = 1.0;

      if (Determinant(lppntdim, JacCopy, NULL, &cond) == SUCCES)
        {
          if (cond > maxcond + FEMTO)
            {
              maxcond = cond;
              maxind  = j;
            }
        }
    }

  if (maxind == -1)
    {
      ErrorMsg(__FILE__, __LINE__, "No non-singular matrix found in LPcondition()");
      free(basemem);
      LPImmediateReturn = 0;
      return FAILURE;
    }

  Jac[maxind*lppntdim + rhsdim] = 1.0;
  COPY(lppntdim2, Jac, 1, JacCopy, 1);
  memset((void *)tanvec, 0, lppntdim*sizeof(double));
  tanvec[rhsdim] = 1.0;

  // Solve the linear system
  retcode = SolveLinearSystem(lppntdim, JacCopy, tanvec, ytol);
  if (retcode != SUCCES)
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to solve for tangent vector in LPcondition()");
      free(basemem);
      LPImmediateReturn = 0;
      return retcode;
    }

  norm = NRM2(lppntdim, tanvec, 1);                                                 // Normalize and store
  SCAL(lppntdim, 1.0/norm, tanvec, 1);

  *curvedir = tanvec[0];
  free(basemem);
  LPImmediateReturn = 0;

  return SUCCES;
}


/*==================================================================================================================================*/

int Determinant(const int N, double *M, double *det, double *cond)

{
  char          whichnorm;
  int           j;
  int           retval = SUCCES, memneeded;
  double        *A, *work, norm;
  LAPACK_SIZE_T nc = N, lwork = 4*N, *ipiv, *iwork, liwork = N, info;

  // Allocate temporarily minimally allowed size for workspace arrays
  memneeded = N*N + lwork;
  if (memneeded > dDETBaseMemDim)
    {
      dDETBaseMem = realloc(dDETBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dDETBaseMem == NULL)
        {
          if (dDETBaseMem) free(dDETBaseMem);
          dDETBaseMem    = NULL;
          dDETBaseMemDim = 0;
          if (iDETBaseMem) free(iDETBaseMem);
          iDETBaseMem    = NULL;
          iDETBaseMemDim = 0;
          return ReportMemError("Determinant");
        }
      dDETBaseMemDim = memneeded;
    }

  memneeded = N + liwork;
  if (memneeded > iDETBaseMemDim)
    {
      iDETBaseMem = realloc(iDETBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iDETBaseMem == NULL)
        {
          if (dDETBaseMem) free(dDETBaseMem);
          dDETBaseMem    = NULL;
          dDETBaseMemDim = 0;
          if (iDETBaseMem) free(iDETBaseMem);
          iDETBaseMem    = NULL;
          iDETBaseMemDim = 0;
          return ReportMemError("Determinant");
        }
      iDETBaseMemDim = memneeded;
    }

  A    = dDETBaseMem;
  work = A + N*N;

  // Copy the matrix
  COPY(N*N, M, 1, A, 1);

  memset((void *)iDETBaseMem, 0, iDETBaseMemDim*sizeof(LAPACK_SIZE_T));
  ipiv  = iDETBaseMem;
  iwork = ipiv + N;

  dgetrf(&nc, &nc, A, &nc, ipiv, &info);
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgetrf", abs((int)info));
      return ILLEGAL_INPUT;
    }

  if (det)
    {
      *det = 1.0;
      if (!info)
        {
          for (j = 0; j < N; j++)
            {
              if (ipiv[j] != (j + 1))
                *det *= -A[j*N + j];
              else
                *det *= A[j*N + j];
            }
        }
    }

  if (info > 0) return SINGULARITY;

  if (cond)
    {
      norm      = anorm(N, N, M);
      whichnorm = '1';
      dgecon(&whichnorm, &nc, A, &nc, &norm, cond, work, iwork, &info);
      if (info < 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in DGECON", abs((int)info));
          return ILLEGAL_INPUT;
        }
    }

  return retval;
}


/*==================================================================================================================================*/

int Eigenval(const int N, double *X, const int symmetric, double *eigval, const int eigenvaltype, double *righteigvec, double *lefteigvec,
             const double tol)

{
  /*
   * This function calculates the dominant, real eigenvalue of the N*N general
   * matrix X and the corresponding right and left eigenvector.
   *
   * The matrix has to be in normal C-wise row format and is transferred in this
   * routine to Format vector format.
   *
   * The eigenvalue is returned in *eigval, whereas the vector eigevec (length N)
   * contains the eigenvector. The content of X is not changed.
   *
   * This function first queries the Lapack routines for optimal workspace sizes.
   * These memoryblocks are then allocated and the decomposition is calculated
   * using the Lapack function "dgeevx". The allocated memory is preserved for
   * subsequent calls.
   */

  // balanc = 'B'           : Do scaling and permutation to make computation more accurate
  // jobvl  = jobvr  = 'V'  : Both eigenvectors are needed to compute error estimates of eigenvalues
  // sense  = 'E'           : Compute error estimates on eigenvalues only
  char          balanc = 'B', jobvl = 'V', jobvr = 'V', sense = 'E';
  char          jobz = 'V', range = 'I', uplo = 'U';
  int           i, j, dbasemem, ibasemem, memneeded;
  double        *Xc, *wr, *wi, *vr, *vl, *scale, *rconde, *rcondv, *work;
  double        abnrm, abstol, ddummy;
  int           retval = SUCCES, eigvalindx;
  LAPACK_SIZE_T nc     = N, ilo, ihi, lwork, info, *iwork, liwork, nfound, *isuppz;

  // Allocate temporarily minimally allowed size for workspace arrays
  if (symmetric)
    {
      dbasemem = 2*N*N + N;
      lwork    = 26*N;
    }
  else
    {
      dbasemem = 3*N*N + 5*N;
      lwork    = N*(N + 6);
    }

  memneeded = dbasemem + lwork;
  if (memneeded > dEVBaseMemDim)
    {
      dEVBaseMem = realloc(dEVBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      dEVBaseMemDim = memneeded;
    }

  if (symmetric)
    {
      ibasemem = 2*N;
      liwork   = 10*N;
    }
  else
    {
      ibasemem = 0;
      liwork   = (2*N - 2);
    }

  memneeded = ibasemem + liwork;
  if (memneeded > iEVBaseMemDim)
    {
      iEVBaseMem = realloc(iEVBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      iEVBaseMemDim = memneeded;
    }

  memset((void *)dEVBaseMem, 0, dEVBaseMemDim*sizeof(double));
  Xc = dEVBaseMem;
  vr = Xc + N*N;
  wr = vr + N*N;
  if (symmetric)
    work = wr + N;
  else
    {
      vl     = wr + N;
      wi     = vl + N*N;
      scale  = wi + N;
      rconde = scale + N;
      rcondv = rconde + N;
      work   = rcondv + N;
    }

  memset((void *)iEVBaseMem, 0, iEVBaseMemDim*sizeof(LAPACK_SIZE_T));
  iwork  = iEVBaseMem;
  isuppz = iwork + liwork;                                                          // Only used in case of symmetric matrix

  // Get the machine precisions
  abstol = dlamch("Safe minimum");

  // Rewrite C-style row vector format to Fortran-style column vector format
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) Xc[j*N + i] = X[i*N + j];

  ilo = ihi = N;
  // Query the Lapack routine for optimal sizes for workspace arrays
  if (symmetric)
    {
      lwork = liwork = -1;
      dsyevr(&jobz, &range, &uplo, &nc, Xc, &nc, &ddummy, &ddummy, &ilo, &ihi, &abstol, &nfound, wr, vr, &nc, isuppz, work, &lwork, iwork, &liwork, &info);
      lwork  = (int)work[0];
      liwork = (int)iwork[0];
    }
  else
    {
      lwork = -1;
      dgeevx(&balanc, &jobvl, &jobvr, &sense, &nc, Xc, &nc, wr, wi, vl, &nc, vr, &nc, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);
      lwork = (int)work[0];
    }

  // Free previous allocation and reallocate preferable workspaces, Check result
  memneeded = dbasemem + lwork;
  if (memneeded > dEVBaseMemDim)
    {
      dEVBaseMem = realloc(dEVBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      dEVBaseMemDim = memneeded;
    }

  memneeded = ibasemem + liwork;
  if (memneeded > iEVBaseMemDim)
    {
      iEVBaseMem = realloc(iEVBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iEVBaseMem == NULL)
        {
          if (dEVBaseMem) free(dEVBaseMem);
          dEVBaseMem    = NULL;
          dEVBaseMemDim = 0;
          if (iEVBaseMem) free(iEVBaseMem);
          iEVBaseMem    = NULL;
          iEVBaseMemDim = 0;
          return ReportMemError("Eigenval");
        }
      iEVBaseMemDim = memneeded;
    }

  memset((void *)dEVBaseMem, 0, dEVBaseMemDim*sizeof(double));
  Xc = dEVBaseMem;
  vr = Xc + N*N;
  wr = vr + N*N;
  if (symmetric)
    work = wr + N;
  else
    {
      vl     = wr + N;
      wi     = vl + N*N;
      scale  = wi + N;
      rconde = scale + N;
      rcondv = rconde + N;
      work   = rcondv + N;
    }

  memset((void *)iEVBaseMem, 0, iEVBaseMemDim*sizeof(LAPACK_SIZE_T));
  iwork  = iEVBaseMem;
  isuppz = iwork + liwork;                                                          // Only used in case of symmetric matrix

  // Rewrite C-style row vector format to Fortran-style column vector format
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) Xc[j*N + i] = X[i*N + j];

  // Now calculate the eigenvalues and vectors using optimal workspaces
  if (symmetric)                                                                    // Symmetric matrices
    {
      dsyevr(&jobz, &range, &uplo, &nc, Xc, &nc, &ddummy, &ddummy, &ilo, &ihi, &abstol, &nfound, wr, vr, &nc, isuppz, work, &lwork, iwork, &liwork, &info);

      // Check for convergence
      if (info < 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dsyevr()", abs((int)info));
          retval = ILLEGAL_INPUT;
        }
      else if (info > 0)
        {
          ErrorMsg(__FILE__, __LINE__, "The algorithm failed to compute eigenvalues!");
          retval = NO_CONVERGENCE;
        }
      else                                                                          // General matrices
        {
          if (eigenvaltype == DOMINANT) eigvalindx = 0;
          else                                                                      // if (eigenvaltype == MINIMUMNORM)
            {
              // Find the real eigenvalue with the smallest norm
              eigvalindx = -1;
              for (i = 0; i < N; i++)
                if ((eigvalindx == -1) || (fabs(wr[i]) < fabs(wr[eigvalindx]))) eigvalindx = i;

              if (eigvalindx == -1)
                {
                  ErrorMsg(__FILE__, __LINE__, "Did not find any real eigenvalue with smallest norm in Eigenval!");
                  retval = FAILURE;
                }
            }
          if (retval == SUCCES)
            {
              *eigval = wr[eigvalindx];
              COPY(nc, vr + eigvalindx*nc, 1, righteigvec, 1);
            }
        }
    }
  else
    {
      dgeevx(&balanc, &jobvl, &jobvr, &sense, &nc, Xc, &nc, wr, wi, vl, &nc, vr, &nc, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work, &lwork, iwork, &info);

      // Check for convergence
      if (info < 0)
        {
          ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgeevx()", abs((int)info));
          retval = ILLEGAL_INPUT;
        }
      else if (info > 0)
        {
          ErrorMsg(__FILE__, __LINE__, "The algorithm failed to compute eigenvalues!");
          retval = NO_CONVERGENCE;
        }
      else
        {
          if (eigenvaltype == DOMINANT)
            {
              eigvalindx = -1;
              for (i = 0; i < N; i++)
                if ((eigvalindx == -1) || (wr[i] > wr[eigvalindx])) eigvalindx = i;
            }
          else                                                                      // if (eigenvaltype == MINIMUMNORM)
            {
              // Find the real eigenvalue with the smallest norm
              eigvalindx = -1;
              for (i = 0; i < N; i++)
                if (wi[i] == (double)0.0)
                  if ((eigvalindx == -1) || (fabs(wr[i]) < fabs(wr[eigvalindx]))) eigvalindx = i;

              if (eigvalindx == -1)
                {
                  ErrorMsg(__FILE__, __LINE__, "Did not find any real eigenvalue with smallest norm in Eigenval!");
                  retval = FAILURE;
                }
            }
          if ((epsMach*abnrm/rconde[eigvalindx]) > tol)
            {
              if (eigenvaltype == DOMINANT)
                ErrorMsg(__FILE__, __LINE__, "The estimated error bound in the largest eigenvalue (%G) exceeds the tolerance level %G!",
                         (epsMach*abnrm/rconde[eigvalindx]), tol);
              else if (eigenvaltype == MINIMUMNORM)
                ErrorMsg(__FILE__, __LINE__, "The estimated error bound in the minimum eigenvalue (%G) exceeds the tolerance level %G!",
                         (epsMach*abnrm/rconde[eigvalindx]), tol);
              retval = NO_CONVERGENCE;
            }
          if (retval == SUCCES)
            {
              *eigval = wr[eigvalindx];
              COPY(nc, vr + eigvalindx*nc, 1, righteigvec, 1);
              // Left and right eigenvectors are the same
              if (lefteigvec) COPY(nc, vl + eigvalindx*nc, 1, lefteigvec, 1);
            }
        }
    }

  return retval;
}


/*==================================================================================================================================*/

int SolveLinearSystem(const int N, double *A, double *B, double tol)

{
  /*
   * This function solves the linear equation system A*x = B, where A is a NxN
   * matrix and B a N-dimensional vector.
   *
   * The matrix A has to be in Fortran column-wise format
   *
   * The solution is returned in *B, whereas the content of A is not changed.
   */

  char          fact, trans, equed, errstr[MAX_STR_LEN];
  int           memneeded;
  double        *Ac, *Af, *r, *c, *Bc, *x, *work, rcond, ferr = 0, berr;
  int           retval = SUCCES;
  LAPACK_SIZE_T nc = N, nrhs = 1, lwork = 4*N, *ipiv, *iwork, info;

  // Allocate temporarily minimally allowed size for workspace arrays
  memneeded = 2*N*N + 4*N + lwork;
  if (memneeded > dSLSBaseMemDim)
    {
      dSLSBaseMem = realloc(dSLSBaseMem, memneeded*sizeof(double));
      // Check for NULL-pointers
      if (dSLSBaseMem == NULL)
        {
          if (dSLSBaseMem) free(dSLSBaseMem);
          dSLSBaseMem    = NULL;
          dSLSBaseMemDim = 0;
          if (iSLSBaseMem) free(iSLSBaseMem);
          iSLSBaseMem    = NULL;
          iSLSBaseMemDim = 0;
          return ReportMemError("SolveLinearSystem");
        }
      dSLSBaseMemDim = memneeded;
    }

  memneeded = 2*N;
  if (memneeded > iSLSBaseMemDim)
    {
      iSLSBaseMem = realloc(iSLSBaseMem, memneeded*sizeof(LAPACK_SIZE_T));
      // Check for NULL-pointers
      if (iSLSBaseMem == NULL)
        {
          if (dSLSBaseMem) free(dSLSBaseMem);
          dSLSBaseMem    = NULL;
          dSLSBaseMemDim = 0;
          if (iSLSBaseMem) free(iSLSBaseMem);
          iSLSBaseMem    = NULL;
          iSLSBaseMemDim = 0;
          return ReportMemError("SolveLinearSystem");
        }
      iSLSBaseMemDim = memneeded;
    }

  memset((void *)dSLSBaseMem, 0, dSLSBaseMemDim*sizeof(double));
  Ac   = dSLSBaseMem;
  Af   = Ac + N*N;
  r    = Af + N*N;
  c    = r + N;
  Bc   = c + N;
  x    = Bc + N;
  work = x + N;

  memset((void *)iSLSBaseMem, 0, iSLSBaseMemDim*sizeof(LAPACK_SIZE_T));
  ipiv  = iSLSBaseMem;
  iwork = ipiv + N;

  fact  = 'E';
  trans = 'N';

  // Fill the matrix and the right-hand side vector
  COPY(N*N, A, 1, Ac, 1);
  COPY(N, B, 1, Bc, 1);

  dgesvx(&fact, &trans, &nc, &nrhs, Ac, &nc, Af, &nc, ipiv, &equed, r, c, Bc, &nc, x, &nc, &rcond, &ferr, &berr, work, iwork, &info);

  // Check for singularity of the matrix
  if (info < 0)
    {
      ErrorMsg(__FILE__, __LINE__, "Illegal value for parameter %d in dgesvx()", abs((int)info));
      retval = ILLEGAL_INPUT;
    }
  else if (info > 0)
    {
      ErrorMsg(__FILE__, __LINE__, "(Nearly) Singular matrix in SolveLinearSystem()!");
      retval = SINGULARITY;
    }
  else
    {
      if (ferr > 10*tol)
        {
          // Matlab does not handle correctly the direct printing of the double number via ErrorMsg
          sprintf(errstr, "Warning: The estimated error bound in the solution of the linear system A*x = B (%G) exceeds the tolerance level %G by an order of magnitude!\n",
                  ferr, tol);
          ErrorMsg(__FILE__, __LINE__, errstr);
        }
      COPY(N, x, 1, B, 1);
    }

  return retval;
}


/*==================================================================================================================================*/

void ResetCurve(void)
{
  if (oldscale) free(oldscale);
  oldscale            = NULL;

  if (dDETBaseMem) free(dDETBaseMem);
  dDETBaseMem         = NULL;
  dDETBaseMemDim      = 0;
  if (iDETBaseMem) free(iDETBaseMem);
  iDETBaseMem         = NULL;
  iDETBaseMemDim      = 0;

  if (dEVBaseMem) free(dEVBaseMem);
  dEVBaseMem          = NULL;
  dEVBaseMemDim       = 0;
  if (iEVBaseMem) free(iEVBaseMem);
  iEVBaseMem          = NULL;
  iEVBaseMemDim       = 0;

  if (dSLSBaseMem) free(dSLSBaseMem);
  dSLSBaseMem         = NULL;
  dSLSBaseMemDim      = 0;
  if (iSLSBaseMem) free(iSLSBaseMem);
  iSLSBaseMem         = NULL;
  iSLSBaseMemDim      = 0;

  fast_iters          = 0;
  slow_iters          = 0;
  FirstTangent        = 1;

  LPImmediateReturn   = 0;
  BPImmediateReturn   = 0;
  EXTImmediateReturn  = 0;

  return;
}


/*==================================================================================================================================*/

int BPconditions(const int pntdim, double *y, int (*fnc)(double *, double *), int method, const int bpcurve, double *res, const double tol)

/*
   * BPconditions - Routine computes the additional conditions determining the location of
   *                a branching point, see the Matcont documentation (Branch point locator,
   *                page 36, eq. 41)
   *
   * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
   *                        equals 2 (for the bifurcation parameters) plus the dimension of
   *                        the vector of state variables, when bpcurve = 1, but equals
   *                        1 (for the single bifurcation parameter) plus the dimension of
   *                        the vector of state variables, when bpcurve = 0
   *              y       : Pointer to an array containing as first element the value
   *                        of the parameter p and as subsequent elements the values of
   *                        the state variables y. The last element is assumed to be the
   *                        second parameter in the BP continuation
   *              fnc     : Pointer to function specifying the system of
   *                        equations. The function must have a (double)
   *                        pointer as first argument, containing the point
   *                        in which to evaluate the system and a (double)
   *                        pointer as second argument, containing the
   *                        results after evaluation of the equations.
   *              method  : Method to use for differential computation: FORWARD or CENTRAL
   *              bpcurve : Routine called for detection of BP in EQ curve (0) or
   *                        during computation of BP curve (1)
   *              res     : Pointer to the result vector to adjust and/or compute. If NULL
   *                        the routine is used to initialize the additional variables for
   *                        the BP continuation
   */
{
  register int  i, j;
  const int     bppntdim = pntdim - bpcurve;
  int           rhsdim = bppntdim - 1, retval = SUCCES;
  double        *basemem, *rhs, *Jac;
  double        eval, *evec;
  double        *jc, *evalr, *evalc;

  // Prevent recurrence
  if (BPImmediateReturn) return retval;

  // Some of the following vectors are over-sized, but that's OK
  rhs = basemem = calloc(pntdim + pntdim*pntdim + 2*rhsdim + rhsdim*rhsdim, sizeof(double));
  if (!basemem) return ReportMemError("BPcondition");

  BPImmediateReturn = 1;

  Jac   = rhs + pntdim;
  evalr = Jac + pntdim*pntdim;
  evalc = evalr + rhsdim;
  jc    = evalc + rhsdim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  // Notice that when continuing a BP curve (bpcurve = 1) we have to call Jacobian() with the
  // full dimension (pntdim) of y to pass the entire argument vector to Equation().
  // As a result, the Jacobian matrix will have one row more than we really need,
  // representing the derivatives w.r.t. to the 2nd parameter of the BP continuation.
  // This additional row, with index pntdim-1 will be ignored.
  // Also notice that the last column will be unassigned, as the systems of equations returned
  // only has size rhsdim = pntdim-bpcurve-1. This last column is hence explicitly set to 0
  Jacobian(pntdim, y, bppntdim, Jac, fnc, method);
  for (j = 0; j < bppntdim; j++) Jac[j*bppntdim + rhsdim] = 0.0;

  /*
    * The resulting Jacobian equals the following matrix of partial derivatives:
    *
    *           |dF1/dp1 ... dFn/dp1  0|
    *           |dF1/dy1 ... dFn/dy1  0|
    *           |   .           .     0|
    *      Df = |   .           .     0|
    *           |   .           .     0|
    *           |dF1/dyn ... dFn/dyn  0|
    *           |dF1/dp2 ... dFn/dp2  0|
    *
    * when bpcurve = 1, and otherwise (bpcurve = 0)
    *
    *           |dF1/dp1 ... dFn/dp1  0|
    *           |dF1/dy1 ... dFn/dy1  0|
    *           |   .           .     0|
    *      Df = |   .           .     0|
    *           |   .           .     0|
    *           |dF1/dyn ... dFn/dyn  0|
    *
    * In which n = pntdim-bpcurve-1. Notice that all coefficients pertaining to yi are to be found
    * in ROW i (as opposed to column i). The Jacobian matrix is hence stored in column-wise
    * (fortran-style) form. However, the routine Eigenval() expects it to be in row-wise (C-style)
    * form. Hence, we need below the transpose of the Jacobian matrix computed above.
    */

  // Additional call to reset global variables
  (*fnc)(y, rhs);

  /*
   *  The extended system to solve for is:
   *
   *     F(x, p) + b*v   = 0
   *     (F_x(x, p))^T v = 0
   *     v^T F_p(x, p)   = 0
   *     v^T v - 1       = 0
   *
   *     with initial conditions b = 0 and v the eigenvector of the matrix (F_x(x, p))^T pertaining to
   *     the eigenvalue with the smallest norm. The unknowns are:
   *
   *     p:  the bifurcation parameter
   *     x:  the solution point
   *     b:  an additional value
   *     v:  the eigenvector (same dimension as x)
   *
   */
  if (res)
    {
      // Adjust the base equations
      eval = y[pntdim];
      evec = y + pntdim + 1;

      //  F(x, p) + b*v   = 0
      for (j = 0; j < rhsdim; j++) res[j] += eval*evec[j];

      // (F_x(x, p))^T v = 0
      for (i = 1; i < rhsdim + 1; i++, j++) res[j] = DOT(rhsdim, evec, 1, Jac + i*bppntdim, 1);

      // v^T F_p(x, p)   = 0
      res[j++] = DOT(rhsdim, evec, 1, Jac, 1);

      // v^T v - 1       = 0
      res[j++] = DOT(rhsdim, evec, 1, evec, 1) - 1;
    }
  else
    {
      // Extract the restricted Jacobian in transposed (row-wise, C-style) form
      for (i = 0; i < rhsdim; i++)
        for (j = 0; j < rhsdim; j++) jc[i*rhsdim + j] = Jac[(j + 1)*bppntdim + i];

      // Find the eigenvector pertaining to the eigenvalue with smallest absolute value
      retval = Eigenval(rhsdim, jc, 0, y + pntdim, MINIMUMNORM, y + pntdim + 1, NULL, tol);
      // Initialize the eigenvalue estimate to 0
      y[pntdim] = 0;
    }

  free(basemem);

  BPImmediateReturn = 0;

  return retval;
}


/*==================================================================================================================================*/

int CurveFuncDeriv(const int pntdim, double *pnt, double *dy, const int fncdim, double *dfuncdp, int (*fnc)(double *, double *), int method,
                   const double dytol)

/*
   * CurveFuncDeriv  -  Routine determines the derivative of the variables y w.r.t.
   *                    to the curve parameter p along the curve defined by the
   *                    system of equations
   *
   *                        F(y,p) = 0
   *
   *                    The vector y should have the same dimension as the number
   *                    of equations (i.e. the dimension of F(y,p)), which should
   *                    be equal to pntdim-1. The vector argument 'pnt' contains
   *                    the value of y and as last element the value of the parameter
   *                    p and hence has a dimension equal to 'pntdim'.
   *                    In addition, the routine determines the derivative w.r.t.
   *                    the curve parameter p of functions G(y, p) of the
   *                    variables y and the parameter p. The dimension of the
   *                    function vector G is determined by 'fncdim'
   *
   * Arguments -  pntdim  : The dimension of the argument vector 'pnt'.
   *              pnt     : Pointer to an array containing as first elements the
   *                        the state variables y and as the last element the
   *                        value of the parameter p. Together 'pnt' contains
   *                        the coordinates of a point on the curve determined
   *                        by F(y, p).
   *              dy      : Pointer to return vector with dy/dp values
   *              fncdim  : The dimension of the function vector G(y, p)
   *              dfuncdp : Point to return vector with dG/dp values
   *              fnc     : Pointer to function specifying the system of
   *                        equations. The function must have a (double)
   *                        pointer as first argument, containing the point
   *                        in which to evaluate the system and a (double)
   *                        pointer as second argument, containing the
   *                        results after evaluation of the equations.
   *                        This array of results should contain as first elements
   *                        the values of F(y,p), following by the values of G(y, p).
   *              method  : FORWARD, CENTRAL or RICHARDSON, determining the manner to
   *                        compute the derivative
   */

{
  register int  j, k;
  int           rhsdim = pntdim - 1, retcode;
  double        *basemem, *J;
  double        *Jac;
  double        *dfuncdy;

  // Prevent recurrence
  if (EXTImmediateReturn) return 0.0;

  J = basemem = calloc(pntdim*((pntdim - 1) + fncdim) + pntdim*(pntdim - 1) + pntdim*fncdim, sizeof(double));
  if (!basemem) return ReportMemError("CurveFuncDeriv");

  EXTImmediateReturn = 1;

  Jac     = J + pntdim*((pntdim - 1) + fncdim);
  dfuncdy = Jac + pntdim*(pntdim - 1);

  // Determine the Jacobian of the extended system (variable plus parameter
  // dependence). Use central differencing for the derivatives.
  // The derivatives w.r.t. to the parameter p are stored in the last row, because
  // the last vector element of the argument 'pnt' is the parameter.
  if (Jacobian(pntdim, pnt, rhsdim + fncdim, J, fnc, method) == FAILURE)
    {
      free(basemem);
      EXTImmediateReturn = 0;
      return FAILED_EVALUATION;
    }

#if (DEBUG == 1)
  int m, n;
  fprintf(stderr, "\nJacobian:\n=========\n");
  for (m = 0; m < pntdim; m++)
    {
      for (n = 0; n < (rhsdim + fncdim); n++) fprintf(stderr, "\t%14.8G", J[m*(rhsdim + fncdim) + n]);
      fprintf(stderr, "\n");
    }
#endif

  // The value dG(y(p), p)/dp is computed by:
  //
  //    dG(y(p), p)/dp = del(G(y(p), p))/del(y) dy/dp + del(G(y(p), p))/del(p)
  //
  // where dy/dp is solved from:
  //
  //   del(F(y(p),p))/del(y) dy/dp + del(F(y(p), p))/del(p) = 0
  //
  for (j = 0; j < pntdim; j++)
    {
      COPY(rhsdim, J + j*(rhsdim + fncdim), 1, Jac + j*rhsdim, 1);                  // Extract dF/dy & dF/dp
      COPY(fncdim, J + j*(rhsdim + fncdim) + rhsdim, 1, dfuncdy + j*fncdim, 1);     // Extract dG/dy & dG/dp
    }
  COPY(rhsdim, Jac + (pntdim - 1)*rhsdim, 1, dy, 1);                                // Store -dF/dp in dy
  SCAL(rhsdim, -1.0, dy, 1);

  // dy/dp is solved from:
  //
  //    del(F(y(p),p))/del(y) dy/dp + del(F(y(p), p))/del(p) = 0
  //
  // in code:
  //
  //    DF[1..n][1..n].dy = -DF[n+1][1..n] with dy = (dy1/dp ... dyn/dp)

  // Solve the linear system
  retcode = SolveLinearSystem(rhsdim, Jac, dy, dytol);
  if (retcode != SUCCES)
    {
      ErrorMsg(__FILE__, __LINE__, "Failed to solve linear system in TangentVec()");
      free(basemem);
      EXTImmediateReturn = 0;
      return retcode;
    }

  /*
   * Now dy contains the derivatives (dy1/dp, ..., dyn/dp) of the variables w.r.t. to the
   * curve parameter p. To compute dG/dp the matrix DG[1..n][1..m] is multiplied with the
   * vector dy and the last row of the matrix DG[n+1][1..m] is added to this result
   */
#if (DEBUG == 1)
  fprintf(stderr, "\ndY/dP:\n======\n");
  for (m = 0; m < rhsdim; m++) fprintf(stderr, "\t%14.8G", dy[m]);
  fprintf(stderr, "\n");
#endif

  // The value dG(y(p), p)/dp is computed by:
  //
  //    dG(y(p), p)/dp = del(G(y(p), p))/del(y) dy/dp + del(G(y(p), p))/del(p)
  //
  for (k = 0; k < fncdim; k++)
    {
      dfuncdp[k] = 0.0;
      for (j = 0; j < rhsdim; j++)
        {
          dfuncdp[k] += dfuncdy[j*fncdim + k]*dy[j];
        }
      dfuncdp[k] += dfuncdy[j*fncdim + k];                                          // j has value rhsdim, the index of dGi/dp
    }

#if (DEBUG == 1)
  fprintf(stderr, "\ndG/dy:\n======\n");
  for (m = 0; m < rhsdim; m++) fprintf(stderr, "\t%14.8G", dfuncdy[m]);
  fprintf(stderr, "\n");

  fprintf(stderr, "\ndG/dP:\n======\n");
  for (m = 0; m < fncdim; m++) fprintf(stderr, "\t%14.8G", dfuncdp[m]);
  fprintf(stderr, "\n\n");
#endif

  free(basemem);
  EXTImmediateReturn = 0;
  return SUCCES;
}


/*==================================================================================================================================*/
