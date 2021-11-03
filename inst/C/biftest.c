/***
   NAME
     biftest
   DESCRIPTION
     This module implements routines that locate bifurcation points along
     the curve that is being continued

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

    Last modification: AMdR - Jun 01, 2020
***/
#ifndef BIFTEST
#define BIFTEST
#endif
#include "globals.h"


/*
 *====================================================================================================================================
 *  Routines to locate bifurcation points
 *====================================================================================================================================
 */

int LocateLP(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol)

/*
 * LocateLP - Routine locates a limit point in an equilibrium bifurcation curve.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y = (p,x)'.
 *              y       : Pointer to an array containing as first element the value
 *                        of the parameter p and as subsequent elements the values of
 *                        the state variables x.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              ytol    : Tolerance determining when change in y equals zero.
 *              rhstol  : Tolerance determining when RHS equals zero.
 */
{
  int     i, retval, outmax;
  double  *basemem, *lppoint;

  lppoint = basemem = calloc(pntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateLP");

  COPY(pntdim, y, 1, lppoint, 1);

  ReportMsg("\nStarting LP location from :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", lppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  LocalizeType = LP;
  retval       = FindPoint(pntdim, lppoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  LocalizeType = UNDEFINED;

  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New LP point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", lppoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      STDOUT("%16.8E", lppoint[0]*pnt_scale[0]);
      for (i = 1; i < pntdim; i++) STDOUT(", %16.8E", lppoint[i]*pnt_scale[i]);
      STDOUT("  ****  LP   ****");
      STDOUT("\n");
      DoOutput = 1;
      outmax = DefineOutput(lppoint, Output);
      DoOutput = 0;
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  LP   ****");
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
    }
  else
    {
      ReportMsg("\nFailed to locate LP bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate LP bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate LP bifurcation point ****\n");
#endif
    }
  fflush(NULL);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/
#if (1)
int LocateBP(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol)

/*
   * LocateBP - Routine locates a transcritical bifurcation in the equation system. See the
   *            Matcont documentation (Branch point locator, page 36, eq. 41)
   *
   *   The extended system to solve for is:
   *
   *        F(x, p) + b*v   = 0
   *        (F_x(x, p))^T v = 0
   *        v^T F_p(x, p)   = 0
   *        v^T v - 1       = 0
   *
   *   with initial conditions b = 0 and v the eigenvector of the matrix (F_x(x, p))^T pertaining to
   *   the eigenvalue with the smallest norm. The unknowns are:
   *
   *   p: the bifurcation parameter
   *   x: the solution point
   *   b: an additional value
   *   v: the eigenvector (same dimension as x)
   *
   * Arguments -  pntdim  : The dimension of the argument vector 'y = (p,x)'.
   *              y       : Pointer to an array containing as first element the value
   *                        of the parameter p and as subsequent elements the values of
   *                        the state variables x.
   *              fnc     : Pointer to function specifying the system of
   *                        equations. The function must have a (double)
   *                        pointer as first argument, containing the point
   *                        in which to evaluate the system and a (double)
   *                        pointer as second argument, containing the
   *                        results after evaluation of the equations.
   *              ytol    : Tolerance determining when change in y equals zero.
   *              rhstol  : Tolerance determining when RHS equals zero.
   */
{
  int     i, retval, outmax;
  int     unknowns;
  double  *basemem, *bppoint;

  unknowns = 2*pntdim;                                                              // p, x, beta, v
  for (i = pntdim; i < unknowns; i++) pnt_scale[i] = 1.0;                           // Set the scale for the additional variables beta and v

  bppoint = basemem = calloc(unknowns, sizeof(double));
  if (!basemem) return ReportMemError("LocateBP");

  COPY(pntdim, y, 1, bppoint, 1);

  // Initialize additional variables in BP continuation
  BPconditions(pntdim, bppoint, fnc, CENTRAL, 0, NULL, rhstol);

  ReportMsg("\n\nStarting BP location from :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", bppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  LocalizeType = BP;
  retval       = FindPoint(unknowns, bppoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  LocalizeType = UNDEFINED;
  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New BP point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", bppoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      STDOUT("%16.8E", bppoint[0]*pnt_scale[0]);
      for (i = 1; i < pntdim; i++) STDOUT(", %16.8E", bppoint[i]*pnt_scale[i]);
      STDOUT("  ****  BP   ****");
      STDOUT("\n");
      DoOutput = 1;
      outmax = DefineOutput(bppoint, Output);
      DoOutput = 0;
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  BP   ****");
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
    }
  else
    {
      ReportMsg("\nFailed to locate BP bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate BP bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate BP bifurcation point ****\n");
#endif
    }

  fflush(NULL);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  free(basemem);

  return retval;
}
#else
int LocateBP(int *dimpntr, double *y, int (*fnc)(double *, double *), double dytol, double rhstol, double *par2)

/*
   * LocateBP - Routine locates a transcritical bifurcation in the equation system
   *
   * Arguments -  pntdim  : The dimension of the argument vector 'y = (p,x)'.
   *              y       : Pointer to an array containing as first element the value
   *                        of the parameter p and as subsequent elements the values of
   *                        the state variables x.
   *              fnc     : Pointer to function specifying the system of
   *                        equations. The function must have a (double)
   *                        pointer as first argument, containing the point
   *                        in which to evaluate the system and a (double)
   *                        pointer as second argument, containing the
   *                        results after evaluation of the equations.
   *              ytol    : Tolerance determining when change in y equals zero.
   *              rhstol  : Tolerance determining when RHS equals zero.
   *              par2    : Pointer to the value of the 2nd bifurcation parameter
   */
{
  int     i, retval, outmax;
  int     bppntdim, unknowns;
  int     OldCurveType, oldpntdim;
  double  *basemem, *bppoint, oldpar;
  double  *tv;
  double  *tmpscales, *savedpntr;

  bppntdim =*dimpntr + 1;
  unknowns = 2*(*dimpntr) + 1;

  bppoint = basemem = calloc(3*unknowns, sizeof(double));
  if (!basemem) return ReportMemError("LocateBP");

  tmpscales = bppoint + unknowns;
  tv        = tmpscales + unknowns;

  COPY(*dimpntr, pnt_scale, 1, tmpscales, 1);
  tmpscales[bppntdim - 1] = 1.0;
  tv[bppntdim - 1]        = 1.0;

  COPY(*dimpntr, y, 1, bppoint, 1);
  oldpar = bppoint[bppntdim - 1] =*par2;

  OldCurveType = CurveType;
  CurveType    = BP;

  savedpntr = pnt_scale;
  pnt_scale = tmpscales;

  oldpntdim =*dimpntr;
  *dimpntr  = bppntdim;

  // Initialize additional variables in BP continuation
  BPconditions(bppntdim, bppoint, fnc, CENTRAL, 1, NULL, dytol);

  ReportMsg("\n\nStarting BP location from :\t");
  for (i = 0; i < bppntdim; i++) ReportMsg("%10.5E\t", bppoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  retval = FindPoint(unknowns, bppoint, NULL, tv, dytol, rhstol, MAXITER, fnc);
  if (retval == SUCCES)
    {
      STDOUT("%16.8E", bppoint[0]*pnt_scale[0]);
      for (i = 1; i < oldpntdim; i++) STDOUT(", %16.8E", bppoint[i]*pnt_scale[i]);
      STDOUT("  ****  BP   ****");
      STDOUT("\n");
      DoOutput = 1;
      outmax = DefineOutput(bppoint, Output);
      DoOutput = 0;
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  BP   ****");
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }
  else
    {
      ReportMsg("\nFailed to locate BP bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate BP bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate BP bifurcation point ****\n");
#endif
      fflush(NULL);
#if (defined(R_PACKAGE))
      R_FlushConsole();
      R_ProcessEvents();
#endif
    }

  *par2     = oldpar;
  CurveType = OldCurveType;
  pnt_scale = savedpntr;
  *dimpntr  = oldpntdim;

  free(basemem);

  return retval;
}

#endif

/*==================================================================================================================================*/

int LocateEXT(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol)

/*
   * LocateEXT -  Routine locates a maximum or minimum value of a function of the
   *              solution variables along a bifurcation curve.
   *
   * Arguments -  pntdim  : The dimension of the argument vector 'y = (p,x)'.
   *              y       : Pointer to an array containing as first element the value
   *                        of the parameter p and as subsequent elements the values of
   *                        the state variables x.
   *              fnc     : Pointer to function specifying the system of
   *                        equations. The function must have a (double)
   *                        pointer as first argument, containing the point
   *                        in which to evaluate the system and a (double)
   *                        pointer as second argument, containing the
   *                        results after evaluation of the equations.
   *              ytol    : Tolerance determining when change in y equals zero.
   *              rhstol  : Tolerance determining when RHS equals zero.
   */
{
  int     i, retval, outmax;
  double  *basemem, *extpoint;

  extpoint = basemem = calloc(pntdim, sizeof(double));
  if (!basemem) return ReportMemError("LocateEXT");

  COPY(pntdim, y, 1, extpoint, 1);

  ReportMsg("\n\nStarting EXT location from :\t");
  for (i = 0; i < pntdim; i++) ReportMsg("%10.5E\t", extpoint[i]*pnt_scale[i]);
  ReportMsg("\n");

  LocalizeType = EXT;
  retval       = FindPoint(pntdim, extpoint, NULL, NULL, dytol, rhstol, MAXITER, fnc);
  LocalizeType = UNDEFINED;

  if (retval == SUCCES)
    {
      // Report on located point to stderr file
      ReportMsg("New EXT point :\t");
      for (i = 0; i < pntdim; i++) ReportMsg("%16.8E  ", extpoint[i]*pnt_scale[i]);
      ReportMsg("\n");

      STDOUT("%16.8E", extpoint[0]*pnt_scale[0]);
      for (i = 1; i < pntdim; i++) STDOUT(", %16.8E", extpoint[i]*pnt_scale[i]);
      STDOUT("  ****  EXT  ****");
      STDOUT("\n");
      DoOutput = 1;
      outmax = DefineOutput(extpoint, Output);
      DoOutput = 0;
      if (biffile && outmax)
        {
          for (i = 0; i < outmax; i++) fprintf(biffile, "%16.8E", Output[i]);
          fprintf(biffile, "  ****  EXT  ****");
          fprintf(biffile, "\n");
        }
      if (outfile && outmax) PrettyPrintArray(outfile, outmax, Output);
    }
  else
    {
      ReportMsg("\nFailed to locate EXT bifurcation point\n\n");
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      STDOUT("**** Failed to locate EXT bifurcation point ****\n");
#else
      fprintf(stderr, "**** Failed to locate EXT bifurcation point ****\n");
#endif
    }

  fflush(NULL);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  free(basemem);

  return retval;
}


/*==================================================================================================================================*/
