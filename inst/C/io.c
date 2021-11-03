/***
   NAME
     io
   DESCRIPTION
     Implements all low-level I/O routines

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

    Last modification: AMdR - Mar 19, 2020
***/
#ifndef IO
#define IO
#endif
#include "globals.h"

/*==================================================================================================================================*/
/*
 * Start of function implementations.
 */

void ReportMsg(const char *fmt, ...)

{
  va_list argpnt;

  va_start(argpnt, fmt);
  if (errfile)
    {
      vfprintf(errfile, fmt, argpnt);
      (void)fflush(errfile);
    }
  va_end(argpnt);

  return;
}


/*==================================================================================================================================*/

void ErrorMsg(const char *name, const int line, const char *fmt, ...)

{
  va_list argpnt;

#if (defined(R_PACKAGE))
  if (!ConsoleSilent) REprintf("\n** %-12s (%3d): ", name, line);
#elif defined(MATLAB_MEX_FILE)
  if (!ConsoleSilent) mexPrintf("\n** %-12s (%3d): ", name, line);
#else
  if (!ConsoleSilent) (void)fprintf(stderr, "\n** %-12s (%3d): ", name, line);
#endif
  va_start(argpnt, fmt);
#if (defined(R_PACKAGE))
  if (!ConsoleSilent) REvprintf(fmt, argpnt);
#elif defined(MATLAB_MEX_FILE)
  if (!ConsoleSilent) mexPrintf(fmt, argpnt);
#else
  if (!ConsoleSilent) vfprintf(stderr, fmt, argpnt);
#endif
  va_end(argpnt);

#if (defined(R_PACKAGE))
  if (!ConsoleSilent) REprintf("\n");
#elif defined(MATLAB_MEX_FILE)
  if (!ConsoleSilent) mexPrintf("\n");
#else
  if (!ConsoleSilent) (void)fprintf(stderr, "\n");
#endif

  if (errfile)
    {
      (void)fprintf(errfile, "\n** %-12s (%3d): ", name, line);
      va_start(argpnt, fmt);
      vfprintf(errfile, fmt, argpnt);
      va_end(argpnt);
      (void)fprintf(errfile, "\n\n");
    }
  (void)fflush(NULL);
#ifdef MATLAB_MEX_FILE
  mexEvalString("pause(0.0001);");
#elif (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  return;
}

/*==================================================================================================================================*/

int ReportMemError(const char *name)

{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexWarnMsgIdAndTxt("MATLAB:MemoryError", "\nMemory allocation error in %s(). Exiting MEX module.....\n\n", name);
#else
  warning("\nMemory allocation error in %s(). Exiting module.....\n\n", name);
  R_FlushConsole();
  R_ProcessEvents();
#endif
  CtrlCPressed = TRUE;
#else
  ErrorMsg(__FILE__, __LINE__, "Memory allocation error in %s()!", name);
  exit(1);                                                                          // Only executed when in command-line mode
#endif

  return FAILURE;
}


/*==================================================================================================================================*/

void NumProcError(const char *name, const int line, const int NumProcErrorCode)

{
  switch (NumProcErrorCode)
    {
      case SINGULARITY:
        ErrorMsg(name, line, "%-45s", "Singular matrix encountered");
        break;
      case NORM_OVERFLOW:
        ErrorMsg(name, line, "%-45s", "Norm overflow in Newton iteration");
        break;
      case NO_CONVERGENCE:
        ErrorMsg(name, line, "%-45s", "No convergence in Newton iteration");
        break;
      case ILLEGAL_INPUT:
      case FAILED_EVALUATION:
        break;
      default:
        ErrorMsg(name, line, "%-45s%d", "Unknown numerical error code: ", NumProcErrorCode);
        break;
    }

  return;
}


/*==================================================================================================================================*/

void PrettyPrintArray(FILE *fp, const int dim, double *vec)

{
  register int  i;
  double        tmp;

  for (i = 0; i < dim; i++)
    {
      tmp = vec[i];
      if (((fabs(tmp) <= 1.0E4) && (fabs(tmp) >= 1.0E-3)) || (tmp == 0.))
        {
          if (fp)
            (void)fprintf(fp, "%16.8f", tmp);
          else
            (void)STDOUT("%16.8f", tmp);
        }
      else
        {
          if (fp)
            (void)fprintf(fp, "%16.8E", tmp);
          else
            (void)STDOUT("%16.8E", tmp);
        }
    }
  if (fp)
    (void)fprintf(fp, "\n");
  else
    (void)STDOUT("\n");
  (void)fflush(fp);
#if (defined(R_PACKAGE))
  R_FlushConsole();
  R_ProcessEvents();
#endif

  return;
}


/*==================================================================================================================================*/
// Ctrl-C detection

#if defined(MATLAB_MEX_FILE)
int checkInterrupt()
{
  int         pressed;
  extern bool utIsInterruptPending();
  extern void utSetInterruptPending(bool);

  // check for a Ctrl-C event
  pressed = utIsInterruptPending();
  if (pressed)
    {
      utSetInterruptPending(false);
      mexPrintf("\n\nCtrl-C detected. Stopping computation\n\n");
      CtrlCPressed = TRUE;
    }

  return (pressed || CtrlCPressed);
}

#elif defined(OCTAVE_MEX_FILE)

int checkInterrupt()
{
  int pressed = 0;

  // checking of Ctrl-C event not implemented

  return pressed;
}

#elif defined(R_PACKAGE)

static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }

// this will call the above in a top-level context so it won't longjmp-out of your context
int checkInterrupt()
{
  Rboolean R_ToplevelExec(void (*fun)(void *), void *data);

  if (CtrlCPressed) return CtrlCPressed;

#if (defined(OPENMP))
  if (omp_get_thread_num() == 0)
#endif
    {
      CtrlCPressed = (R_ToplevelExec(chkIntFn, NULL) == FALSE);
      if (CtrlCPressed) REprintf("\n\nUser interrupt detected. Stopping computation\n\n");
    }

  return (CtrlCPressed);
}

#endif

/*==================================================================================================================================*/
