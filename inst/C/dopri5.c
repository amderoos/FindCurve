/*
  NAME
    dopri5
  DESCRIPTION
    DOPRI5 solver of an ODE system

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

  Last modification: AMdR - May 06, 2020
*/
/*
 *====================================================================================================================================
 * Macro definitions that can be (p)re-defined by the user
 *====================================================================================================================================
 */

#ifndef ODE_DOPRI
#define ODE_DOPRI
#include "globals.h"


/*
 *====================================================================================================================================
 * Macro definitions needed in integration routine (system independent)
 *====================================================================================================================================
 */

#ifndef BETAFAC
#define BETAFAC                   0.04
#endif
#ifndef FAC1
#define FAC1                      0.2
#endif
#ifndef FAC2
#define FAC2                      10.0
#endif
#ifndef SAFETY
#define SAFETY                    0.9
#endif
#ifndef FACOLD
#define FACOLD                    1.0E-4
#endif


/*
 *====================================================================================================================================
 * Global variables needed in integration routine (system independent)
 *====================================================================================================================================
 */

#define VECTOR_COPIES             15

static double                     old_x, xval, del_h;
static double                     *y, *yy1, *ytmp, *old_y;
static double                     *k1, *k2, *k3, *k4, *k5, *k6;
static double                     *rcont1, *rcont2, *rcont3, *rcont4, *rcont5;

static int                        Ode_Dim;


/*==================================================================================================================================*/

static void dopri5(void (*rhs)(double, double *, double *), double dt)

/*
 * dopri5 - Routine performs an integration step of the system specified in the
 *          function (*rhs)()using the DOPRI5 integration method. The code is
 *          adapted from the original C source code by E. Hairer & G. Wanner.
 */

{
  register int        i;
  static const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0, a21 = 0.2, a31 = 3.0/40.0, a32 = 9.0/40.0, a41 = 44.0/45.0, a42 = -56.0/15.0,
                      a43 = 32.0/9.0, a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0, a61 = 9017.0/3168.0,
                      a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0, a71 = 35.0/384.0, a73 = 500.0/1113.0,
                      a74 = 125.0/192.0, a75 = -2187.0/6784.0, a76 = 11.0/84.0, d1 = -12715105075.0/11282082432.0, d3 = 87487479700.0/32700410799.0,
                      d4 = -10690763975.0/1880347072.0, d5 = 701980252875.0/199316789632.0, d6 = -1453857185.0/822651844.0,
                      d7 = 69997945.0/29380423.0, e1 = 71.0/57600.0, e3 = -71.0/16695.0, e4 = 71.0/1920.0, e5 = -17253.0/339200.0, e6 = 22.0/525.0,
                      e7 = -1.0/40.0;

  (*rhs)(xval, y, k1);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*a21*k1[i];

  (*rhs)(xval + c2*dt, yy1, k2);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*(a31*k1[i] + a32*k2[i]);

  (*rhs)(xval + c3*dt, yy1, k3);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*(a41*k1[i] + a42*k2[i] + a43*k3[i]);

  (*rhs)(xval + c4*dt, yy1, k4);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);

  (*rhs)(xval + c5*dt, yy1, k5);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);

  (*rhs)(xval + dt, yy1, k6);

  for (i = 0; i < Ode_Dim; i++) yy1[i] = y[i] + dt*(a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);

  (*rhs)(xval + dt, yy1, k2);

  for (i = 0; i < Ode_Dim; i++) rcont5[i] = dt*(d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);

  for (i = 0; i < Ode_Dim; i++) k4[i] = dt*(e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k2[i]);

  return;
} /* dopri5 */


/*==================================================================================================================================*/
#define ITMAX                     500
#define EPS                       1.0e-16

static double zbrent(double *left, double *right, double (*stop)(double, double *))

{
  register int  i;
  int           iter;
  double        a, b, c = 0.0, d = 0.0, e = 0.0, min1, min2;
  double        fa, fb, fc, p, q, r, s, tol1, xm;

  a  = 0.0;
  fa = (*stop)(old_x, left);
  b  = 1.0;
  fb = (*stop)(old_x + del_h, right);

  if (fb*fa > 0.0) return -1.0;

  fc = fb;
  for (iter = 0; iter < ITMAX; iter++)
    {
      if (fb*fc > 0.0)
        {
          c  = a;
          fc = fa;
          e = d = b - a;
        }
      if (fabs(fc) < fabs(fb))
        {
          a  = b;
          fa = fb;
          b  = c;
          fb = fc;
          c  = a;
          fc = fa;
        }
      tol1 = 2.0*EPS*fabs(b) + 0.5*Odesolve_Func_Tol;
      xm   = 0.5*(c - b);

      if ((fabs(xm) <= tol1 && fb*fc <= 0.0) || fb == 0.0) return b;

      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
          s = fb/fa;
          if (a == c)
            {
              p = 2.0*xm*s;
              q = 1.0 - s;
            }
          else
            {
              q = fa/fc;
              r = fb/fc;
              p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
              q = (q - 1.0)*(r - 1.0)*(s - 1.0);
            }
          if (p > 0.0) q = -q;
          p              = fabs(p);
          min1           = 3.0*xm*q - fabs(tol1*q);
          min2           = fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2))
            {
              e = d;
              d = p/q;
            }
          else
            {
              d = xm;
              e = d;
            }
        }
      else
        {
          d = xm;
          e = d;
        }
      a  = b;
      fa = fb;
      if (fabs(d) > tol1)
        b += d;
      else
        b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
      for (i = 0; i < Ode_Dim; i++) ytmp[i] = rcont1[i] + b*(rcont2[i] + (1.0 - b)*(rcont3[i] + b*(rcont4[i] + (1.0 - b)*rcont5[i])));
      fb                                    = (*stop)(old_x + b*del_h, ytmp);
    }

  return -1.0;
}


#undef ITMAX
#undef EPS
/*==================================================================================================================================*/

int odesolve(double *yvec, int vecdim, double *xinit, double xmax, void (*rhs)(double, double *, double *), double (*stop)(double, double *))

/*
 * odesolve - Solve the system of ODEs specified in the function (*rhs)() up to
 *            a maximum time xmax or until the function (*stop)() returns a 0 result
 */

{
  register int  i;
  double        oldstop = -1.0, newstop = -1.0;
  double        err, sk, sqr, errmax;
  double        fac, fac11, hnew;
  double        yd0, ydiff, bspl, facold = FACOLD;
  double        htmp, xnext;
  int           adjust_step = 1, stopped = 0, ierr;
  static long   odesolveMemAllocated = 0L;

  Ode_Dim = vecdim;
  if ((VECTOR_COPIES*Ode_Dim) > odesolveMemAllocated)
    {
      odesolveMem = (double *)realloc(odesolveMem, VECTOR_COPIES*Ode_Dim*sizeof(double));
      if (!odesolveMem) return ReportMemError("odesolve");
      odesolveMemAllocated = VECTOR_COPIES*Ode_Dim;
    }

  y      = odesolveMem;
  yy1    = y + Ode_Dim;
  ytmp   = yy1 + Ode_Dim;
  old_y  = ytmp + Ode_Dim;
  k1     = old_y + Ode_Dim;
  k2     = k1 + Ode_Dim;
  k3     = k2 + Ode_Dim;
  k4     = k3 + Ode_Dim;
  k5     = k4 + Ode_Dim;
  k6     = k5 + Ode_Dim;
  rcont1 = k6 + Ode_Dim;
  rcont2 = rcont1 + Ode_Dim;
  rcont3 = rcont2 + Ode_Dim;
  rcont4 = rcont3 + Ode_Dim;
  rcont5 = rcont4 + Ode_Dim;

  xval =*xinit;
  (void)memcpy((void *)y, (void *)yvec, Ode_Dim*sizeof(double));
  if (stop) newstop = oldstop = (*stop)(xval, yvec);
  hnew                        = Odesolve_Init_Step;

  while (xval < xmax - Odesolve_Abs_Err)
    {
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
      if (checkInterrupt()) return FAILURE;
#endif
      del_h = min(hnew, xmax - xval);
      xnext = Odesolve_Fixed_Step*(1 + floor((1 + Odesolve_Rel_Err)*xval/Odesolve_Fixed_Step));
      if (xnext < xval + del_h)
        {
          del_h       = xnext - xval;
          adjust_step = 0;
        }
      dopri5(rhs, del_h);                                                           // Do an integration step

      errmax = err = 0.0;                                                           // error estimation
      for (i = 0, ierr = 0; i < Ode_Dim; i++)
        {
          sk  = Odesolve_Abs_Err + Odesolve_Rel_Err*max(fabs(y[i]), fabs(yy1[i]));
          sqr = k4[i]/sk;
          if ((sqr * sqr) > errmax)
            {
              errmax = sqr * sqr;
              ierr = i;
            }
          err += sqr*sqr;
        }
      errmax /= err;
      errmax *= 100.0;
      err = sqrt(err/(double)Ode_Dim);

      if (err > 1.0)                                                                // Step rejected
        {
          if (hnew <= Odesolve_Min_Step)                                            // Minimum step has failed
            {
            ErrorMsg(__FILE__, __LINE__, "At t = %f %s\n%-23s%s (%.1f%%) %s %d.",
                     xval, "the error test failed repeatedly with the minimum step size.",
                     " ", "Maximum contribution", errmax, "by component", ierr);
            return FAILURE;
            }
          // If bigger than accuracy take smaller step and restart
          fac11       = pow(err, 0.2 - BETAFAC*0.75);
          hnew        = del_h/min(1.0/FAC1, fac11/SAFETY);
          hnew        = max(hnew, Odesolve_Min_Step);
          adjust_step = 0;
        }
      else                                                                          // Step accepted
        {
          // Step size adjustment using Lund-stabilization: we require
          //      fac1 <=  hnew/h <= fac2
          // No increase if failed just before
          if (adjust_step)
            {
              fac11  = pow(err, 0.2 - BETAFAC*0.75);
              fac    = fac11/pow(facold, BETAFAC);
              fac    = max(1.0/FAC2, min(1.0/FAC1, fac/SAFETY));
              hnew   = del_h/fac;
              facold = max(err, FACOLD);
              hnew   = min(hnew, Odesolve_Max_Step);
              hnew   = max(hnew, Odesolve_Min_Step);
            }
          adjust_step = 1;

          // Update variables for event location and continuous output
          if (stop)
            {
              for (i = 0; i < Ode_Dim; i++)
                {
                  yd0       = y[i];
                  ydiff     = yy1[i] - yd0;
                  bspl      = del_h*k1[i] - ydiff;
                  rcont1[i] = y[i];
                  rcont2[i] = ydiff;
                  rcont3[i] = bspl;
                  rcont4[i] = -del_h*k2[i] + ydiff - bspl;
                }
            }
          // Store previous solution point and update the basic data copy
          old_x = xval;
          xval += del_h; /* Update time value        */
          (void)memcpy((void *)old_y, (void *)y, Ode_Dim*sizeof(double));
          (void)memcpy((void *)y, (void *)yy1, Ode_Dim*sizeof(double));

          if (stop)
            {
              newstop = (*stop)(xval, y);
              stopped = ((oldstop*newstop) <= 0.0);
              if (stopped) break;
              oldstop = newstop;
            }
        }
    }

  if (stop && stopped)                                                              // Stop through treshold
    {
      htmp = zbrent(old_y, y, stop);
      if (htmp < 0.0)
        {
          ErrorMsg(__FILE__, __LINE__, "Root location in zbrent() failed. Last solution point kept.");
          return FAILURE;
        }

      *xinit = old_x + htmp*del_h;
      for (i    = 0; i < Ode_Dim; i++)
        yvec[i] = rcont1[i] + htmp*(rcont2[i] + (1.0 - htmp)*(rcont3[i] + htmp*(rcont4[i] + (1.0 - htmp)*rcont5[i])));
      return SUCCES_ODE_EVENT;
    }

  // Stop through maximum time
  *xinit = xval;
  (void)memcpy((void *)yvec, (void *)y, Ode_Dim*sizeof(double));

  return SUCCES_ODE_MAXTIME;
}


/*==================================================================================================================================*/
#undef BETAFAC
#undef FAC1
#undef FAC2
#undef SAFETY
#undef FACOLD
#undef VECTOR_COPIES

#endif
