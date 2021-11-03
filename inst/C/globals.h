/***
   NAME
     globals
   DESCRIPTION
     Header file with global definitions, such as error codes, return
     values and function prototypes.

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

    Last modification: AMdR - Aug 11, 2021
***/
#ifndef GLOBALS
#define GLOBALS

#if defined(__APPLE__) && !defined(__clang__)
// Work-arounds for using GCC with Accelerate framework on Mac OS Yosemite
// First address a bug in <os/base.h>. It should protect against __has_extension() being undefined. Provide a phony definition of __has_extension()
#ifndef __has_extension
#define __has_extension(x)        0
#endif
// Second, GCC doesn't support blocks. The following prevents using vImage features (simplified interoperability with Core Graphics and Core Video) that can not be used by GCC anyway
#define vImage_Utilities_h
#define vImage_CVUtilities_h
#endif

#include "ctype.h"
#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdarg.h"
#include "stdint.h"
#include "string.h"
#include <sys/stat.h>

#if defined(__linux__) || defined(__APPLE__)
#include "sys/param.h"                                                              // Defines MAXPATHLEN

#else
#if (!defined(MAXPATHLEN) && defined(PATH_MAX))
#define MAXPATHLEN PATH_MAX
#else
#define MAXPATHLEN 1024
#endif
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

#if defined(_MSC_VER)
#define INLINE __forceinline                                                        // use __forceinline (VC++ specific)
#else
#define INLINE                    static inline                                     // use standard inline
#endif

#if !defined(TRUE)
#if defined(true)
#define TRUE                      true
#else
#define TRUE                      1
#endif
#endif

/*
 *====================================================================================================================================
 *  Various other macro definitions
 *====================================================================================================================================
 */
#define MAX_STR_LEN               1024
#define MAX_STEPREDUCE            2048
#define MICRO                     1.0E-6

#define max(a, b)                 (((a) > (b)) ? (a) : (b))
#define min(a, b)                 (((a) < (b)) ? (a) : (b))
#define sign(a)                   (((a) < 0.0) ? (-1.0) : (1.0))

#if (defined(BIFTEST) || defined(CURVE) || defined(IO) || defined(ODE_DOPRI) || defined(NLEQ_H))
#undef EXTERN
#define EXTERN                    extern
#else
#undef EXTERN
#define EXTERN

#endif

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#define STDOUT                    mexPrintf
#elif defined(R_PACKAGE)
#define R_USE_C99_IN_CXX
#include "R.h"
#include <Rinternals.h>

#include "R_ext/Print.h"
#include <R_ext/Utils.h>

#define STDOUT                    Rprintf
#else
#define STDOUT                    printf
#endif

/*
 *====================================================================================================================================
 *  Macros used in curve.c
 *====================================================================================================================================
 */
#ifndef MAXITER
#define MAXITER                   40
#endif
#ifndef MAX_OUTPUTDIM
#define MAX_OUTPUTDIM             100
#endif

#define FAILURE                   0
#define SUCCES                    1

#define SUCCES_ODE_MAXTIME        10
#define SUCCES_ODE_EVENT          11

#define SINGULARITY               100
#define NORM_OVERFLOW             101
#define NO_CONVERGENCE            102
#define ILLEGAL_INPUT             103
#define FAILED_EVALUATION         104

#define FORWARD                   0                                                 // Computational methods for derivatives
#define CENTRAL                   1

#define UNDEFINED                 1000
#define BP                        1001
#define EQ                        1002
#define LP                        1003
#define EXT                       1004

#define DOMINANT                  10000
#define MINIMUMNORM               10001

EXTERN double                     Jacobian_Min_Step;
EXTERN int                        Jacobian_Updates;
EXTERN double                     Jacobian_Step;

EXTERN double                     Odesolve_Init_Step;                               // Dopri5 initial step size
EXTERN double                     Odesolve_Fixed_Step;                              // Dopri5 fixed step size
EXTERN double                     Odesolve_Min_Step;                                // Dopri5 minimum step size
EXTERN double                     Odesolve_Max_Step;                                // Dopri5 maximum step size
EXTERN double                     Odesolve_Abs_Err;                                 // Dopri5 absolute error in step
EXTERN double                     Odesolve_Rel_Err;                                 // Dopri5 relative error in size
EXTERN double                     Odesolve_Func_Tol;                                // Dopri5 function tolerance

EXTERN int                        Bifparone;
EXTERN int                        Bifpartwo;


/*
 *====================================================================================================================================
 *  Global variable definitions
 *====================================================================================================================================
 */

EXTERN double                     *odesolveMem;
EXTERN double                     *pnt_scale;
EXTERN FILE                       *biffile;
EXTERN FILE                       *errfile;
EXTERN FILE                       *outfile;
EXTERN double                     epsMach;
EXTERN double                     curvestep, Maxcurvestep;
EXTERN int                        Stepchange;
EXTERN int                        Stepreduce;
EXTERN int                        CurveType;
EXTERN int                        CurveEnd;
EXTERN int                        ConsoleSilent;
EXTERN int                        LocalizeType;
EXTERN int                        DoOutput;
EXTERN double                     Output[MAX_OUTPUTDIM];

EXTERN int                        EXTfun, EXTpar;


/*
 *====================================================================================================================================
 *  Function prototypeing: in main source file
 *====================================================================================================================================
 */

int                               DefineOutput(double *x, double *output);


/*
 *====================================================================================================================================
 *  Function prototypeing: biftest.c
 *====================================================================================================================================
 */

int                               LocateLP(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol);
//int                             LocateBP(int *dimpntr, double *y, int (*fnc)(double *, double *), double dytol, double rhstol, double *par2);
int                               LocateBP(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol);
int                               LocateEXT(const int pntdim, double *y, int (*fnc)(double *, double *), double dytol, double rhstol);


/*
 *====================================================================================================================================
 *  Function prototypeing: curve.c
 *====================================================================================================================================
 */

double                            anorm(int, int, double *);
int                               SetScales(double *point, int basedim);
int                               FindPoint(const int pntdim, double *guess, double *JacImport, double *tanvec, double ytol, double rhstol,
                                            const int max_iter, int (*fnc)(double *, double *));
int                               TangentVec(const int pntdim, double *sol, double *JacExport, double *tanvec, int (*fnc)(double *, double *),
                                             double *det, const double ytol);
int                               Jacobian(const int pntdim, double *pnt, const int fncdim, double *jac, int (*fnc)(double *, double *), int method);
int                               LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), const int method, const int lpcurve,
                                              double *curvedir, const double ytol);
int                               Determinant(const int N, double *M, double *det, double *cond);
int                               Eigenval(const int N, double *X, const int symmetric, double *eigval, const int eigenvaltype, double *righteigvec,
                                           double *lefteigvec, const double tol);
int                               SolveLinearSystem(const int N, double *A, double *B, double tol);
void                              ResetCurve(void);
int                               BPconditions(const int pntdim, double *y, int (*fnc)(double *, double *), int method, const int bpcurve, double *res,
                                               const double tol);
int                               CurveFuncDeriv(const int pntdim, double *pnt, double *dy, const int fncdim, double *dfuncdp,
                                                 int (*fnc)(double *, double *), int method, const double dytol);


/*
 *====================================================================================================================================
 *  Function prototypeing: dopri5.c
 *====================================================================================================================================
 */

int                               odesolve(double *, int, double *, double, void (*)(double, double *, double *), double (*)(double, double *));

/*
 *====================================================================================================================================
 *  Function prototypeing: io.c
 *====================================================================================================================================
 */

void                              ReportMsg(const char *, ...);
void                              ErrorMsg(const char *, const int, const char *, ...);
int                               ReportMemError(const char *name);
void                              NumProcError(const char *, const int, const int);
void                              PrettyPrintArray(FILE *, const int, double *);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE) || defined(R_PACKAGE)
// For Ctrl-C detection
int                               checkInterrupt(void);

EXTERN int                        CtrlCPressed;
#endif


/*
 *====================================================================================================================================
 *  Single/double precision and BLAS/LAPACK function mapping
 *====================================================================================================================================
 */

#if !defined(MATLAB_MEX_FILE) && !defined(R_PACKAGE)                                // Command-line usage on Mac OS & Linux or Octave on Mac OS

#ifdef OCTAVE_MEX_FILE
#include "mex.h"
#endif

#if defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#else
#include "cblas.h"
#endif

// BLAS Level 1 functions

#define ASUM                      cblas_dasum
#define AXPY                      cblas_daxpy
#define COPY                      cblas_dcopy
#define DOT                       cblas_ddot
#define NRM2                      cblas_dnrm2
#define SCAL                      cblas_dscal

// Lapack index type and function equivalences

#define LAPACK_SIZE_T             int

#else                                                                               // Matlab or R

#if defined(MATLAB_MEX_FILE)                                                        // Matlab on Windows or Mac OS
#include "mex.h"
#include "blas.h"
#include "lapack.h"

// Blas & Lapack index type

#define BLAS_SIZE_T               ptrdiff_t
#define LAPACK_SIZE_T             ptrdiff_t


#elif defined(R_PACKAGE)                                                            // R under Mac OS or under Windows using gcc from Rtools
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"

// Blas & Lapack index type

#define BLAS_SIZE_T               const int
#define LAPACK_SIZE_T             int

// Normalise BLAS Level 1 function names

#define dasum                     dasum_
#define daxpy                     daxpy_
#define dcopy                     dcopy_
#define ddot                      ddot_
#define dnrm2                     dnrm2_
#define dscal                     dscal_

#endif

// BLAS Level 1 functions

INLINE double ASUM(int a, double *b, int c)
{
  BLAS_SIZE_T aa = a, cc = c;
  return dasum(&aa, b, &cc);
}

INLINE void AXPY(int a, double b, double *c, int d, double *e, int f)
{
  BLAS_SIZE_T aa = a, dd = d, ff = f;
  double      bb = b;
  daxpy(&aa, &bb, c, &dd, e, &ff);
}

INLINE void COPY(int a, double *b, int c, double *d, int e)
{
  BLAS_SIZE_T aa = a, cc = c, ee = e;
  dcopy(&aa, b, &cc, d, &ee);
}

INLINE double DOT(int a, double *b, int c, double *d, int e)
{
  BLAS_SIZE_T aa = a, cc = c, ee = e;
  return ddot(&aa, b, &cc, d, &ee);
}

INLINE double NRM2(int a, double *b, int c)
{
  BLAS_SIZE_T aa = a, cc = c;
  return dnrm2(&aa, b, &cc);
}

INLINE void SCAL(int a, double b, double *c, int d)
{
  BLAS_SIZE_T aa = a, dd = d;
  double      bb = b;
  dscal(&aa, &bb, c, &dd);
}

#endif

/*
 * Routines and operations for which Lapack functions are used:
 *
 * Solving a system of Linear equations:
 *        FindPoint()   - Compute adjustment in Newton iteration
 *        TangentVec()  - Compute tangent vector of solution curve using the full Jacobian
 *        LPcondition() - Compute tangent vector of solution curve using the full Jacobian
 *
 * Determine the largest eigenvalue:
 *        Eigenval()    - To determine the largest eigenvalue and corresponding eigenvector of next generation matrix
 *                        To determine the largest eigenvalue of the Hermitian and Jacobian matrix of selection matrix
 *
 * Determine matrix determinant and/or condition:
 *        LPcondition() - Determine matrix condition to find most non-singular matrix to solve for tangent
 *        TangentVec()  - Determine matrix determinant for BP detection
 */

#if !defined(MATLAB_MEX_FILE)

#define dgetrf                    dgetrf_                                           // Used in: Determinant()
#define dgecon                    dgecon_                                           // Used in: Determinant()
#define dgeevx                    dgeevx_                                           // Used in: Eigenval()
#define dgesvx                    dgesvx_                                           // Used in: SolveLinearSystem()
#define dlamch                    dlamch_                                           // Used in: Eigenval()
#define dsyevr                    dsyevr_                                           // Used in: Eigenval()

#if defined(__linux__) && !defined(R_PACKAGE)
extern void                       dgetrf(LAPACK_SIZE_T *m, LAPACK_SIZE_T *n, double *a, LAPACK_SIZE_T *lda, LAPACK_SIZE_T *ipiv, LAPACK_SIZE_T *info);
extern void                       dgecon(char *norm, LAPACK_SIZE_T *n, double *a, LAPACK_SIZE_T *lda, double *anorm, double *rcond, double *work,
                                         LAPACK_SIZE_T *iwork, LAPACK_SIZE_T *info);
extern void                       dgeevx(char *balanc, char *jobvl, char *jobvr, char *sense, LAPACK_SIZE_T *n, double *a, LAPACK_SIZE_T *lda,
                                         double *wr, double *wi, double *vl, LAPACK_SIZE_T *ldvl, double *vr, LAPACK_SIZE_T *ldvr, LAPACK_SIZE_T *ilo,
                                          LAPACK_SIZE_T *ihi, double *scale, double *abnrm, double *rconde, double *rcondv, double *work,
                                          LAPACK_SIZE_T *lwork, LAPACK_SIZE_T *iwork, LAPACK_SIZE_T *info);
extern void                       dgesvx(char *fact, char *trans, LAPACK_SIZE_T *n, LAPACK_SIZE_T *nrhs, double *a, LAPACK_SIZE_T *lda, double *af,
                                         LAPACK_SIZE_T *ldaf, LAPACK_SIZE_T *ipiv, char *equed, double *r, double *c, double *b, LAPACK_SIZE_T *ldb,
                                         double *x, LAPACK_SIZE_T *ldx, double *rcond, double *ferr, double *berr, double *work, LAPACK_SIZE_T *iwork,
                                         LAPACK_SIZE_T *info);
extern double                     dlamch(char *cmach);
extern void                       dsyevr(char *jobz, char *range, char *uplo, LAPACK_SIZE_T *n, double *a, LAPACK_SIZE_T *lda, double *vl, double *vu,
                                         LAPACK_SIZE_T *il, LAPACK_SIZE_T *iu, double *abstol, LAPACK_SIZE_T *m, double *w, double *z, LAPACK_SIZE_T *ldz,
                                         LAPACK_SIZE_T *isuppz, double *work, LAPACK_SIZE_T *lwork, LAPACK_SIZE_T *iwork, LAPACK_SIZE_T *liwork,
                                         LAPACK_SIZE_T *info);
#endif

#endif


/*==================================================================================================================================*/

#endif /* GLOBALS */
