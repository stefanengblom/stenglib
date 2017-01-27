/* powerseries.c */

/* S. Engblom 2007-08-17 (Minor revision) */
/* S. Engblom 2007-06-13 (Minor revision) */
/* S. Engblom 2007-01-22 */

#include <math.h>

#include "mex.h"
#include "matrix.h"

/* forward declaration */
void powerseries(double *yr,double *yi,
		 const double *cr,const double *ci,int Nc,
		 const double *xr,const double *xi,int Nx,
		 double tol);

#define ISDOUBLEMATRIX(A) (mxIsDouble(A) && !mxIsSparse(A) && \
			   mxGetNumberOfDimensions(A) == 2)
#define ISREALSCALAR(A) (!mxIsCell(A) && !mxIsStruct(A) && \
			 !mxIsComplex(A) && mxGetNumberOfElements(A) == 1)

/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  if (1 < nlhs || nrhs != 3)
    mexErrMsgTxt("Expecting one output and 3 inputs.");

  /* inputs C and X */
  if (!ISDOUBLEMATRIX(prhs[0]) || !ISDOUBLEMATRIX(prhs[1]))
    mexErrMsgTxt("First two arguments must be double matrices.");

  /* input tol */
  if (!ISREALSCALAR(prhs[2]))
    mexErrMsgTxt("Expecting a real scalar tolerance.");

  /* create output */
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]),
				 mxGetN(prhs[1]),
				 mxIsComplex(prhs[0]) || 
				 mxIsComplex(prhs[1]) ? mxCOMPLEX : mxREAL);

  /* evaluate the result */
  powerseries(mxGetPr(plhs[0]),mxGetPi(plhs[0]),
	      mxGetPr(prhs[0]),mxGetPi(prhs[0]),
	      mxGetNumberOfElements(prhs[0]),
	      mxGetPr(prhs[1]),mxGetPi(prhs[1]),
	      mxGetNumberOfElements(prhs[1]),
	      *mxGetPr(prhs[2]));
}
/*------------------------------------------------------------------------*/
void powerseries(double *yr,double *yi,
		 const double *cr,const double *ci,int Nc,
		 const double *xr,const double *xi,int Nx,
		 double tol)
/* Summation of power series. Computes y = sum_(k >= 0) c[k]*x^k where
   y = (yr,yi), c = (cr,ci) and x = (xr,xi) and where any of the
   imaginary parts may be NULL. The coefficient vector has length Nc
   and the coordinate vector length Nx. It is required that y be
   allocated and cleared before the call.*/
{
  bool warn = false;

  if (yi == NULL)
    for (int i = 0; i < Nx; i++) {
      int j;
      double powr = 1.0;
      for (j = 0; j < Nc; j++) {
	double termr = cr[j]*powr;
	yr[i] += termr;
	if (fabs(termr) < fabs(yr[i])*tol) break;
	powr *= xr[i];
      }
      warn = warn || j == Nc;
    }
  else if (ci == NULL)
    for (int i = 0; i < Nx; i++) {
      int j;
      double powr = 1.0,powi = 0.0;
      for (j = 0; j < Nc; j++) {
	double termr = cr[j]*powr;
	double termi = cr[j]*powi;
	yr[i] += termr;
	yi[i] += termi;
	if (fabs(termr)+fabs(termi) < (fabs(yr[i])+fabs(yi[i]))*tol)
	  break;
	double powr_ = powr*xr[i]-powi*xi[i];
	double powi_ = powr*xi[i]+powi*xr[i];
	powr = powr_;
	powi = powi_;
      }
      warn = warn || j == Nc;
    }
  else if (xi == NULL)
    for (int i = 0; i < Nx; i++) {
      int j;
      double powr = 1.0;
      for (j = 0; j < Nc; j++) {
	double termr = cr[j]*powr;
	double termi = ci[j]*powr;
	yr[i] += termr;
	yi[i] += termi;
	if (fabs(termr)+fabs(termi) < (fabs(yr[i])+fabs(yi[i]))*tol)
	  break;
	powr *= xr[i];
      }
      warn = warn || j == Nc;
    }
  else
    for (int i = 0; i < Nx; i++) {
      int j;
      double powr = 1.0,powi = 0.0;
      for (j = 0; j < Nc; j++) {
	double termr = cr[j]*powr-ci[j]*powi;
	double termi = cr[j]*powi+ci[j]*powr;
	yr[i] += termr;
	yi[i] += termi;
	if (fabs(termr)+fabs(termi) < (fabs(yr[i])+fabs(yi[i]))*tol)
	  break;
	double powr_ = powr*xr[i]-powi*xi[i];
	double powi_ = powr*xi[i]+powi*xr[i];
	powr = powr_;
	powi = powi_;
      }
      warn = warn || j == Nc;
    }

  if (warn)
    mexWarnMsgIdAndTxt("powerseries:w1",
		       "Series did not converge to the prescribed "
		       "accuracy for all arguments.");
}
/*------------------------------------------------------------------------*/
