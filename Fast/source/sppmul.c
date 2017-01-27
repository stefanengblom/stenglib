/* sppmul.c */
/* S. Engblom 2005-10-27 */

#include "mex.h"
#include "matrix.h"

#define MXTYPE(A,B) (mxIsComplex(A) || mxIsComplex(B) ? \
		     mxCOMPLEX : mxREAL)

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check of syntax */
  if (nrhs != 3 || nlhs > 1)
    mexErrMsgTxt("Expecting three inputs and one output.");

  /* input */
  if (!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Expecting a sparse array as the first input.");
  const int M = mxGetM(prhs[0]),N = mxGetN(prhs[0]);
  const mwSize *jc = mxGetJc(prhs[0]);
  const mwIndex *ir = mxGetIr(prhs[0]);

  if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != jc[N])
    mexErrMsgTxt("Expecting a matching double array as "
		 "the second input.");
  const double *prv = mxGetPr(prhs[1]),*piv = mxGetPi(prhs[1]);

  if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) > 2 || 
      mxGetM(prhs[2]) != N)
    mexErrMsgTxt("Expecting a double matrix with matching first "
		 "dimension as the third input.");
  const int K = mxGetN(prhs[2]);
  const double *prX = mxGetPr(prhs[2]),*piX = mxGetPi(prhs[2]);

  /* allocate output */
  plhs[0] = mxCreateDoubleMatrix(M,K,MXTYPE(prhs[1],prhs[2]));
  double *prY = mxGetPr(plhs[0]),*piY = mxGetPi(plhs[0]);

  /* evaluate the product */
  if (piY == NULL)
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < N; i++)
	for (mwSize j = jc[i]; j < jc[i+1]; j++)
	  prY[ir[j]] += prv[j]*prX[i];
      prY += M;
      prX += N;
    }
  else if (piX == NULL)
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < N; i++)
	for (mwSize j = jc[i]; j < jc[i+1]; j++) {
	  prY[ir[j]] += prv[j]*prX[i];
	  piY[ir[j]] += piv[j]*prX[i];
	}
      prY += M; piY += M;
      prX += N;
    }
  else if (piv == NULL)
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < N; i++)
	for (mwSize j = jc[i]; j < jc[i+1]; j++) {
	  prY[ir[j]] += prv[j]*prX[i];
	  piY[ir[j]] += prv[j]*piX[i];
	}
      prY += M; piY += M;
      prX += N; piX += N;
    }
  else
    for (int k = 0; k < K; k++) {
      for (int i = 0; i < N; i++)
	for (mwSize j = jc[i]; j < jc[i+1]; j++) {
	  prY[ir[j]] += prv[j]*prX[i]-piv[j]*piX[i];
	  piY[ir[j]] += piv[j]*prX[i]+piv[j]*prX[i];
	}
      prY += M; piY += M;
      prX += N; piX += N;
    }
}
/*-----------------------------------------------------------------------*/
