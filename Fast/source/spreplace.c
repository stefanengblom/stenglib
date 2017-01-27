/* spreplace.c */
/* S. Engblom 2016-11-23 */

#include <string.h>

#include "mex.h"
#include "matrix.h"

#define MXTYPE(A) (mxIsComplex(A) ? mxCOMPLEX : mxREAL)

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check of syntax */
  if (nrhs != 2 || nlhs > 1)
    mexErrMsgTxt("Expecting two inputs and one output.");

  /* input */
  if (!mxIsSparse(prhs[0]))
    mexErrMsgTxt("Expecting a sparse array as the first input.");
  const int M = mxGetM(prhs[0]),N = mxGetN(prhs[0]);
  const mwSize *jc = mxGetJc(prhs[0]);
  const mwIndex *ir = mxGetIr(prhs[0]);

  if (!mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1]) != jc[N])
    mexErrMsgTxt("Expecting a matching double array as "
		 "the second input.");
  const double *prx = mxGetPr(prhs[1]),*pix = mxGetPi(prhs[1]);

  /* allocate and create output */
  plhs[0] = mxCreateSparse(M,N,jc[N],MXTYPE(prhs[1]));
  memcpy(mxGetJc(plhs[0]),jc,(N+1)*sizeof(jc[0]));
  memcpy(mxGetIr(plhs[0]),ir,jc[N]*sizeof(ir[0]));
  memcpy(mxGetPr(plhs[0]),prx,jc[N]*sizeof(double));
  if (pix)
    memcpy(mxGetPi(plhs[0]),pix,jc[N]*sizeof(double));
}
/*-----------------------------------------------------------------------*/
