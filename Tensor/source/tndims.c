/* tndims.c */
/* S. Engblom 2005-04-10 */

#include "mex.h"
#include "matrix.h"

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check of syntax */
  if (nrhs != 1 || nlhs > 1)
    mexErrMsgIdAndTxt("tndims:e1","Expecting one input and one output.");

  /* input */
  mwSize ndimA = mxGetNumberOfDimensions(prhs[0]);

  /* adjust */
  if (ndimA == 2)
    ndimA -= (mxGetN(prhs[0]) == 1)*(1+(mxGetM(prhs[0]) == 1));

  /* output */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(plhs[0]) = (double)ndimA;
}
/*-----------------------------------------------------------------------*/
