/* tsize.c */
/* S. Engblom 2005-04-10 */

#include <math.h>

#include "mex.h"
#include "matrix.h"

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check of syntax */
  if (nrhs != 1 && nrhs != 2 || nlhs > 1)
    mexErrMsgIdAndTxt("tsize:e1",
      "Expecting one or two inputs and one output.");

  /* input */
  mwSize ndimA = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *sizA = mxGetDimensions(prhs[0]);

  /* adjust the number of dimensions */
  if (ndimA == 2)
    ndimA -= (sizA[1] == 1)*(1+(sizA[0] == 1));

  /* full or indexed version */
  if (nrhs == 1) {
    /* output */
    plhs[0] = mxCreateDoubleMatrix(1,ndimA,mxREAL);
    double *siz = mxGetPr(plhs[0]);
    for (int i = 0; i < ndimA; i++)
      siz[i] = (double)sizA[i];
  }
  else {
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
        mxIsSparse(prhs[1]))
      mexErrMsgIdAndTxt("tsize:e2",
        "Dimension argument must be real, double and non-sparse.");

    /* input DIMS */
    const mwSize lenDIMS = mxGetNumberOfElements(prhs[1]);
    const double *prDIMS = mxGetPr(prhs[1]);

    /* check input DIMS */
    for (int i = 0; i < lenDIMS; i++)
      if (prDIMS[i] < 1.0 || prDIMS[i] != ceil(prDIMS[i]))
        mexErrMsgIdAndTxt("tsize:e3",
          "Size argument must be nonnegative integers.");

    /* output */
    plhs[0] = mxCreateDoubleMatrix(1,lenDIMS,mxREAL);
    double *siz = mxGetPr(plhs[0]);
    for (int i = 0; i < lenDIMS; i++) {
      if (prDIMS[i] <= ndimA)
        siz[i] = (double)sizA[(int)prDIMS[i]-1];
      else
        siz[i] = 1.0;
    }
  }
}
/*-----------------------------------------------------------------------*/
