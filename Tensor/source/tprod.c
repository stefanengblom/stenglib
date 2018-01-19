/* tprod.c */

/* S. Engblom 2010-09-23 (Minor revision) */
/* S. Engblom 2010-02-02 (Minor revision) */
/* S. Engblom 2005-04-29 */

#include <string.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"

/* Note: All integer arguments in calls to the BLAS function dgemm()
   will be casted according to (int *)&BLASINT. This is to provide for
   the (bad!) platforms where sizeof(int) < sizeof(size_t) and the
   binary library uses size_t while the declaration is erroneously
   int. For such cases you will have to manually define BLASINT as
   size_t when calling the compiler. */
#ifndef BLASINT
#define BLASINT int
#endif

// declaration of BLAS-kernel
void dgemm_(const char *,const char *,
	    const int *,const int *,const int *,const double *,
	    const double *,const int *,
	    const double *,const int *,
	    const double *,double *,const int *);

// forward declarations
void checkix(const double *ix,mwSize len,int *min,int *max);
void sortix(const double *ix,mwSize len,int edg,int *in);
void logicalsize(const mxArray *A,mwSize nda,mwSize *sizA);
void permute(const mxArray *A,
	     const int *p1,const int*p2,const int *p3,int len,
	     mxArray **Ap);

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  // check of syntax
  if (nrhs != 4 || nlhs > 1)
    mexErrMsgIdAndTxt("tprod:e9",
		      "Expecting four inputs and one output.");

  // check of input
  if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]) ||
      !mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]))
    mexErrMsgIdAndTxt("tprod:e8",
		      "Expecting a double, non-sparse array.");
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || 
      mxIsSparse(prhs[2]) || !mxIsDouble(prhs[3]) || 
      mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]))
    mexErrMsgIdAndTxt("tprod:e10",
		      "Index vector must be real, double and non-sparse.");

  // inputs IA and IB in *unit offset*
  const double *ia = mxGetPr(prhs[2])-1;
  const double *ib = mxGetPr(prhs[3])-1;
  const mwSize nda = mxGetNumberOfElements(prhs[2]);
  const mwSize ndb = mxGetNumberOfElements(prhs[3]);
  int min = 0,max = 0;

  // perform some checks and determine the span of the indices
  checkix(ia,nda,&min,&max);
  checkix(ib,ndb,&min,&max);

  // sanity check
  if (max-min > nda+ndb)
    mexErrMsgIdAndTxt("tprod:e7","Index vector is not contiguous or "
                      "a summation index is missing.");

  /* Sort the indices out. The vectors ina and inb points back into ia
     and ib (in *unit offset*) according to the absolute ordering. */
  const int edg = min;
  const int len = max-min+1;
  const int ndc = max; // the number of dimensions of the output
  int ina[len],inb[len];

  sortix(ia,nda,edg,memset(ina,0,len*sizeof(ina[0])));
  sortix(ib,ndb,edg,memset(inb,0,len*sizeof(inb[0])));

  /* Set operations: the permutations contain indices pointing back
     into ia and ib in *unit offset* for the private parts (pa1 and
     pb1), the negative common parts (pa2 and pb2) and the positive
     common parts (pa3 and pb3). */
  int pa1[len],pa2[len],pa3[len];
  int pb1[len],pb2[len],pb3[len];

  memset(pa1,0,len*sizeof(pa1[0])); memset(pb1,0,len*sizeof(pb1[0]));
  memset(pa2,0,len*sizeof(pa2[0])); memset(pb2,0,len*sizeof(pb2[0]));
  memset(pa3,0,len*sizeof(pa3[0])); memset(pb3,0,len*sizeof(pb3[0]));
  for (int i = 0,ia1 = 0,ib1 = 0,i2 = 0,i3 = 0; i < len; i++) {
    const int ix = i+edg;

    // a setdiff produces the private part
    if (ina[i] && !inb[i]) {
      if (ix < 0)
        mexErrMsgIdAndTxt("tprod:e3","Summation index is missing.");
      pa1[ia1++] = ina[i];
    }
    // setdiff...
    else if (!ina[i] && inb[i]) {
      if (ix < 0)
        mexErrMsgIdAndTxt("tprod:e3","Summation index is missing.");
      pb1[ib1++] = inb[i];
    }
    // an intersection produces the (positive and negative) common part
    else if (ina[i] && inb[i]) {
      if (ix < 0) {
        pa2[i2] = ina[i]; pb2[i2++] = inb[i];
      }
      else {
        pa3[i3] = ina[i]; pb3[i3++] = inb[i];
      }
    }
    else if (ix != 0)
      mexErrMsgIdAndTxt("tprod:e4","Index vector must be contiguous.");
  }

  // determine the logical size of the inputs A and B
  mwSize sizA[nda],sizB[ndb];
  logicalsize(prhs[0],nda,sizA);
  logicalsize(prhs[1],ndb,sizB);

  // check common sizes
  for (int i = 0; pa2[i] != 0; i++)
    if (sizA[pa2[i]-1] != sizB[pb2[i]-1])
      mexErrMsgIdAndTxt("tprod:e6","Dimensions must agree.");
  for (int i = 0; pa3[i] != 0; i++)
    if (sizA[pa3[i]-1] != sizB[pb3[i]-1])
      mexErrMsgIdAndTxt("tprod:e6","Dimensions must agree.");

  // permute to an order suitable for computations
  mxArray *Ap,*Bp;

  permute(prhs[0],pa1,pa2,pa3,nda,&Ap);
  const double *prA = mxGetPr(Ap);
  const double *piA = mxGetPi(Ap);

  permute(prhs[1],pb1,pb2,pb3,ndb,&Bp);
  const double *prB = mxGetPr(Bp);
  const double *piB = mxGetPi(Bp);

  // allocate the (permuted) output C
  mwSize sizC[ndc > 2 ? ndc : 2];
  sizC[0] = sizC[1] = 1; /* singletons */
  BLASINT m = 1,n = 1,k = 1,l = 1;
  {
    // the size of C, including the sizes of the various blocks
    int j = 0;
    for (int i = 0; pa1[i] != 0; i++) m *= (sizC[j++] = sizA[pa1[i]-1]);
    for (int i = 0; pa2[i] != 0; i++) k *= sizA[pa2[i]-1];
    for (int i = 0; pb1[i] != 0; i++) n *= (sizC[j++] = sizB[pb1[i]-1]);
    for (int i = 0; pa3[i] != 0; i++) l *= (sizC[j++] = sizA[pa3[i]-1]);
  }
  mxArray *Cp = mxCreateNumericArray(ndc > 2 ? ndc : 2,sizC,
				     mxDOUBLE_CLASS,
				     piA != NULL || piB != NULL ? 
				     mxCOMPLEX : mxREAL);
  double *prC = mxGetPr(Cp);
  double *piC = mxGetPi(Cp);

  /* Evaluate the product as a series of matrix products, C(:,:,i) =
     A(:,:,i).'*B(:,:,i), i = 1..l. This corresponds to the canonical
     form C = tprod(A,B,[-1 1 3],[-1 2 3]) where A is k-by-m-by-l, B
     is k-by-n-by-l and C is m-by-n-by-l. */
  const double one = 1.0,minusone = -1.0;
  const char *t1 = "T",*t2 = "N";

  // real x real, complex x real, real x complex, complex x complex
  if (piC == NULL)
    for (int jc = 0,ja = 0,jb = 0; jc < m*n*l; ) {
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&prC[jc],(int *)&m);
      jc += m*n; ja += m*k; jb += k*n;
    }
  else if (piB == NULL)
    for (int jc = 0,ja = 0,jb = 0; jc < m*n*l; ) {
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&prC[jc],(int *)&m);
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
             &piA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&piC[jc],(int *)&m);
      jc += m*n; ja += m*k; jb += k*n;
    }
  else if (piA == NULL)
    for (int jc = 0,ja = 0,jb = 0; jc < m*n*l; ) {
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&prC[jc],(int *)&m);
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&piB[jb],(int *)&k,&one,&piC[jc],(int *)&m);
      jc += m*n; ja += m*k; jb += k*n;
    }
  else
    for (int jc = 0,ja = 0,jb = 0; jc < m*n*l; ) {
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&prC[jc],(int *)&m);
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
             &piA[ja],(int *)&k,&prB[jb],(int *)&k,&one,&piC[jc],(int *)&m);
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&one,
	     &prA[ja],(int *)&k,&piB[jb],(int *)&k,&one,&piC[jc],(int *)&m);
      dgemm_(t1,t2,(int *)&m,(int *)&n,(int *)&k,&minusone,
	     &piA[ja],(int *)&k,&piB[jb],(int *)&k,&one,&prC[jc],(int *)&m);
      jc += m*n; ja += m*k; jb += k*n;
    }

  // can deallocate permuted versions of A and B now
  mxDestroyArray(Bp);
  mxDestroyArray(Ap);

  // permute according to input
  mxArray *prhs2[2];
  prhs2[0] = Cp;
  prhs2[1] = mxCreateDoubleMatrix(1,ndc > 2 ? ndc : 2,mxREAL);
  {
    // note: unit offset
    int j = 1;
    double *p = mxGetPr(prhs2[1])-1;
    p[1] = 1.0; p[2] = 2.0; // singletons 

    for (int i = 0; pa1[i] != 0; i++) p[(int)ia[pa1[i]]] = j++;
    for (int i = 0; pb1[i] != 0; i++) p[(int)ib[pb1[i]]] = j++;
    for (int i = 0; pa3[i] != 0; i++) p[(int)ia[pa3[i]]] = j++;
  }
  mexCallMATLAB(1,plhs,2,prhs2,"permute");
  mxDestroyArray(prhs2[1]);
  mxDestroyArray(Cp);
}
/*-----------------------------------------------------------------------*/
void checkix(const double *ix,mwSize len,int *min,int *max)
/* Checks the indices ix[1..len] and adjusts min and max
   accordingly. */
{
  for (int i = 1; i <= len; i++) {
    if (ix[i] == 0.0 || ix[i] != ceil(ix[i]) || isinf(ix[i]))
      mexErrMsgIdAndTxt("tprod:e1",
			"Index vector must contain nonzero integers.");

    if (ix[i] < *min) *min = (int)ix[i];
    else if (*max < ix[i]) *max = (int)ix[i];
  }
}
/*-----------------------------------------------------------------------*/
void sortix(const double *ix,mwSize len,int edg,int *in)
/* Sorts the indices ix[1..len] and determines a rank-table in. The
   smallest index is edg and in is assumed to be allocated and cleared
   prior to call. */
{
  for (int i = 1; i <= len; i++) {
    const int ii =  (int)ix[i]-edg;

    if (in[ii] != 0)
      mexErrMsgIdAndTxt("tprod:e2","Indices must be distinct.");
    in[ii] = i;
  }
}
/*-----------------------------------------------------------------------*/
void logicalsize(const mxArray *A,mwSize nda,mwSize *sizA)
/* Determines the size sizA of the array A. The sizes of all
   dimensions i, 1 <= i <= nda are included in sizA. */
{
  mwSize ndimA = mxGetNumberOfDimensions(A);
  const mwSize *msizA = mxGetDimensions(A);

  // mask for Matlab's stupid convention
  if (ndimA == 2) ndimA -= (msizA[1] == 1)*(1+(msizA[0] == 1));
  if (nda < ndimA)
    mexErrMsgIdAndTxt("tprod:e5","Wrong size of index vector.");

  memcpy(sizA,msizA,ndimA*sizeof(sizA[0]));
  for (int i = ndimA; i < nda; i++) sizA[i] = 1;
}
/*-----------------------------------------------------------------------*/
void permute(const mxArray *A,
	     const int *p1,const int*p2,const int *p3,int len,
	     mxArray **Ap)
/* Permutes the array A according to the ordering p = [p2 p1 p3] which
   must be of total length len. The result is stored in Ap. */
{
  // input to Matlab
  mxArray *prhs[2];
  prhs[0] = (mxArray *)A;
  prhs[1] = mxCreateDoubleMatrix(1,len > 2 ? len : 2,mxREAL);

  int j = 0;
  double *p = mxGetPr(prhs[1]);
  p[0] = 1.0; p[1] = 2.0; // singletons

  // construct the permutation p
  for (int i = 0; p2[i] != 0; i++) p[j++] = (double)p2[i];
  for (int i = 0; p1[i] != 0; i++) p[j++] = (double)p1[i];
  for (int i = 0; p3[i] != 0; i++) p[j++] = (double)p3[i];

  // call Matlab and clean up
  mexCallMATLAB(1,Ap,2,prhs,"permute");
  mxDestroyArray(prhs[1]);
}
/*-----------------------------------------------------------------------*/
