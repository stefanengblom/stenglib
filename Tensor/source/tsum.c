/* tsum.c */
/* S. Engblom 2005-08-26 */

#include <string.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"

/* forward declarations */
void checkix(const double *ix,mwSize len,int *max);
void sortix(const double *ix,mwSize len,int *in);
void logicalsize(const mxArray *A,mwSize nda,mwSize *sizA);
void permute(mxArray **Ap,const mxArray *A,
	     const int *p1,const int *p2,int len);
void sum(int dim,double *prB,double *piB,
	 double *prA,double *piA,double **prC,double **piC,
	 mwSize ndims,mwSize *siz,mwSize *str);

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

#define MXTYPE(A,B) (mxIsComplex(A) || mxIsComplex(B) ? \
		     mxCOMPLEX : mxREAL)

#define ISDOUBLETENSOR(A)     (mxIsDouble(A) && !mxIsSparse(A))
#define ISDOUBLEREALTENSOR(A) (mxIsDouble(A) && !mxIsSparse(A) && \
			       !mxIsComplex(A))

/*-----------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* check of syntax */
  if (nrhs != 2 && nrhs != 4 || nlhs > 1)
    mexErrMsgIdAndTxt("tsum:e1","Expecting two or four inputs and "
		      "one output.");

  /* tensor-sum of two arrays */
  if (nrhs == 4) {
    /* check of input */
    if (!ISDOUBLETENSOR(prhs[0]) || !ISDOUBLETENSOR(prhs[1]))
      mexErrMsgIdAndTxt("tsum:e2","Expecting a double, "
			"non-sparse array.");
    if (!ISDOUBLEREALTENSOR(prhs[2]) || !ISDOUBLEREALTENSOR(prhs[3]))
      mexErrMsgIdAndTxt("tsum:e3","Index vector must be real, "
			"double and non-sparse.");

    /* inputs IA and IB in *unit offset* */
    const double *ia = mxGetPr(prhs[2])-1;
    const double *ib = mxGetPr(prhs[3])-1;
    const mwSize nda = mxGetNumberOfElements(prhs[2]);
    const mwSize ndb = mxGetNumberOfElements(prhs[3]);
    int max = 0; /* the number of dimensions of the output */

    /* perform some checks and determine the span of the indices */
    checkix(ia,nda,&max); checkix(ib,ndb,&max);

    /* sanity check */
    if (max > nda+ndb)
      mexErrMsgIdAndTxt("tsum:e5","Index vector must be contiguous.");

    /* Sort the indices out. The vectors ina and inb points back into ia
       and ib (in *unit offset*) according to the absolute ordering. */
    int ina[max],inb[max];

    sortix(ia,nda,memset(ina,0,max*sizeof(ina[0])));
    sortix(ib,ndb,memset(inb,0,max*sizeof(inb[0])));

    /* Set operations: the permutations contain indices pointing back
       into ia and ib in *unit offset* for the private parts (pa1 and
       pb1) and the common parts (pa2 and pb2). */
    int pa1[max+1],pa2[max+1],pb1[max+1],pb2[max+1];

    memset(pa1,0,(max+1)*sizeof(pa1[0])); memset(pa2,0,(max+1)*sizeof(pa2[0]));
    memset(pb1,0,(max+1)*sizeof(pb1[0])); memset(pb2,0,(max+1)*sizeof(pb2[0]));
    for (int i = 0,ia1 = 0,ib1 = 0,i2 = 0; i < max; i++) {

      /* a setdiff produces the private part */
      if (ina[i] && !inb[i])
	pa1[ia1++] = ina[i];
      /* setdiff... */
      else if (!ina[i] && inb[i])
	pb1[ib1++] = inb[i];
      /* an intersection produces the common part */
      else if (ina[i] && inb[i]) {
	pa2[i2] = ina[i]; pb2[i2++] = inb[i];
      }
      else
	mexErrMsgIdAndTxt("tsum:e5","Index vector must be contiguous.");
    }

    /* determine the logical size of the inputs A and B */
    mwSize sizA[nda],sizB[ndb];
    logicalsize(prhs[0],nda,sizA); logicalsize(prhs[1],ndb,sizB);

    /* permute to an order suitable for computations */
    mxArray *Ap,*Bp;

    permute(&Ap,prhs[0],pa2,pa1,nda);
    const double *prA = mxGetPr(Ap),*piA = mxGetPi(Ap);

    permute(&Bp,prhs[1],pb2,pb1,ndb);
    const double *prB = mxGetPr(Bp),*piB = mxGetPi(Bp);

    /* allocate the (permuted) output C (here prhs2[0] as it will be
       permuted by prhs2[1] in the end) */
    mxArray *prhs2[2];
    prhs2[1] = mxCreateDoubleMatrix(1,MAX(2,max),mxREAL);
    mwSize sizC[MAX(2,max)]; sizC[0] = sizC[1] = 1; /* singletons */
    int k = 1,m = 1,n = 1;
    {
      /* build the dimensions of C, including the sizes of the various
	 blocks along with the final permutation at the same time */
      int j = 0;
      double *p = mxGetPr(prhs2[1])-1; /* unit offset */
      p[1] = 1.0; p[2] = 2.0; /* singletons */

      for (int i = 0; pa2[i] != 0; i++) {
	/* check common sizes */
	if (sizA[pa2[i]-1] != sizB[pb2[i]-1])
	  mexErrMsgIdAndTxt("tsum:e8","Dimensions must agree.");
	k *= (sizC[j++] = sizA[pa2[i]-1]);
	p[(int)ia[pa2[i]]] = j;
      }
      for (int i = 0; pb1[i] != 0; i++) {
	m *= (sizC[j++] = sizB[pb1[i]-1]);
	p[(int)ib[pb1[i]]] = j;
      }
      for (int i = 0; pa1[i] != 0; i++) {
	n *= (sizC[j++] = sizA[pa1[i]-1]);
	p[(int)ia[pa1[i]]] = j;
      }
    }
    prhs2[0] = mxCreateNumericArray(MAX(2,max),sizC,mxDOUBLE_CLASS,
				    MXTYPE(prhs[0],prhs[1]));
    double *prC = mxGetPr(prhs2[0]),*piC = mxGetPi(prhs2[0]);

    /* Evaluate the sum on the canonical form C = tsum(A,B,[1 3],[1
       2]), where A is k-by-n, B is k-by-m and C is k-by-m-by-n. */

    /* real + real, complex + real, real + complex, complex + complex */
    if (piA == NULL || piB == NULL) {
      for (int ja = 0; ja < n; ja++, prA += k, prB -= k*m)
	for (int jb = 0; jb < m; jb++, prA -= k) {
	  for (int ii = 0; ii < k; ii++)
	    *prC++ = (*prA++)+(*prB++);
	}

      /* straightforward copy of imaginary parts, if any */
      if (piA != NULL)
	for (int ja = 0; ja < n; ja++, piA += k)
	  for (int jb = 0; jb < m; jb++, piC += k)
	    memcpy(piC,piA,k*sizeof(double));
      else if (piB != NULL)
	for (int ja = 0; ja < n; ja++, piC += k*m)
	  memcpy(piC,piB,k*m*sizeof(double));
    }
    else
      /* same loop as in the real case, rewritten for typographical
	 reasons */
      for (int ja = 0; ja < n; ja++) {
	for (int jb = 0; jb < m; jb++) {
	  for (int ii = 0; ii < k; ii++) {
	    *prC++ = (*prA++)+(*prB++);
	    *piC++ = (*piA++)+(*piB++);
	  }
	  prA -= k; piA -= k;
	}
	prA += k; piA += k;
	prB -= k*m; piB -= k*m;
      }

    /* deallocate permuted versions of A and B */
    mxDestroyArray(Bp); mxDestroyArray(Ap);

    /* permute according to input */
    mexCallMATLAB(1,plhs,2,prhs2,"permute");
    mxDestroyArray(prhs2[1]); mxDestroyArray(prhs2[0]);
  }
  /* tensor-sum of one array */
  else {
    /* check of input */
    if (!ISDOUBLETENSOR(prhs[0]))
      mexErrMsgIdAndTxt("tsum:e2","Expecting a double, "
			"non-sparse array.");
    if (!ISDOUBLEREALTENSOR(prhs[1]))
      mexErrMsgIdAndTxt("tsum:e3","Index vector must be real, "
			"double and non-sparse.");

    /* input IA */
    const double *ia = mxGetPr(prhs[1]);
    const mwSize nda = mxGetNumberOfElements(prhs[1]);

    for (int i = 0; i < nda; i++)
      if (ia[i] == 0.0 || ia[i] != ceil(ia[i]))
	mexErrMsgIdAndTxt("tsum:e9","Index vector must contain "
			  "nonzero integers.");

    if (mxGetNumberOfElements(prhs[0]) == 0) {
      /* empty dimensions is special since the output might be larger
	 than the input (a sum over an empty dimension produces a
	 singleton dimension) */
      const mwSize ndimA = mxGetNumberOfDimensions(prhs[0]);
      mwSize sizA[ndimA];
      memcpy(sizA,mxGetDimensions(prhs[0]),ndimA*sizeof(sizA[0]));

      for (int i = 0; i < nda; i++) {
	const int dim = abs((int)ia[i]);
	/* fix for a GCC-bug in abs() -- see sum() below */
	if (0 < dim && dim <= ndimA) sizA[dim-1] = 1;
      }
      plhs[0] = mxCreateNumericArray(ndimA,sizA,mxDOUBLE_CLASS,
				     mxIsComplex(prhs[1]) ? mxCOMPLEX 
				     : mxREAL);
      return;
    }
    else if (nda == 0) {
      plhs[0] = mxDuplicateArray(prhs[0]);
      return;
    }

    /* input and working arrays A and B */
    double *prA = mxGetPr(prhs[0]),*piA = mxGetPi(prhs[0]);
    double *prB,*piB = NULL;

    /* dimensions and stride, to be updated */
    const mwSize ndimA = mxGetNumberOfDimensions(prhs[0]);
    mwSize sizA[ndimA],strA[ndimA+1];
    memcpy(sizA,mxGetDimensions(prhs[0]),ndimA*sizeof(sizA[0]));
    strA[0] = 1;
    for (int i = 0; i < ndimA; i++) strA[i+1] = strA[i]*sizA[i];

    /* B is A summed along the first dimension in IA */
    {
      const int dim = abs((int)ia[0]);
      int remove = 1;
      /* fix for a GCC-bug in abs() -- see sum() below */
      if (0 < dim && dim <= ndimA) remove = sizA[dim-1];
      prB = mxMalloc(strA[ndimA]/remove*sizeof(double));
      if (piA != NULL) piB = mxMalloc(strA[ndimA]/remove*sizeof(double));
    }

    /* must be saved: */
    double *prB_ = prB,*piB_ = piB;

    /* sum over all given dimensions */
    sum((int)ia[0],prB,piB,prA,piA,NULL,NULL,ndimA,sizA,strA);
    for (int i = 1; i < nda; i++)
      sum((int)ia[i],NULL,NULL,prB,piB,&prB,&piB,ndimA,sizA,strA);

    /* account for negative indices and reallocate the final sum */
    memmove(prB_,prB,strA[ndimA]*sizeof(double)); prB = prB_;
    prB = mxRealloc(prB,strA[ndimA]*sizeof(double));    
    if (piB != NULL) {
      memmove(piB_,piB,strA[ndimA]*sizeof(double)); piB = piB_;
      piB = mxRealloc(piB,strA[ndimA]*sizeof(double));
    }
    mxArray *B = mxCreateDoubleMatrix(0,0,piB != NULL ? mxCOMPLEX 
				      : mxREAL);
    mxFree(mxGetPr(B)); mxSetPr(B,prB);
    if (piB != NULL) {
      mxFree(mxGetPi(B)); mxSetPi(B,piB);
    }
    mxSetDimensions(B,sizA,ndimA);
    plhs[0] = B;
  }
}
/*-----------------------------------------------------------------------*/
void checkix(const double *ix,mwSize len,int *max)
/* Checks the indices ix[1..len] and adjusts max accordingly. */
{
  for (int i = 1; i <= len; i++) {
    if (ix[i] <= 0.0 || ix[i] != ceil(ix[i]) || isinf(ix[i]))
      mexErrMsgIdAndTxt("tsum:e4","Index vector must contain "
			"positive integers.");

    if (*max < ix[i]) *max = (int)ix[i];
  }
}
/*-----------------------------------------------------------------------*/
void sortix(const double *ix,mwSize len,int *in)
/* Sorts the indices ix[1..len] and determines a rank-table, in, which
   is assumed to be allocated and cleared prior to call. */
{
  for (int i = 1; i <= len; i++) {
    const int ii = (int)ix[i]-1;

    if (in[ii] != 0)
      mexErrMsgIdAndTxt("tsum:e6","Indices must be distinct.");
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

  /* mask for Matlab's stupid convention */
  if (ndimA == 2) ndimA -= (msizA[1] == 1)*(1+(msizA[0] == 1));
  if (nda < ndimA)
    mexErrMsgIdAndTxt("tsum:e7","Wrong size of index vector.");

  memcpy(sizA,msizA,ndimA*sizeof(sizA[0]));
  for (int i = ndimA; i < nda; i++) sizA[i] = 1;
}
/*-----------------------------------------------------------------------*/
void permute(mxArray **Ap,const mxArray *A,
	     const int *p1,const int *p2,int len)
/* Permutes the array A according to the ordering p = [p1 p2] which
   must be of total length len. The result is stored in Ap. */
{
  /* input to Matlab */
  mxArray *prhs[2];
  prhs[0] = (mxArray *)A;
  prhs[1] = mxCreateDoubleMatrix(1,MAX(2,len),mxREAL);

  int j = 0;
  double *p = mxGetPr(prhs[1]);
  p[0] = 1.0; p[1] = 2.0; /* singletons */

  /* construct the permutation p */
  for (int i = 0; p1[i] != 0; i++) p[j++] = (double)p1[i];
  for (int i = 0; p2[i] != 0; i++) p[j++] = (double)p2[i];

  /* call Matlab and clean up */
  mexCallMATLAB(1,Ap,2,prhs,"permute");
  mxDestroyArray(prhs[1]);
}
/*-----------------------------------------------------------------------*/
void sum(int dim,double *prB,double *piB,
	 double *prA,double *piA,double **prC,double **piC,
	 mwSize ndims,mwSize *siz,mwSize *str)
/* Sums the array (prA,piA) along the dimension dim (dim must be
   nonzero and the corresponding dimension is not allowed to be of
   vanishing width). The resulting array is placed in (prB,piB) and
   contains on exit the newly created sum with siz[dim-1] = 1 (siz and
   str are thus updated). Use prB = piB = NULL to perform the sum
   within the original memory (prA,piA). In this case, the result
   starts at (prC,piC) which is somewhere inside (prA,piA). The input
   dim may be negative, indicating summation backwards. */
{
  const bool fwd = dim >= 0;
  dim = abs(dim);

  /* early returns (plus a small fix for GCC-problems with abs(the
     most negative number)) */
  if (dim < 0 || ndims < dim || siz[dim-1] <= 1) {
    if (prB == NULL) {
      *prC = prA;
      *piC = piA;
    }
    else {
      memcpy(prB,prA,str[ndims]*sizeof(double));
      if (piA != NULL) memcpy(piB,piA,str[ndims]*sizeof(double));
    }
    return;
  }

  if (piA == NULL) {
    if (fwd) {
      /* forward summation */
      if (prB == NULL) *prC = prB = prA;

      /* summing along the first non-singleton dimension is faster */
      if (str[dim-1] == 1)
	for (int i = 0; i < str[ndims]; i += str[dim], prB++) {
	  *prB = *prA++;
	  for (int j = 1; j < siz[dim-1]; j++)
	    *prB += *prA++;
	}
      else
	for (int i = 0; i < str[ndims]; i += str[dim], prB += str[dim-1]) {
	  memmove(prB,prA,str[dim-1]*sizeof(double));
	  prA += str[dim-1];
	  for (int j = 1; j < siz[dim-1]; j++)
	    for (int k = 0; k < str[dim-1]; k++)
	      prB[k] += *prA++;
	}
    }
    else {
      /* backward summation */
      const bool inside = prB == NULL;
      if (inside)
	prB = prA += str[ndims];
      else {
	prB += str[ndims]/siz[dim-1];
	prA += str[ndims];
      }

      /* it's still easier to sum along the first non-singleton... */
      if (str[dim-1] == 1) {
	for (int i = 0; i < str[ndims]; i += str[dim]) {
	  *--prB = *--prA;
	  for (int j = 1; j < siz[dim-1]; j++)
	    *prB += *--prA;
	}
	/* the sum has been computed in the wrong memory area */
	if (inside) *prC = prB;
      }
      else {
	for (int i = 0; i < str[ndims]; i += str[dim]) {
	  prB -= str[dim-1];
	  prA -= str[dim-1];
	  memmove(prB,prA,str[dim-1]*sizeof(double));
	  for (int j = 1; j < siz[dim-1]; j++)
	    for (int k = str[dim-1]; --k >= 0; )
	      prB[k] += *--prA;
	}
	if (inside) *prC = prB;
      }
    }
  }
  else {
    /* the complex case follows analogously */
    if (fwd) {
      if (prB == NULL) {
	*prC = prB = prA;
	*piC = piB = piA;
      }

      if (str[dim-1] == 1)
	for (int i = 0; i < str[ndims]; i += str[dim], prB++, piB++) {
	  *prB = *prA++; *piB = *piA++;
	  for (int j = 1; j < siz[dim-1]; j++) {
	    *prB += *prA++;
	    *piB += *piA++;
	  }
	}
      else
	for (int i = 0; i < str[ndims]; i += str[dim], 
	       prB += str[dim-1], piB += str[dim-1]) {
	  memmove(prB,prA,str[dim-1]*sizeof(double));
	  prA += str[dim-1];
	  memmove(piB,piA,str[dim-1]*sizeof(double));
	  piA += str[dim-1];
	  for (int j = 1; j < siz[dim-1]; j++)
	    for (int k = 0; k < str[dim-1]; k++) {
	      prB[k] += *prA++;
	      piB[k] += *piA++;
	    }
	}
    }
    else {
      const bool inside = prB == NULL;
      if (inside) {
	prB = prA += str[ndims];
	piB = piA += str[ndims];
      }
      else {
	prB += str[ndims]/siz[dim-1];
	prA += str[ndims];
	piB += str[ndims]/siz[dim-1];
	piA += str[ndims];
      }

      if (str[dim-1] == 1) {
	for (int i = 0; i < str[ndims]; i += str[dim]) {
	  *--prB = *--prA;
	  *--piB = *--piA;
	  for (int j = 1; j < siz[dim-1]; j++) {
	    *prB += *--prA;
	    *piB += *--piA;
	  }
	}
	if (inside) {
	  *prC = prB;
	  *piC = piB;
	}
      }
      else {
	for (int i = 0; i < str[ndims]; i += str[dim]) {
	  prB -= str[dim-1];
	  prA -= str[dim-1];
	  memmove(prB,prA,str[dim-1]*sizeof(double));
	  piB -= str[dim-1];
	  piA -= str[dim-1];
	  memmove(piB,piA,str[dim-1]*sizeof(double));
	  for (int j = 1; j < siz[dim-1]; j++)
	    for (int k = str[dim-1]; --k >= 0; ) {
	      prB[k] += *--prA;
	      piB[k] += *--piA;
	    }
	}
	if (inside) {
	  *prC = prB;
	  *piC = piB;
	}
      }
    }
  }

  /* update the dimensions */
  siz[dim-1] = 1;
  for (int i = dim-1; i < ndims; i++) str[i+1] = str[i]*siz[i];
}
/*-----------------------------------------------------------------------*/
