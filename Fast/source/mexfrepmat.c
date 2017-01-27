/* mexfrepmat.c */

/* S. Engblom 2010-02-02 (Minor revision) */
/* S. Engblom 2007-05-04 (Revision) */
/* S. Engblom 2005-05-06 (Revision) */
/* S. Engblom 2004-10-22 */

#include <math.h>
#include <string.h>
#include <limits.h>

#include "mex.h"
#include "matrix.h"

// rmemcpy() below assumes this:
#if CHAR_BIT != 8
#error "Code assumes that 1 char is 1 byte."
#endif

/*------------------------------------------------------------------------*/

// forward declarations
void replicate(const char *prA,const char *piA,
	       int ndimA,const mwSize *sizA,int lenA,
	       char *prB,char *piB,
	       int ndimB,const mwSize *sizB,int lenB,int nbytes);
void sparse_replicate(const char *prA,const char *piA,
		      const mwIndex *irA,const mwSize *jcA,const mwSize *sizA,
		      char *prB,char *piB,
		      mwIndex *irB,mwSize *jcB,const mwSize *sizB,
		      int rM,int rN,int nnzB,int nbytes);
void rmemcpy(void *s,size_t block,size_t len);
void r2memcpy(void *s1,void *s2,size_t block,size_t len);

/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  // check of syntax
  if (nrhs != 2 || nlhs > 1)
    mexErrMsgIdAndTxt("frepmat:e1","Expecting two inputs and one output.");

  if (!(mxIsNumeric(prhs[0]) || mxIsChar(prhs[0]) || mxIsLogical(prhs[0])))
    mexErrMsgIdAndTxt("frepmat:e2","Only numerical, character and "
		      "logical arrays supported.");

  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]))
    mexErrMsgIdAndTxt("frepmat:e3","Size argument must be real, double "
		      "and non-sparse.");

  // input A
  const int ndimA = mxGetNumberOfDimensions(prhs[0]);
  const mwSize *sizA = mxGetDimensions(prhs[0]);
  const int lenA = mxGetNumberOfElements(prhs[0]);
  const char *prA = mxGetData(prhs[0]);
  const char *piA = mxGetImagData(prhs[0]);
  const bool complx = piA != NULL;
  const bool sparse = mxIsSparse(prhs[0]);

  // input REP
  const int lenREP = mxGetNumberOfElements(prhs[1]);
  const double *prREP = mxGetPr(prhs[1]);

  // output B
  int ndimB = ndimA >= lenREP ? ndimA : lenREP;
  mwSize sizB[ndimB];

  // check input REP
  for (int i = 0; i < lenREP; i++)
    if (prREP[i] < 0.0 || prREP[i] != ceil(prREP[i]))
      mexErrMsgIdAndTxt("frepmat:e4","Size argument must be nonnegative "
			"integers.");

  if (sparse && ndimB > 2)
    mexErrMsgIdAndTxt("frepmat:e5","Sparse N-dimensional matrices are "
		      "not supported.");
 
  // determine the size of the result B
  if (ndimB == ndimA) {
    int i;
    for (i = 0; i < lenREP; i++) sizB[i] = sizA[i]*prREP[i];
    for ( ; i < ndimA; i++) sizB[i] = sizA[i];
  }
  else {
    int i;
    for (i = 0; i < ndimA; i++) sizB[i] = prREP[i]*sizA[i];
    for ( ; i < lenREP; i++) sizB[i] = prREP[i];
  }

  // allocate and do the job
  if (!sparse) {
    plhs[0] = mxCreateNumericArray(ndimB,sizB,
				   mxGetClassID(prhs[0]),
				   complx ? mxCOMPLEX : mxREAL);
    replicate(prA,piA,ndimA,sizA,lenA,
	      mxGetData(plhs[0]),mxGetImagData(plhs[0]),
	      mxGetNumberOfDimensions(plhs[0]),
	      mxGetDimensions(plhs[0]),
	      mxGetNumberOfElements(plhs[0]),
	      mxGetElementSize(plhs[0]));
  }
  else {
    const int rM = lenREP < 1 ? 1 : prREP[0];
    const int rN = lenREP < 2 ? 1 : prREP[1];
    const int nnzB = mxGetJc(prhs[0])[sizA[1]]*rM*rN;

    if (mxGetClassID(prhs[0]) == mxLOGICAL_CLASS)
      plhs[0] = mxCreateSparseLogicalMatrix(sizB[0],sizB[1],nnzB);
    else
      plhs[0] = mxCreateSparse(sizB[0],sizB[1],nnzB,
			       complx ? mxCOMPLEX : mxREAL);

    sparse_replicate(prA,piA,mxGetIr(prhs[0]),mxGetJc(prhs[0]),sizA,
		     mxGetData(plhs[0]),mxGetImagData(plhs[0]),
		     mxGetIr(plhs[0]),mxGetJc(plhs[0]),
		     mxGetDimensions(plhs[0]),
		     rM,rN,nnzB,mxGetElementSize(plhs[0]));
  }
}
/*------------------------------------------------------------------------*/
void replicate(const char *prA,const char *piA,
	       int ndimA,const mwSize *sizA,int lenA,
	       char *prB,char *piB,
	       int ndimB,const mwSize *sizB,int lenB,int nbytes)
/* Replicates the matrix A within the matrix B by filling B in all its
   dimensions with repeated copies of A. Each element of both matrices
   is assumed to occupy nbytes bytes.

   Note: using char * instead of void * as done above means that it is
   implicitly assumed that sizeof(char) == 1 (byte). If this is not
   true, the code won't work. */
{
  // empty case
  if (lenB == 0) return;

  // various
  const bool complx = piA != NULL;
  int d1,dimcount[ndimA];
  int strB[ndimB+1],strBsizA[ndimA];

  // construct stride for B
  strB[0] = 1;
  for (int i = 0; i < ndimB; i++) strB[i+1] = strB[i]*sizB[i];

  // frequently used product
  for (int i = 0; i < ndimA; i++) strBsizA[i] = strB[i]*sizA[i];

  // find first nonmatching dimension
  for (d1 = 0; d1 < ndimA-1; d1++)
    if (sizA[d1] != sizB[d1]) break;

  // a single initial copy of A is placed in B
  memset(dimcount,0,ndimA*sizeof(dimcount[0]));
  for (int i = 0, j = 0, block = strBsizA[d1]; ; ) {
    // copy as much as possible
    memcpy(&prB[j*nbytes],&prA[i*nbytes],block*nbytes);
    if (complx) memcpy(&piB[j*nbytes],&piA[i*nbytes],block*nbytes);

    // increase pointer in A
    if ((i += block) == lenA) break;

    // increase pointer in B
    for (int k = d1+1; j += strB[k],++dimcount[k] == sizA[k]; k++) {
      dimcount[k] = 0;
      j -= strBsizA[k];
    }
  }

  // replicate over the first dimensions
  memset(dimcount,0,ndimA*sizeof(dimcount[0]));  
  for (int r = d1; r < ndimA-1; r++) {
    const int block = strBsizA[r];
    for (int i = 0; i < lenB; ) {
      // copy in a recursive fashion
      if (complx)
	r2memcpy(&prB[i*nbytes],&piB[i*nbytes],
		 block*nbytes,strB[r+1]*nbytes);
      else
	rmemcpy(&prB[i*nbytes],block*nbytes,strB[r+1]*nbytes);

      // determine next adress in B
      for (int k = r+1; i += strB[k],++dimcount[k] == sizA[k]; ) {
	dimcount[k] = 0;
	i -= strBsizA[k];
	if (++k == ndimA) {
	  i = lenB;
	  break;
	}
      }
    }
  }

  // the final dimensions are easier
  if (complx)
    r2memcpy(prB,piB,strBsizA[ndimA-1]*nbytes,lenB*nbytes);
  else
    rmemcpy(prB,strBsizA[ndimA-1]*nbytes,lenB*nbytes);
}
/*------------------------------------------------------------------------*/
void sparse_replicate(const char *prA,const char *piA,
		      const mwIndex *irA,const mwSize *jcA,
		      const mwSize *sizA,
		      char *prB,char *piB,
		      mwIndex *irB,mwSize *jcB,const mwSize *sizB,
		      int rM,int rN,int nnzB,int nbytes)
/* Sparse replication of matrix. Works as replicate() above (and with
   the same limitations on sizeof(char)) but for sparse matrices A and
   B instead. */
{
  // empty case
  if (nnzB == 0) return;

  const bool complx = piA != NULL;

  /* copy jcA to jcB and then replicate the copy rN-1 times in
     dimension 2 */
  memcpy(jcB,jcA,(sizA[1]+1)*sizeof(jcA[0]));
  for (int i = 1,jj = sizA[1]+1; i < rN; i++) {
    // add new offset for each round
    const int offset = jcB[jj-1];
    for (int ii = 1; ii <= sizA[1]; ii++,jj++)
      jcB[jj] = jcA[ii]+offset;
  }

  if (rM > 1) {
    /* replicate rM times in dimension 1 by scaling the cumulative
       pointer to account for added elements in each column */
    for (int jj = 1; jj <= sizB[1]; jj++) jcB[jj] *= rM;

    // construct irB columnwise and copy values at the same time
    for (int j = 0; j < sizA[1]; j++) {
      // copy and replicate the values
      memcpy(&prB[jcB[j]*nbytes],&prA[jcA[j]*nbytes],
	     (jcA[j+1]-jcA[j])*nbytes);
      if (complx) {
	memcpy(&piB[jcB[j]*nbytes],&piA[jcA[j]*nbytes],
	       (jcA[j+1]-jcA[j])*nbytes);
	r2memcpy(&prB[jcB[j]*nbytes],&piB[jcB[j]*nbytes],
		 (jcA[j+1]-jcA[j])*nbytes,(jcB[j+1]-jcB[j])*nbytes);
      }
      else
	rmemcpy(&prB[jcB[j]*nbytes],(jcA[j+1]-jcA[j])*nbytes,
		(jcB[j+1]-jcB[j])*nbytes);

      /* copy irA to irB and then replicate the copy rM-1 times in
	 dimension 1 */
      memcpy(&irB[jcB[j]],&irA[jcA[j]],(jcA[j+1]-jcA[j])*sizeof(irA[0]));
      for (int i = 1,jj = jcB[j]+jcA[j+1]-jcA[j]; i < rM; i++) {
	// add new offset for each round
	const int offset = i*sizA[0];
	for (int ii = jcA[j]; ii < jcA[j+1]; ii++,jj++)
	  irB[jj] = irA[ii]+offset;
      }
    }
  }
  else {
    // faster case, no replication in dimension 1
    memcpy(irB,irA,jcA[sizA[1]]*sizeof(irA[0]));
    memcpy(prB,prA,jcA[sizA[1]]*nbytes);
    if (complx)
      memcpy(piB,piA,jcA[sizA[1]]*nbytes);
  }

  // finally replicate the copy rN-1 times in dimension 2
  rmemcpy(irB,jcA[sizA[1]]*rM*sizeof(irB[0]),jcB[sizB[1]]*sizeof(irB[0]));
  if (complx)
    r2memcpy(prB,piB,jcA[sizA[1]]*rM*nbytes,
	     jcB[sizB[1]]*nbytes);
  else
    rmemcpy(prB,jcA[sizA[1]]*rM*nbytes,
	    jcB[sizB[1]]*nbytes);
}
/*------------------------------------------------------------------------*/
void rmemcpy(void *s,size_t block,size_t len)
/* Formally the same as memcpy(s+block,s,len-block), except that
   according to ANSI-C, memcpy() does not allow overlapping
   blocks. Assumes that sizeof(char) == 1 byte. */
{
  // the size of the block is doubled for each round
  while (block < len >> 1) {
    memcpy((char *)s+block,s,block);
    block <<= 1;
  }

  // the remainder
  memcpy((char *)s+block,s,len-block);
}
/*------------------------------------------------------------------------*/
void r2memcpy(void *s1,void *s2,size_t block,size_t len)
/* The same as rmemcpy(s1,block,len) followed by
   rmemcpy(s2,block,len). Used for real and imaginary parts. */
{
  while (block < len >> 1) {
    memcpy((char *)s1+block,s1,block);
    memcpy((char *)s2+block,s2,block);
    block <<= 1;
  }
  memcpy((char *)s1+block,s1,len-block);
  memcpy((char *)s2+block,s2,len-block);
}
/*------------------------------------------------------------------------*/
