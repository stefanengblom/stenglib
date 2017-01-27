/* fsparse.c */

/* S. Engblom 2013-12-02 (OpenMP) */
/* S. Engblom 2010-02-02 (Minor revision) */
/* S. Engblom 2007-05-04 (Revision) */
/* S. Engblom 2005-05-05 (Revision) */
/* S. Engblom 2004-10-29 */

#include <math.h>
#include <string.h>
// temporary fix for CC under Solaris:
#ifndef NO_STDINT
#include <stdint.h>
#endif

#include "mex.h"
#include "matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef FSPARSE_TIME

#define StartTime
#define StopTime
#define GetTime(dest)

#else

#include <sys/time.h>

static double *time_vec; // global dummy for output time vector

static struct timeval TIME_before,TIME_after;
#define StartTime gettimeofday(&TIME_before,NULL)
#define StopTime gettimeofday(&TIME_after,NULL)
#define GetTime(dest) ((dest) = (double)(TIME_after.tv_sec-TIME_before.tv_sec)+ \
		       (TIME_after.tv_usec-TIME_before.tv_usec)/1000000.0)

#endif // FSPARSE_TIME

// print intermediate results:
#undef PRINT_INTERMEDIATE

/*------------------------------------------------------------------------*/

// forward declarations
bool mx_IsInt(const mxArray *array_ptr);

bool getix(int **ix,int M,int N,int *max,bool nocopy,const mxArray *IX);

mxArray *sparse2sparse(const mxArray *S);
mxArray *full2sparse(const mxArray *S);

void squeeze(mxArray *S);

void sparse_insert(mwIndex *irS,double *prS,double *piS,
		   const int *irank,const int *rank,const mwSize *jrS,
		   const int *ii,
		   const double *sr,const double *si,
		   int smod,int sdiv,int len,int M);
void sparse_inserti(mwIndex *irS,double *prS,double *piS,
		    const int *irank,
		    const int *ii,int imod,
		    const double *sr,const double *si,
		    int smod,int sdiv,int len);

mxArray *sparse(const int *ii,const int *jj,
		const double *sr,const double *si,
		int smod,int sdiv,
		int len,int M,int N,int Nzmax);
mxArray *sparse_nosort(const int *ii,const int *jj,
		       const double *sr,const double *si,
		       int smod,int sdiv,
		       int len,int M,int N,int Nzmax);
mxArray *gsparse(const int *ii,int imod,
		 const int *jj,int jdiv,
		 const double *sr,const double *si,
		 int smod,int sdiv,
		 int len,int M,int N,int Nzmax);
mxArray *gsparse_nosort(const int *ii,int imod,
			const int *jj,int jdiv,
			const double *sr,const double *si,
			int smod,int sdiv,
			int len,int M,int N,int Nzmax);

/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
#ifdef FSPARSE_TIME
  // output time vector allocated here already
  plhs[1] = mxCreateDoubleMatrix(1,6,mxREAL);
  time_vec = mxGetPr(plhs[1]);

  // check of syntax
  if (nrhs < 1 || nrhs == 2 || 6 < nrhs  || nlhs > 2)
    mexErrMsgIdAndTxt("fsparse:e1",
		      "Expecting 1, 3..6 inputs and one or two outputs.");
#else
  if (nrhs < 1 || nrhs == 2 || 6 < nrhs  || nlhs > 1)
    mexErrMsgIdAndTxt("fsparse:e1",
		      "Expecting 1, 3..6 inputs and one output.");
#endif

  // special case for one input
  if (nrhs == 1) {
    if (!mxIsDouble(prhs[0]))
      mexErrMsgIdAndTxt("fsparse:e2",
			"Single input argument must be double.");

    if (mxIsSparse(prhs[0])) {
      plhs[0] = sparse2sparse(prhs[0]);
      return;
    }
    else {
      if (mxGetNumberOfDimensions(prhs[0]) > 2)
	mexErrMsgIdAndTxt("fsparse:e3",
			  "Single input argument must be 2-D.");
      plhs[0] = full2sparse(prhs[0]);
      return;
    }
  }

  if (!mx_IsInt(prhs[0]) && !mxIsDouble(prhs[0]) || 
      mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ||
      !mx_IsInt(prhs[1]) && !mxIsDouble(prhs[1]) || 
      mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]))
    mexErrMsgIdAndTxt("fsparse:e4",
		      "Index argument must be real, double or integers " 
		      "and non-sparse.");

  if (!mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]))
    mexErrMsgIdAndTxt("fsparse:e5",
		      "Value argument must be double and non-sparse.");

  if (mxGetNumberOfDimensions(prhs[0]) > 2 ||
      mxGetNumberOfDimensions(prhs[1]) > 2 ||
      mxGetNumberOfDimensions(prhs[2]) > 2)
    mexErrMsgIdAndTxt("fsparse:e6","Input arguments must be 2-D.");

  if (nrhs > 3) {
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]))
      mexErrMsgIdAndTxt("fsparse:e7","Size argument must be real, "
			"double and non-sparse.");
    if (nrhs > 4) {
      if (mxGetNumberOfElements(prhs[4]) != 0 && !mxIsChar(prhs[4]))
	mexErrMsgIdAndTxt("fsparse:e8","Flag argument must be "
			  "character array.");
      if (nrhs > 5) {
	// must determine the number of threads here (before any parallel region)
	if (!mxIsDouble(prhs[5]) || 
	    mxGetNumberOfElements(prhs[5]) != 1)
	  mexErrMsgIdAndTxt("fsparse:e17","Number of threads must be a "
			    "real scalar.");
	const int nthreads = (int)mxGetScalar(prhs[5]);
	if (nthreads < 1)
	  mexErrMsgIdAndTxt("fsparse:e18","Number of threads must be >= 1.");
#ifdef _OPENMP
	omp_set_num_threads(nthreads);
#endif // ignored otherwise
      }
    }
  }

  // input
  const int Mii = mxGetM(prhs[0]),Nii = mxGetN(prhs[0]);
  const int Mjj = mxGetM(prhs[1]),Njj = mxGetN(prhs[1]);
  const bool nocopyii = mx_IsInt(prhs[0]),nocopyjj = mx_IsInt(prhs[1]);  
  int *ii,*jj;
  const int Mss = mxGetM(prhs[2]),Nss = mxGetN(prhs[2]);
  const double *sr = mxGetPr(prhs[2]);
  const double *si = mxGetPi(prhs[2]);
  const int len = Mii*Njj;
  int M = 0,N = 0,Nzmax = -1,sort = 1;

  // check of 'assembly' syntax
  if (Nii != Njj && Nii != 1 ||
      Mjj != Mii && Mjj != 1 ||
      Mss != Mii && Mss != 1 ||
      Nss != Njj && Nss != 1)
    mexErrMsgIdAndTxt("fsparse:e9","Sizes mismatch.");

  // input ii and jj
  StartTime;
  bool ok1,ok2;
#ifndef _OPENMP
  ok1 = getix(&ii,Mii,Nii,&M,nocopyii,prhs[0]);
  ok2 = getix(&jj,Mjj,Njj,&N,nocopyjj,prhs[1]);
#else
  // independent calls
#pragma omp single nowait
  ok1 = getix(&ii,Mii,Nii,&M,nocopyii,prhs[0]);
#pragma omp single nowait
  ok2 = getix(&jj,Mjj,Njj,&N,nocopyjj,prhs[1]);
#pragma omp barrier
#endif // _OPENMP
  if (!ok1 || !ok2)
    mexErrMsgIdAndTxt("fsparse:e10","Index argument must be "
		      "nonnegative integers.");
  StopTime;
  GetTime(time_vec[0]);

  // determine the input dimensions [M N Nzmax] of the output
  if (nrhs > 3) {
    const int szlen = mxGetNumberOfElements(prhs[3]);
    const double *szval = mxGetPr(prhs[3]);

    if (szlen > 0) {
      if (szval[0] < 0.0 || szval[0] != ceil(szval[0]))
	mexErrMsgIdAndTxt("fsparse:e11","Size argument must be "
			  "nonnegative integer.");
      if (M > szval[0])
	mexErrMsgIdAndTxt("fsparse:e12","Index exceeds matrix dimensions.");
      M = szval[0];
      if (szlen > 1) {
	if (szval[1] < 0.0 || szval[1] != ceil(szval[1]))
	  mexErrMsgIdAndTxt("fsparse:e11","Size argument must be "
			    "nonnegative integer.");
	if (N > szval[1])
	  mexErrMsgIdAndTxt("fsparse:e12","Index exceeds "
			    "matrix dimensions.");
	N = szval[1];
	if (szlen > 2) {
	  if (szval[2] < 0.0 || szval[2] != ceil(szval[2]))
	    mexErrMsgIdAndTxt("fsparse:e13","Nzmax argument must be "
			      "nonnegative integer.");
	  Nzmax = szval[2];
	  if (szlen > 3)
	    mexErrMsgIdAndTxt("fsparse:e14","Size argument must "
			      "contain 3 elements or less.");
	}
      }
    }
  }

  // sorted/not sorted output
  if (nrhs > 4) {
    if (mxGetNumberOfElements(prhs[4]) == 0)
      sort = 1;
    else {
      char buf[15]; // read 14 characters at most
      if (mxGetString(prhs[4],buf,15) != 0)
	mexErrMsgIdAndTxt("fsparse:e15","Unrecognized flag argument.");
      if (strcmp(buf,"nosort") == 0)
	sort = 0;
      else if (strcmp(buf,"sort") != 0)
	mexErrMsgIdAndTxt("fsparse:e15","Unrecognized flag argument.");
    }
  }

  // empty case
  if (len == 0) {
    plhs[0] = mxCreateSparse(M,N,Nzmax == -1 ? 0 : Nzmax,
			     si == NULL ? mxREAL : mxCOMPLEX);
    return;
  }

  if (Nii == Njj && Mjj == Mii) {
    /* cases when ii and jj have the same shape but (sr,si) have one
       of 4 different shapes */
    const int smod = Nss != Njj ? Mii : len;
    const int sdiv = Mss != Mii ? Mii : 1;
    if (sort)
      plhs[0] = sparse(ii,jj,sr,si,smod,sdiv,len,
		       M,N,Nzmax);
    else   
      plhs[0] = sparse_nosort(ii,jj,sr,si,smod,sdiv,len,
			      M,N,Nzmax);
  }
  else {
    // fully general case
    const int imod = Nii != Njj ? Mii : len;
    const int jdiv = Mjj != Mii ? Mii : 1;
    const int smod = Nss != Njj ? Mii : len;
    const int sdiv = Mss != Mii ? Mii : 1;
    if (sort)
      plhs[0] = gsparse(ii,imod,jj,jdiv,sr,si,smod,sdiv,len,
			M,N,Nzmax);
    else   
      plhs[0] = gsparse_nosort(ii,imod,jj,jdiv,sr,si,smod,sdiv,len,
			       M,N,Nzmax);
  }

  // deallocate
  if (!nocopyii) mxFree(ii);
  if (!nocopyjj) mxFree(jj);

  // squeeze out zero elements
  squeeze(plhs[0]);
}
/*------------------------------------------------------------------------*/
bool mx_IsInt(const mxArray *array_ptr)
/* Returns logical 1 (true) if array_ptr is a numeric array containing
   integers (int8, int16, int32 or int64 depending on the platform),
   and logical 0 (false) otherwise. This is useful since the test is
   not provided in MEX.

   In the name of the function, an extra underscore is used in order
   to avoid confusing it with true MEX-functions. */
{
  const int id = mxGetClassID(array_ptr);
  const size_t siz = mxGetElementSize(array_ptr);
  
  /* check that the class is an integer and that its size matches the
     size of an int */
  return (id == mxINT8_CLASS || id == mxINT16_CLASS ||
	  id == mxINT32_CLASS || id == mxINT64_CLASS) && siz == sizeof(int);
}
/*------------------------------------------------------------------------*/
bool getix(int **ix,int M,int N,int *max,bool nocopy,const mxArray *IX)
/* Gets indices ix from mxArray IX. The dimensions are M-by-N, max is
   set to the maximum index and nocopy defines the type of IX (double
   or int). */
{
  bool ok = true;
  int mx = *max;

#ifndef _OPENMP
  if (nocopy) {
    // no copy
    const int *iix = (*ix = (int *)mxGetData(IX));
    for (int i = 0; i < M*N; i++) {
      if (iix[i] < 1)
	return false;
      if (iix[i] > mx) mx = iix[i];
    }
  }
  else {
    // typecast copy
    const double *ival = mxGetPr(IX);
    int *iix = (*ix = mxMalloc(M*N*sizeof(int)));
    for (int i = 0; i < M*N; i++) {
      if (ival[i] < 1.0 || ival[i] != ceil(ival[i]))
	return false;
      if ((iix[i] = ival[i]) > mx) mx = ival[i];
    }
  }
#else // _OPENMP
  if (nocopy) {
    // no copy
    const int *iix = (*ix = (int *)mxGetData(IX));
#pragma omp parallel shared (mx)
{
    int mymx = mx; // local version of mx
#pragma omp for
    for (int i = 0; i < M*N; i++) {
      if (iix[i] > mymx)
	mymx = iix[i];
      else if (iix[i] < 1)
	ok = false; // no harm in continuing
    }

    if (mx < mymx)
#pragma omp critical
      // ensure nothing changed, then make the swap:
      if (mx < mymx) mx = mymx;
} // end omp parallel
  }
  else {
    // typecast copy
    const double *ival = mxGetPr(IX);
    int *iix;
#pragma omp critical
     // not thread-safe:
    iix = (*ix = mxMalloc(M*N*sizeof(int)));
#pragma omp parallel shared (mx)
{
    int mymx = mx; // local version of mx
#pragma omp for
    for (int i = 0; i < M*N; i++) {
      if (ival[i] < 1.0 || ival[i] != ceil(ival[i]))
	ok = false; // no harm in continuing
      else if ((iix[i] = ival[i]) > mymx)
	mymx = ival[i];
    }

    if (mx < mymx)
#pragma omp critical
      // ensure nothing changed, then make the swap:
      if (mx < mymx) mx = mymx;
} // end omp parallel
  }
#endif // _OPENMP
  *max = mx;
  return ok;
}
/*------------------------------------------------------------------------*/
mxArray *sparse2sparse(const mxArray *S)
/* Returns a deep copy T of a sparse matrix S. The allocation is
   exact. */
{
  const mwSize *jcS = mxGetJc(S);
  const mwIndex *irS = mxGetIr(S);
  const double *prS = mxGetPr(S),*piS = mxGetPi(S);
  const int N = mxGetN(S),M = mxGetM(S),Nnz = jcS[N];
  const bool real = piS == NULL;
  mxArray *T;

  /* straightforward */
  T = mxCreateSparse(M,N,Nnz,real ? mxREAL : mxCOMPLEX);
  memcpy(mxGetJc(T),jcS,(N+1)*sizeof(jcS[0]));
  memcpy(mxGetIr(T),irS,Nnz*sizeof(irS[0]));
  memcpy(mxGetPr(T),prS,Nnz*sizeof(prS[0]));
  if (!real) memcpy(mxGetPi(T),piS,Nnz*sizeof(piS[0]));

  return T;
}
/*------------------------------------------------------------------------*/
mxArray *full2sparse(const mxArray *A)
/* Constructs a sparse matrix S from a full matrix A. */
{
  const double *prA = mxGetPr(A),*piA = mxGetPi(A);
  const int N = mxGetN(A),M = mxGetM(A);
  const bool real = piA == NULL;
  mxArray *S;

  if (real) {
    mwSize *jcS = mxCalloc(N+1,sizeof(mwSize));
    mwIndex *irS;
    double *prS;

    /* determine the column pointer */
    for (int c = 1,k = 0; c <= N; c++,k += M)
      for (int i = 0; i < M; i++)
	if (prA[k+i] != 0.0) jcS[c]++;
    for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

    /* allocate */
    S = mxCreateSparse(0,0,jcS[N],mxREAL);
    mxSetM(S,M);
    mxSetN(S,N);
    irS = mxGetIr(S);
    prS = mxGetPr(S);

    /* set the column pointer */
    mxFree(mxGetJc(S));
    mxSetJc(S,jcS);

    /* copy data */
    for (int c = 1,k = 0,dest = 0; c <= N; c++,k += M)
      for (int i = 0; i < M; i++)
	if (prA[k+i] != 0.0) {
	  irS[dest] = i;
	  prS[dest++] = prA[k+i];
	}
  }
  else {
    mwSize *jcS = mxCalloc(N+1,sizeof(mwSize));
    mwIndex *irS;
    double *prS,*piS;

    for (int c = 1,k = 0; c <= N; c++,k += M)
      for (int i = 0; i < M; i++)
	if (prA[k+i] != 0.0 && piA[k+i] != 0.0) jcS[c]++;
    for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

    S = mxCreateSparse(0,0,jcS[N],mxCOMPLEX);
    mxSetM(S,M);
    mxSetN(S,N);
    irS = mxGetIr(S);
    prS = mxGetPr(S);
    piS = mxGetPi(S);

    mxFree(mxGetJc(S));
    mxSetJc(S,jcS);

    for (int c = 1,k = 0,dest = 0; c <= N; c++,k += M)
      for (int i = 0; i < M; i++)
	if (prA[k+i] != 0.0 && piA[k+i] != 0.0) {
	  irS[dest] = i;
	  prS[dest] = prA[k+i];
	  piS[dest++] = piA[k+i];
	}
  }

  return S;
}
/*------------------------------------------------------------------------*/
void squeeze(mxArray *S)
/* Removes any zero elements explicitly stored in the sparse matrix
   S. No reallocation is performed.

   There is a quite complicated and potentially faster algorithm based
   on memmove() that performs the same operation. However, benchmark
   tests indicate that the following simple code optimizes better. */
{
  const int N = mxGetN(S);
  mwSize *jcS = mxGetJc(S);
  mwIndex *irS = mxGetIr(S);
  double *prS = mxGetPr(S),*piS = mxGetPi(S);

  if (piS == NULL) {
    int c,i;

    /* find the first zero, if any */
    for (i = c = 0; c < N; c++)
      for ( ; i < jcS[c+1]; i++)
	if (prS[i] == 0.0)
	  goto rfound0; /* a 'double break' */
    return;

  rfound0:
    /* copy in a conservative fashion */
    for (int dest = i++; c < N; c++) {
      for ( ; i < jcS[c+1]; i++)
	if (prS[i] != 0.0) {
	  irS[dest] = irS[i];
	  prS[dest++] = prS[i];
	}
      jcS[c+1] = dest;
    }
  }
  else {
    int c,i;
    for (i = c = 0; c < N; c++)
      for ( ; i < jcS[c+1]; i++)
	if (prS[i] == 0.0 && piS[i] == 0.0)
	  goto zfound0;
    return;

  zfound0:
    for (int dest = i++; c < N; c++) {
      for ( ; i < jcS[c+1]; i++)
	if (prS[i] != 0.0 || piS[i] != 0.0) {
	  irS[dest] = irS[i];
	  prS[dest] = prS[i];
	  piS[dest++] = piS[i];
	}
      jcS[c+1] = dest;
    }
  }
}
/*------------------------------------------------------------------------*/
void sparse_insert(mwIndex *irS,double *prS,double *piS,
		   const int *irank,const int *rank,const mwIndex *jrS,
		   const int *ii,const double *sr,const double *si,
		   int smod,int sdiv,int len,int M)
/* Inserts elements into sparse matrix. Input is the sparse matrix
   itself (irS,prS,piS), an index-table irank, the rowindices ii and
   the values of the elements (sr,si). Four different formats of the
   values are allowed as indicated by the parameters smod and sdiv.

   Currently, input rank, jrS, and M are only used #ifdef _OPENMP and
   for the full case 3 below. */
{
  const bool real = si == NULL;

  switch (2*(smod == len)+(sdiv == 1)) {
  case 3 : /* full case */
#ifndef _OPENMP
    if (real)
      for (int i = 0; i < len; i++) {
      	irS[irank[i]] = ii[i]-1;
      	prS[irank[i]] += sr[i];
      }
    else
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i]-1;
	prS[irank[i]] += sr[i];
	piS[irank[i]] += si[i];
      }
#else // _OPENMP
    if (real) {

    if (rank != NULL) {
      /* needed since not all cases respond to _OPENMP and uses a
	 different syntax */
#pragma omp parallel
{
      const int nThreads = omp_get_num_threads();
      const int myId = omp_get_thread_num();
      const int rstart = 1+M*myId/nThreads;
      const int rend = M*(myId+1)/nThreads;
      int istart;
      if (rstart == 1)
	istart = 0;
      else
	istart = jrS[rstart-1];

      if (rend >= 1) {
        for (int i = istart; i < jrS[rend]; i++) irS[irank[i]] = ii[rank[i]]-1;
        for (int i = istart; i < jrS[rend]; i++) prS[irank[i]] += sr[rank[i]];
      }
} // end parallel
    }
    else {
#pragma omp single nowait
      for (int i = 0; i < len; i++) irS[irank[i]] = ii[i]-1;
#pragma omp single nowait
      for (int i = 0; i < len; i++) prS[irank[i]] += sr[i];
    }
    }
    else {
#pragma omp parallel
{
  if (rank != NULL) {
    const int nThreads = omp_get_num_threads();
    const int myId = omp_get_thread_num();
    const int rstart = 1+M*myId/nThreads;
    const int rend = M*(myId+1)/nThreads;
    int istart;
    if (rstart == 1)
      istart = 0;
    else
      istart = jrS[rstart-1];

    if (rend >= 1)
      for (int i = istart; i < jrS[rend]; i++) {
	irS[irank[i]] = ii[rank[i]]-1;
	prS[irank[i]] += sr[rank[i]];
	piS[irank[i]] += si[rank[i]];
      }
  }
  else {
#pragma omp single nowait
      for (int i = 0; i < len; i++) irS[irank[i]] = ii[i]-1;
#pragma omp single nowait
      for (int i = 0; i < len; i++) prS[irank[i]] += sr[i];
#pragma omp single nowait
      for (int i = 0; i < len; i++) piS[irank[i]] += si[i];
  }
} // end parallel
    }
#endif // _OPENMP
    break;

  case 2 : /* horizontal case */
    if (real)
      for (int j = 0; j < len; j += sdiv) {
	const double ssr = sr[j/sdiv];
	for (int i = j; i < j+sdiv; i++) {
	  irS[irank[i]] = ii[i]-1;
	  prS[irank[i]] += ssr;
	}
      }
    else
      for (int j = 0; j < len; j += sdiv) {
	const double ssr = sr[j/sdiv];
	const double ssi = si[j/sdiv];
	for (int i = j; i < j+sdiv; i++) {
	  irS[irank[i]] = ii[i]-1;
	  prS[irank[i]] += ssr;
	  piS[irank[i]] += ssi;
	}
      }
    break;

  case 1 : /* vertical case */
    if (real)
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i]-1;
	prS[irank[i]] += sr[i%smod];
      }
    else
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i]-1;
	prS[irank[i]] += sr[i%smod];
	piS[irank[i]] += si[i%smod];
      }
    break;

  case 0 : /* scalar case */
#ifndef _OPENMP
    if (real) {
      const double ssr = sr[0];
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i]-1;
	prS[irank[i]] += ssr;
      }
    }
    else {
      const double ssr = sr[0];
      const double ssi = si[0];
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i]-1;
	prS[irank[i]] += ssr;
	piS[irank[i]] += ssi;
      }
    }
#else // _OPENMP
    if (real) {
#pragma omp parallel
{
      const double ssr = sr[0];
#pragma omp single nowait
      for (int i = 0; i < len; i++) irS[irank[i]] = ii[i]-1;
#pragma omp single nowait
      for (int i = 0; i < len; i++) prS[irank[i]] += ssr;
} // end omp parallel
    }
    else {
#pragma omp parallel
{
      const double ssr = sr[0];
      const double ssi = si[0];
#pragma omp single nowait
      for (int i = 0; i < len; i++) irS[irank[i]] = ii[i]-1;
#pragma omp single nowait
      for (int i = 0; i < len; i++) prS[irank[i]] += ssr;
#pragma omp single nowait
      for (int i = 0; i < len; i++) piS[irank[i]] += ssi;
} // end omp parallel
    }
#endif // _OPENMP
    break;
  }
}
/*------------------------------------------------------------------------*/
void sparse_inserti(mwIndex *irS,double *prS,double *piS,
		    const int *irank,
		    const int *ii,int imod,
		    const double *sr,const double *si,
		    int smod,int sdiv,int len)
/* Same as sparse_insert() above except that ii is assumed to be
   vertical. */
{
  const bool real = si == NULL;

  switch (2*(smod == len)+(sdiv == 1)) {
  case 3 : /* full case */
    if (real)
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i%imod]-1;
	prS[irank[i]] += sr[i];
      }
    else
      for (int i = 0; i < len; i++) {
	irS[irank[i]] = ii[i%imod]-1;
	prS[irank[i]] += sr[i];
	piS[irank[i]] += si[i];
      }
    break;
  case 2 : /* horizontal case */
    if (real)
      for (int j = 0; j < len; j += sdiv) {
	const double ssr = sr[j/sdiv];
	for (int i = 0; i < imod; i++) {
	  irS[irank[j+i]] = ii[i]-1;
	  prS[irank[j+i]] += ssr;
	}
      }
    else
      for (int j = 0; j < len; j += sdiv) {
	const double ssr = sr[j/sdiv];
	const double ssi = si[j/sdiv];
	for (int i = 0; i < imod; i++) {
	  irS[irank[j+i]] = ii[i]-1;
	  prS[irank[j+i]] += ssr;
	  piS[irank[j+i]] += ssi;
	}
      }
    break;
  case 1 : /* vertical case */
    if (real)
      for (int i = 0; i < imod; i++) {
	const double ssr = sr[i];
	const int iii = ii[i]-1;
	for (int j = i; j < len; j += imod) {
	  irS[irank[j]] = iii;
	  prS[irank[j]] += ssr;
	}
      }
    else
      for (int i = 0; i < imod; i++) {
	const double ssr = sr[i];
	const double ssi = si[i];
	const int iii = ii[i]-1;
	for (int j = i; j < len; j += imod) {
	  irS[irank[j]] = iii;
	  prS[irank[j]] += ssr;
	  piS[irank[j]] += ssi;
	}
      }
    break;
  case 0 : /* scalar case */
    if (real) {
      const double ssr = sr[0];
      for (int j = 0; j < len; j += imod)
	for (int i = 0; i < imod; i++) {
	  irS[irank[j+i]] = ii[i]-1;
	  prS[irank[j+i]] += ssr;
	}
    }
    else {
      const double ssr = sr[0];
      const double ssi = si[0];
      for (int j = 0; j < len; j += imod)
	for (int i = 0; i < imod; i++) {
	  irS[irank[j+i]] = ii[i]-1;
	  prS[irank[j+i]] += ssr;
	  piS[irank[j+i]] += ssi;
	}
    }
    break;
  }
}
/*------------------------------------------------------------------------*/
#ifndef _OPENMP
mxArray *sparse(const int *ii,const int *jj,
		const double *sr,const double *si,
		int smod,int sdiv,
		int len,int M,int N,int Nzmax)
/* Constructs a sparse matrix in Compressed Column Storage (CCS) from
   triplet format [ii,jj,sr(si)]. An ordinary Matlab sparse matrix is
   thus constructed.

   Input (smod,sdiv) determine the shape of the value array sr(si) and
   len is the length of the index arrays which both must have the same
   shape. [M N Nzmax] determine the dimensions of the resulting
   matrix. If Nzmax = -1, then the allocation is exact (i.e. nnz(S) =
   nzmax(S)). Otherwise, Nzmax must be greater than or equal to the
   number of nonzeros needed to be stored.

   The memory demand of the algorithm is (at peak and for sufficiently
   large output) <input>+int[len]+<output>. */
{
  // output
  mxArray *S;
  mwSize *jcS; // column pointer for sparse matrix S
  int *irank;  // inverse rank array of length len

  mwSize *jrS; // accumulated "pessimistic" row counter
  int *rank;   // rank-array for rows
  int *hcol;   // cache memory for columns

  // Part 1: count and accumulate indices to rows
  StartTime;
  jrS = mxCalloc(M+1,sizeof(jrS[0]));
  for (int i = 0; i < len; i++) jrS[ii[i]]++;
  for (int r = 2; r <= M; r++) jrS[r] += jrS[r-1];
  StopTime;
  GetTime(time_vec[1]);
#ifdef PRINT_INTERMEDIATE
  mexPrintf("jrS = [");
  for (int r = 0; r <= M; r++) mexPrintf("%d,",jrS[r]);
  mexPrintf("]\n\n");
#endif

  // Part 2: build rank with the active use of jrS
  StartTime;
  rank = mxMalloc(len*sizeof(rank[0]));
  jrS--; /* (unit-offset in ii) */
  for (int i = 0; i < len; i++) rank[jrS[ii[i]]++] = i;
  // rank now allows for row-wise traversal
  StopTime;
  GetTime(time_vec[2]);
#ifdef PRINT_INTERMEDIATE
  mexPrintf("rank = [");
  for (int i = 0; i < len; i++) mexPrintf("%d,",rank[i]);
  mexPrintf("]\n");
  mexPrintf("jrS = [*,");
  for (int r = 1; r <= M+1; r++) mexPrintf("%d,",jrS[r]);
  mexPrintf("]\n\n");
#endif

  /* Part 3: loop over input and make each column unique with respect
     to rowindices, building both an index vector irank and the final
     column pointer at the same time */
  StartTime;
  jcS = mxCalloc(N+1,sizeof(jcS[0]));
  hcol = mxCalloc(N,sizeof(hcol[0]));
  hcol--; /* (unit-offset in jj) */
  irank = mxMalloc(len*sizeof(irank[0]));
  for (int row = 1,i = 0; row <= M; row++)
    for ( ; i < jrS[row]; i++) {
      const int ixijs = rank[i]; // index into input data triplet (ii,jj,sr)
      const int col = jj[ixijs]; // column index

      // new element?
      if (hcol[col] < row) {
	hcol[col] = row; // remembered by the row index
	jcS[col]++;      // count it
      }

      // irank keeps track of where it should go
      irank[ixijs] = jcS[col]-1;
    }
  mxFree(++hcol);
  mxFree(rank);
  mxFree(++jrS);
  StopTime;
  GetTime(time_vec[3]);
#ifdef PRINT_INTERMEDIATE
  mexPrintf("irank = [");
  for (int i = 0; i < len; i++) mexPrintf("%d,",irank[i]);
  mexPrintf("]\n");
  mexPrintf("jcS = [");
  for (int c = 0; c <= N; c++) mexPrintf("%d,",jcS[c]);
  mexPrintf("]\n\n");
#endif
  
  // Part 4: accumulate pointer to columns
  StartTime;
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  // irank must account for the previous accumulation
  jcS--; /* (again, unit-offset in jj) */
  for (int i = 0; i < len; i++) irank[i] += jcS[jj[i]];
  jcS++;
  StopTime;
  GetTime(time_vec[4]);
#ifdef PRINT_INTERMEDIATE
  mexPrintf("irank = [");
  for (int i = 0; i < len; i++) mexPrintf("%d,",irank[i]);
  mexPrintf("]\n");
  mexPrintf("jcS = [");
  for (int c = 0; c <= N; c++) mexPrintf("%d,",jcS[c]);
  mexPrintf("]\n\n");
#endif

  // allocate output
  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  // set the column pointer
  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  // insert the data
  StartTime;
  sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		irank,0,0,ii,sr,si,smod,sdiv,len,M);
  StopTime;
  GetTime(time_vec[5]);

  mxFree(irank);
  return S;
}
/*------------------------------------------------------------------------*/
#else // _OPENMP
/*------------------------------------------------------------------------*/
mxArray *sparse(const int *ii,const int *jj,
		const double *sr,const double *si,
		int smod,int sdiv,
		int len,int M,int N,int Nzmax)
/* This is the OpenMP-version of the sparse function above. */
{
  // output
  mxArray *S;
  mwSize **jcS; // column pointer, one per thread
  mwSize *jcS_; // final column pointer
  int *irank;   // inverse rank array of length len
  int *irankP;  // permuted version of irank

  mwSize **jrS; // accumulated "pessimistic" row counter
  int *rank;    // rank-array for rows

  // Part 1: count and accumulate indices to rows
  StartTime;
  const int nThreads = omp_get_max_threads();
  jrS = mxMalloc((nThreads+1)*sizeof(jrS[0]));
  for (int k = 0; k <= nThreads; k++) {
    jrS[k] = mxCalloc(M+1,sizeof(jrS[k][0]));
    jrS[k]--; /* (unit-offset in ii) */
  }

#pragma omp parallel
{
  const int myId = omp_get_thread_num();
  const int istart = len*myId/nThreads;
  const int iend = len*(myId+1)/nThreads;
  for (int i = istart; i < iend; i++)
    jrS[myId+1][ii[i]]++;

#pragma omp barrier

  // accumulate jrS over the threads
#pragma omp for
  for (int r = 1; r <= M; r++)
    for (int k = 1; k < nThreads; k++)
      jrS[k+1][r] += jrS[k][r];

  // serial accumulation in jrS[0]
#pragma omp single
  for (int r = 1; r <= M; r++)
    jrS[0][r+1] += jrS[0][r]+jrS[nThreads][r];

  // determine a private jrS for each thread
#pragma omp for
  for (int r = 1; r <= M; r++)
    for (int k = 1; k < nThreads; k++)
      jrS[k][r] += jrS[0][r];
} // end parallel
  StopTime;
  GetTime(time_vec[1]);

  // Part 2: build rank with the active use of jrS
  StartTime;
  rank = mxMalloc(len*sizeof(rank[0]));

#pragma omp parallel
{
  const int myId = omp_get_thread_num();
  const int istart = len*myId/nThreads;
  const int iend = len*(myId+1)/nThreads;
  for (int i = istart; i < iend; i++)
    rank[jrS[myId][ii[i]]++] = i;
  // rank now allows for row-wise traversal
} // end parallel
  StopTime;
  GetTime(time_vec[2]);

  /* Part 3: loop over input and make each column unique with respect
     to rowindices, building both a permuted index vector irankP and
     the final column pointer at the same time */
  StartTime;
  jcS = mxMalloc((nThreads+1)*sizeof(jcS[0]));
  for (int k = 0; k <= nThreads; k++)
    jcS[k] = mxCalloc(N+1,sizeof(jcS[k][0]));
  irankP = mxMalloc(len*sizeof(irankP[0]));
  if (2*(smod == len)+(sdiv == 1) != 3)
    // *** case not fully implemented
    irank = mxMalloc(len*sizeof(irank[0]));

#pragma omp parallel
{
  int *hcol; // cache memory for columns
#pragma omp critical
  hcol = mxCalloc(N,sizeof(hcol[0]));
  hcol--; /* (unit-offset in jj) */

  const int myId = omp_get_thread_num();
  const int rstart = 1+M*myId/nThreads;
  const int rend = M*(myId+1)/nThreads;
  int istart;
  if (rstart == 1)
    istart = 0;
  else
    istart = jrS[nThreads-1][rstart-1];

  // loop over segment of row indices
  for (int row = rstart,i = istart; row <= rend; row++)
    // loop over single row
    for ( ; i < jrS[nThreads-1][row]; i++) {
      const int ixijs = rank[i]; // index into input data triplet (ii,jj,sr)
      const int col = jj[ixijs]; // column index

      // new element?
      if (hcol[col] < row) {
      	hcol[col] = row;  // store row index
      	jcS[myId+1][col]++; // count it
      }

      // irankP keeps track of where it should go
      irankP[i] = jcS[myId+1][col]-1;
    }
#pragma omp critical
  mxFree(++hcol);

#pragma omp barrier

  // accumulate jcS over the threads
#pragma omp for
  for (int c = 1; c <= N; c++)
    for (int k = 1; k < nThreads; k++)
      jcS[k+1][c] += jcS[k][c];

  // serial accumulation in jcS[0]
#pragma omp single
{
  for (int c = 1; c <= N; c++)
    jcS[0][c] += jcS[0][c-1]+jcS[nThreads][c];
  jcS[0]--; /* (unit-offset in jj) */
}

  // determine a private jcS for each thread
#pragma omp for
  for (int c = 1; c <= N; c++)
    for (int k = 1; k < nThreads; k++)
      jcS[k][c] += jcS[0][c];

  // irankP must now account to these changes to jcS
  if (rend >= 1)
    for (int i = istart; i < jrS[nThreads-1][rend]; i++)
      irankP[i] += jcS[myId][jj[rank[i]]];

  if (2*(smod == len)+(sdiv == 1) != 3) {
    // *** case not fully implemented
    if (rend >= 1)
      for (int i = istart; i < jrS[nThreads-1][rend]; i++)
	irank[rank[i]] = irankP[i];
  }
} // end parallel
  StopTime;
  GetTime(time_vec[3]);

  // Part 4: final accumulation of jcS
  StartTime;
  /* (nil in this version) */
  StopTime;
  GetTime(time_vec[4]);

  // allocate output
  jcS[0]++;
  if (Nzmax == -1)
    Nzmax = jcS[0][N];
  else if (Nzmax < jcS[0][N])
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  // set the column pointer
  mxFree(mxGetJc(S));
  mxSetJc(S,jcS[0]);
  // free the remaining pointers
  for (int k = 1; k <= nThreads; k++)
    mxFree(jcS[k]);
  mxFree(jcS);

  // insert the data
  StartTime;
  if (2*(smod == len)+(sdiv == 1) == 3)
    sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		  irankP,rank,jrS[nThreads-1],ii,sr,si,smod,sdiv,len,M);
  else
    sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		  irank,0,0,ii,sr,si,smod,sdiv,len,M);
  StopTime;
  GetTime(time_vec[5]);

  // *** case not fully implemented
  if (2*(smod == len)+(sdiv == 1) != 3)
    mxFree(irank);
  mxFree(irankP);

  return S;
}
#endif // _OPENMP
/*------------------------------------------------------------------------*/

/* table is only allowed to become about 70% full */
const double HASH_LOAD = 1.0/0.70;

size_t hashsize(size_t len)
/* Size of hash table. This is simply
   exp2(ceil(log2(HASH_LOAD*len))). */
{
  size_t ans = 1;

  len = ceil(HASH_LOAD*len);
  do ans <<= 1; while (len >>= 1);

  return ans;
}
/*------------------------------------------------------------------------*/
uint32_t hash(int ii,int jj)
/* Computes a hashkey from an integer pair (ii,jj).

   Adapted from a code by Paul Hsieh. See
   http://www.azillionmonkeys.com/qed/hash.html for further details.
*/
{
  uint32_t hash = jj,temp;

  /* main loop executed twice */
  hash += ii&0x0000FFFF;
  temp = ((ii&0xFFFF0000) >> 5)^hash;
  hash = (hash << 16)^temp;
  hash += hash >> 11;

  hash += jj&0x0000FFFF;
  temp = ((jj&0xFFFF0000) >> 5)^hash;
  hash = (hash << 16)^temp;
  hash += hash >> 11;

  /* "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 2;
  hash += hash >> 15;
  hash ^= hash << 10;

  return hash;
}
/*------------------------------------------------------------------------*/
mxArray *sparse_nosort(const int *ii,const int *jj,
		       const double *sr,const double *si,
		       int smod,int sdiv,
		       int len,int M,int N,int Nzmax)
/* This function does the same job as sparse() above except that the
   result is *not* sorted w.r.t. rowindices. The sparse matrix thus
   constructed is therefore a sparse matrix which is valid in many,
   but not all situations in Matlab.

   The allocation is still exact (i.e. nzmax(S) = nnz(S)) and the
   memory demand is (at peak and for sufficiently large output)
   <input>+int[N+2]+int[len]+<output>, where len is the size of the
   input index set. */
{

#define SPS1 /* new version with hashing */
#undef SPS2  /* updated version without hashing */

#ifdef SPS1

  /* output */
  mxArray *S;
  mwSize *jcS;

  /* rank- and hash-table */
  int *irank,*hash_tb;
  size_t hash_sz = hashsize(len);

  jcS = mxCalloc(N+1,sizeof(jcS[0]));
  irank = mxMalloc(len*sizeof(irank[0]));
  hash_tb = mxCalloc(hash_sz,sizeof(hash_tb[0]));

  for (int i = 0; i < len; i++) {
    const int row = ii[i],col = jj[i];
    uint32_t hashkey = hash(row,col)&(hash_sz-1);

    /* collision? */
    if (hash_tb[hashkey] != 0 && 
	!(row == ii[hash_tb[hashkey]-1] && 
	  col == jj[hash_tb[hashkey]-1])) {
      /* a fairly independent and odd (incremental) hashkey */
      const uint32_t hashkey2nd = hash(col,row)|1;

      /* rehash until empty slot found */
      do
	hashkey = (hashkey+hashkey2nd)&(hash_sz-1);
      while (hash_tb[hashkey] != 0 && 
	     !(row == ii[hash_tb[hashkey]-1] && 
	       col == jj[hash_tb[hashkey]-1]));
    }

    /* fill in table */
    if (hash_tb[hashkey] == 0) {
      hash_tb[hashkey] = i+1;
      irank[i] = jcS[col]++;
    }
    else
      irank[i] = irank[hash_tb[hashkey]-1];
  }
  mxFree(hash_tb);
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  /* account for the accumulation of indices */
  jcS--;
  for (int i = 0; i < len; i++) irank[i] += jcS[jj[i]];
  jcS++;

  /* allocate output */
  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  /* set the column pointer */
  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  /* insert the data */
  sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		irank,0,0,ii,sr,si,smod,sdiv,len,M);
  mxFree(irank);

  return S;

#elif defined(SPS2)

  /* output */
  mxArray *S;
  int *jcS;

  /* rank-tables and help-pointer for rowindices */
  int *rank,*irank,*hrow;

  /* determine a "pessimistic" column pointer */
  jcS = mxCalloc(N+2,sizeof(jcS[0]));
  jcS++; /* allows for a fast restore */
  for (int i = 0; i < len; i++) jcS[jj[i]]++;
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  /* produce the rank-table so that the data may be traversed in order
     with respect to columns */
  rank = mxMalloc(len*sizeof(rank[0]));
  jcS--; /* restore jcS and account for unit offset in jj */
  for (int i = 0; i < len; i++) rank[jcS[jj[i]]++] = i;

  /* allocate index- and help-pointers */
  irank = mxMalloc(len*sizeof(irank[0]));
  hrow = mxCalloc(M,sizeof(hrow[0]));
  hrow--; /* unit offset in ii */

  /* loop over input and make each column unique with respect to
     rowindices, building both an index vector irank and the final
     column pointer at the same time */
  for (int col = 1,i = 0; col <= N; col++) {
    const int begin = jcS[col-1],end = jcS[col];
    jcS[col] = begin;
    for ( ; i < end; i++) {
      const int ixijs = rank[i];
      const int row = ii[ixijs];

      /* new element; mark and count it */
      if (hrow[row] <= begin) hrow[row] = ++jcS[col];

      /* remember where it should go */
      irank[ixijs] = hrow[row]-1;
    }
  }
  mxFree(++hrow);
  mxFree(rank);

  /* allocate output */
  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  /* set the column pointer */
  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  /* insert the data */
  sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		irank,0,0,ii,sr,si,smod,sdiv,len,M);
  mxFree(irank);

  return S;

#endif
}
/*------------------------------------------------------------------------*/
/*- General versions -----------------------------------------------------*/
/*------------------------------------------------------------------------*/
mxArray *gsparse(const int *ii,int imod,
		 const int *jj,int jdiv,
		 const double *sr,const double *si,
		 int smod,int sdiv,
		 int len,int M,int N,int Nzmax)
/* Derived from sparse(). */
{

#define GSP1 /* new version */
#undef GSP2  /* (very) old version */

#ifdef GSP1

  mxArray *S;
  mwSize *jcS;
  mwIndex *irS;

  int *rank,*irank,*hcol;

  irS = mxCalloc(M+1,sizeof(irS[0]));
  /* for (int i = 0; i < len; i++) irS[ii[i%imod]]++; */
  if (imod == len)
    for (int i = 0; i < len; i++) irS[ii[i]]++;
  else
    for (int i = 0; i < imod; i++) irS[ii[i]] += len/imod;
  for (int r = 2; r <= M; r++) irS[r] += irS[r-1];

  rank = mxMalloc(len*sizeof(rank[0]));
  irS--;
  for (int i = 0; i < len; i++) rank[irS[ii[i%imod]]++] = i;

  jcS = mxCalloc(N+1,sizeof(jcS[0]));
  hcol = mxCalloc(N,sizeof(hcol[0]));
  hcol--;
  irank = mxMalloc(len*sizeof(irank[0]));
  for (int row = 1,i = 0; row <= M; row++)
    for ( ; i < irS[row]; i++) {
      const int ixijs = rank[i];
      const int col = jj[ixijs/jdiv];

      if (hcol[col] < row) {
	hcol[col] = row;
	jcS[col]++;
      }

      irank[ixijs] = jcS[col]-1;
    }
  mxFree(++irS);
  mxFree(rank);
  mxFree(++hcol);
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  jcS--;
  /* for (int i = 0; i < len; i++) irank[i] += jcS[jj[i/jdiv]]; */
  if (jdiv == 1)
    for (int i = 0; i < len; i++) irank[i] += jcS[jj[i]];
  else
    for (int c = 0, i = 0; c < len/jdiv; c++)
      for ( ; i < (c+1)*jdiv; i++) irank[i] += jcS[jj[c]];
  jcS++;

  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  if (imod == len)
    sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		  irank,0,0,ii,sr,si,smod,sdiv,len,M);
  else
    sparse_inserti(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		   irank,ii,imod,sr,si,smod,sdiv,len);
  mxFree(irank);

  return S;

#elif defined(GSP2)

  mxArray *S;
  int *jcS,*irS;

  int *rank,*irank;

  jcS = mxCalloc(N+2,sizeof(jsC[0]));
  jcS++;
  irS = mxCalloc(M+1,sizeof(irS[0]));

  for (int i = 0; i < len; i++) {
    jcS[jj[i/jdiv]]++;
    irS[ii[i%imod]]++;
  }
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];
  for (int r = 2; r <= M; r++) irS[r] += irS[r-1];

  irank = mxMalloc(len*sizeof(irank[0]));
  irS--;
  for (int i = 0; i < len; i++)
    irank[irS[ii[i%imod]]++] = i;
  mxFree(++irS);

  rank = mxMalloc(len*sizeof(rank[0]));
  jcS--;
  for (int i = 0; i < len; i++)
    rank[jcS[jj[irank[i]/jdiv]]++] = irank[i];

  for (int i = 0,col = jj[rank[0]/jdiv],row; ; ) {
    const int end = jcS[col];

    irank[rank[i]] = jcS[col-1];
    jcS[col] = jcS[col-1]+1;
    row = ii[rank[i]%imod];

    for (i++; i < end; i++) {
      const int ixijs = rank[i];

      if (row < ii[ixijs%imod]) {
	row = ii[ixijs%imod];
	jcS[col]++;
      }

      irank[ixijs] = jcS[col]-1;
    }

    if (i < len)
      for (col++; col < jj[rank[i]/jdiv]; col++)
	jcS[col] = jcS[col-1];
    else {
      for (col++; col <= N; col++) jcS[col] = jcS[col-1];
      break;
    }
  }
  mxFree(rank);

  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  if (imod == len)
    sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		  irank,0,0,ii,sr,si,smod,sdiv,len,M);
  else
    sparse_inserti(mxGetIr(S),mxGetPr(S),mxGetPi(S),
		   irank,ii,imod,sr,si,smod,sdiv,len);
  mxFree(irank);

  return S;

#endif

}
/*------------------------------------------------------------------------*/
mxArray *gsparse_nosort(const int *ii,int imod,
			const int *jj,int jdiv,
			const double *sr,const double *si,
			int smod,int sdiv,
			int len,int M,int N,int Nzmax)
/* Derived from sparse_nosort(). */
{

#define GSPS1 /* new version using hashing */
#undef GSPS2  /* (very) old version without hashing */

#ifdef GSPS1

  /* output */
  mxArray *S;
  mwSize *jcS;

  /* rank- and hash-table */
  int *irank,*hash_tb;
  size_t hash_sz = hashsize(len);

  jcS = mxCalloc(N+1,sizeof(jcS[0]));
  irank = mxMalloc(len*sizeof(irank[0]));
  hash_tb = mxCalloc(hash_sz,sizeof(hash_tb[0]));

  for (int i = 0; i < len; i++) {
    const int row = ii[i%imod],col = jj[i/jdiv];
    uint32_t hashkey = hash(row,col)&(hash_sz-1);

    /* collision? */
    if (hash_tb[hashkey] != 0 && 
	!(row == ii[(hash_tb[hashkey]-1)%imod] && 
	  col == jj[(hash_tb[hashkey]-1)/jdiv])) {
      /* a fairly independent and odd (incremental) hashkey */
      const uint32_t hashkey2nd = hash(col,row)|1;

      /* rehash until empty slot found */
      do
	hashkey = (hashkey+hashkey2nd)&(hash_sz-1);
      while (hash_tb[hashkey] != 0 && 
	     !(row == ii[(hash_tb[hashkey]-1)%imod] && 
	       col == jj[(hash_tb[hashkey]-1)/jdiv]));
    }

    /* fill in table */
    if (hash_tb[hashkey] == 0) {
      hash_tb[hashkey] = i+1;
      irank[i] = jcS[col]++;
    }
    else
      irank[i] = irank[hash_tb[hashkey]-1];
  }
  mxFree(hash_tb);
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  /* account for the accumulation of indices */
  jcS--;
  /* for (int i = 0; i < len; i++) irank[i] += jcS[jj[i/jdiv]]; */
  if (jdiv == 1)
    for (int i = 0; i < len; i++) irank[i] += jcS[jj[i]];
  else
    for (int c = 0, i = 0; c < len/jdiv; c++)
      for ( ; i < (c+1)*jdiv; i++) irank[i] += jcS[jj[c]];
  jcS++;

  /* allocate output */
  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  /* set the column pointer */
  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  /* insert the data */
  if (imod == len)
    sparse_insert(mxGetIr(S), mxGetPr(S),mxGetPi(S),
		  irank,0,0,ii,sr,si,smod,sdiv,len,M);
  else
    sparse_inserti(mxGetIr(S), mxGetPr(S),mxGetPi(S),
		   irank,ii,imod,sr,si,smod,sdiv,len);
  mxFree(irank);

  return S;

#elif defined(GSPS2)

  mxArray *S;
  int *jcS;

  int *rank,*irank,*hrow;

  jcS = mxCalloc(N+2,sizeof(jsC[0]));
  jcS++;
  for (int i = 0; i < len; i++) jcS[jj[i/jdiv]]++;
  for (int c = 2; c <= N; c++) jcS[c] += jcS[c-1];

  rank = mxMalloc(len*sizeof(rank[0]));
  jcS--;
  for (int i = 0; i < len; i++) rank[jcS[jj[i/jdiv]]++] = i;

  irank = mxMalloc(len*sizeof(irank[0]));
  hrow = mxCalloc(M,sizeof(hrow[0]));
  hrow--;

  for (int i = 0,col = jj[rank[0]/jdiv]; ; ) {
    const int begin = jcS[col-1];
    const int end = jcS[col];
    for (jcS[col] = begin; i < end; i++) {
      const int row = ii[rank[i]%imod];

      if (hrow[row] <= begin) hrow[row] = ++jcS[col];

      irank[rank[i]] = hrow[row]-1;
    }

    if (i < len)
      for (col++; col < jj[rank[i]/jdiv]; col++) jcS[col] = jcS[col-1];
    else {
      for (col++; col <= N; col++) jcS[col] = jcS[col-1];
      break;
    }
  }
  mxFree(++hrow);
  mxFree(rank);

  if (Nzmax == -1)
    Nzmax = jcS[N];
  else if (Nzmax < jcS[N]) {
    mxFree(irank);
    mxFree(jcS);
    mexErrMsgIdAndTxt("fsparse:e16","Allocation limited by caller: "
		      "sparse matrix does not fit.");
  }
  S = mxCreateSparse(0,0,Nzmax,si == NULL ? mxREAL : mxCOMPLEX);
  mxSetM(S,M);
  mxSetN(S,N);

  mxFree(mxGetJc(S));
  mxSetJc(S,jcS);

  if (imod == len)
    sparse_insert(mxGetIr(S), mxGetPr(S),mxGetPi(S),
		  irank,0,0,ii,sr,si,smod,sdiv,len,M);
  else
    sparse_inserti(mxGetIr(S), mxGetPr(S),mxGetPi(S),
		   irank,ii,imod,sr,si,smod,sdiv,len);
  mxFree(irank);

  return S;

#endif

}
/*------------------------------------------------------------------------*/
