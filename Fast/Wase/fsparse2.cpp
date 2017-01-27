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

#include "gnuSort.h"
#include <assert.h>

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










			//#include "mex.h"
			//#include "matrix.h"
			//#include <iostream>
			//#include <omp.h>

			//#include <cstdlib>

			#include <parallel/algorithm>

			struct elementGnuSort{
			  int i;
			  int j;
			  double val;
			};

			bool compare_col_2(struct elementGnuSort ele1, struct elementGnuSort ele2){
			  return ele1.i < ele2.i;
			};

			int sumGnuSort(const int* array, const int len){
			  int out = 0;
			  for (size_t i = 0; i < len; i++) {
			    out += array[i];
			  }
			  return out;
			};

			bool operator<(const elementGnuSort& lhs, const elementGnuSort& rhs){ return lhs.j < rhs.j; };

			mxArray* gnuSort(const int *ii_in,const int *jj,
			    const double *sr,
			    int len,int M,int N,int Nzmax){

				printf("%d %d %d %d \n",len, M, N, Nzmax);
				//return;
				int c;
				printf("press enter!");
				c = getchar( );

			      struct elementGnuSort* list;
			      list = (struct elementGnuSort*) malloc(len*sizeof(struct elementGnuSort));
			      #pragma omp parallel for
			      for (int iter = 0; iter < len; iter++) {
			        printf("%d\n", ii_in[iter]);
			        list[iter].i = ii_in[iter];
			        list[iter].j = jj[iter];
			        list[iter].val = sr[iter]; //Let's ignore imaginary data for now. FIX!
			      }

			      //printf("print0");
			      printf("press enter!");
			      c = getchar( );
			      //Sort the entire list according to column index.
			      __gnu_parallel::sort(list, list+len, __gnu_parallel::multiway_mergesort_tag());

			      //printf("print1");
			      printf("press enter!");
			      c = getchar( );
			      	mwSize *jcS_d;
					jcS_d = (mwSize*) mxCalloc(N+1, sizeof(mwSize));
					printf("press enter!");
					c = getchar( );
					  //Calculate jcS_d.
					  //(This is slow and serial! FIX!!)
					  int current_row = 0;
					  for (int k = 0; k < len; k++) {
					  	printf("k: %d \n",k);
					    if (list[k].j != current_row) {
					      printf("current_row: %d \n", current_row);

					      //This should deal with the empty columns
					      while(true){
						      if (list[k].j-1 > current_row)
						      {
						      	assert(current_row>0);
						      	jcS_d[current_row] = jcS_d[current_row-1];
						      	current_row++;
						      }
						      else{
						      	break;
						      }
						  }

					      jcS_d[current_row] = k;
					      current_row++;
					    }
					  }
					  //Add the final element
					  jcS_d[current_row] = len;
					  assert(current_row+1==N+1);

					  for (int i = 0; i < N+1; ++i)
					  {
					  	printf("%d\n", jcS_d[i]);
					  }

					  printf("press enter!");
					  c = getchar( );
					  printf("print2");
					  #pragma omp parallel for
					  for (int rowNr = 0; rowNr < N; rowNr++) {
					    struct elementGnuSort* start = list;
					    start = list+jcS_d[rowNr];

					    struct elementGnuSort* end;
					    if(rowNr == N-1){
					      end = list+len;
					    }
					    else{
					      end = list+jcS_d[rowNr+1];
					    }
					    std::sort(start, end, compare_col_2);
					  }
					  c = getchar( );
					  printf("print3");
					  //Find the doubles!
					  int* doubles;
					  doubles = (int*) calloc(N, sizeof(int));
					  c = getchar( );
					  //... and count the doubles.
					  #pragma omp parallel for schedule(dynamic)
					  for (int rowNr = 0; rowNr < N; rowNr++) {
					    int k_end = rowNr==N-1 ? len-1 : jcS_d[rowNr+1];

					    for (int k = jcS_d[rowNr]; k < k_end; k++) {
					      if(list[k].i == list[k+1].i && list[k].j == list[k+1].j){
					        doubles[rowNr]++;
					      }
					    }
					  }
					  //printf("print4");
					  int doubles_sum = sumGnuSort(doubles, N);
					  //Calculate the culmative sum of doubles.
					  //(This is very serial. Shit.)
					  int cumsum = 0;
					  for (int k = 0; k < N; k++) {
					    cumsum += doubles[k];
					    doubles[k] = cumsum;
					  }

					  //printf("print5");
					  //Copy the values.
					  mwSize* jcS_d_old;
					  jcS_d_old = (mwSize*) mxMalloc(sizeof(jcS_d));
					  #pragma omp parallel for schedule(dynamic)
					  for (int k = 0; k < N; k++) {
					    jcS_d_old[k] = jcS_d[k];
					  }
					  //printf("print6");
					  //Remove a lot of the doubles from jcS_d
					  //#pragma omp parallel for
					  for (int rowNr = 0; rowNr < N-1; rowNr++) {
					    jcS_d[rowNr+1] -= doubles[rowNr];
					  }

					  /*for (size_t k = 0; k < N; k++) {
					    //printf("doubles[] = %d \n", doubles[k] );
					  }*/

					  double* val;
					  mwSize* ii;
					  //printf("sum(doubles,N) = %d \n",sumGnuSort(doubles,N) );
					  val = (double*) mxMalloc( (len - doubles_sum ) * sizeof(double) );
					  ii = (mwSize*) mxMalloc( (len - doubles_sum ) * sizeof(mwSize) );

					  //Sum and copy the doubles
					  //int counter = 0;
					  #pragma omp parallel for schedule(dynamic)
					  for (int rowNr = 0; rowNr < N; rowNr++) {
					    int counter = jcS_d[rowNr];
					    //printf("counter = %d, jcS_d[%d] = %d \n", counter, rowNr, jcS_d[rowNr] );
					    //assert(counter ==  jcS_d[rowNr] );
					    int k_end = rowNr==N-1 ? len-1 : jcS_d_old[rowNr+1];
					    for (int k = jcS_d_old[rowNr]; k < k_end; k++) {
					      //printf(" (k=%d, val[k]=%f ) ",k, list[k].val);
					      if(list[k].i == list[k+1].i && list[k].j == list[k+1].j){
					        //printf(" %f + %f = ", list[k].val, list[k+1].val);
					        (list[k+1].val) += list[k].val;
					        //printf("%f \n", list[k+1].val);

					      }
					      else{
					        //The -1 is to turn in into 0-index.
					        ii[counter] = list[k].i-1;
					        val[counter] = list[k].val;
					        counter++;
					        //printf("counter = %d \n", counter );
					      }
					      //The last element is forgotten!
					      if(k_end == len -1 && k==k_end-1){
					      	//The -1 is to turn in into 0-index.
					        ii[counter] = list[k+1].i-1;
					        val[counter] = list[k+1].val;
					        counter++;
					        jcS_d[N] = counter;
					      }
					    }
					  }

					  mxArray* S = mxCreateSparse(0,0,len-doubles_sum,mxREAL);
					  mxSetNzmax(S,len-doubles_sum);
					  mxSetM(S,M);
					  mxSetN(S,N);

					  // set the column pointer
					  mxSetJc(S,jcS_d);

					  mxSetPr(S,val);

					  mxSetIr(S,ii);

					// insert the data
					//  sparse_insert(mxGetIr(S),mxGetPr(S),mxGetPi(S),
					//		irank,0,0,ii,sr,si,smod,sdiv,len,M);

					  //mxFree(irank);


					  //free(val);
					  //free(ii);
					  
					  //free(jcS_d_old);
					  //free(list);

					  //free(jcS_d);

					  return S;




			      
			    }

			void gnuSortSimple(){
			}
















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
  //printf("This happens, yeah?\n");

#pragma omp single nowait
  //printf("Nii: %d and Mii: %d\n",Nii,Mii );
  ok1 = getix(&ii,Mii,Nii,&M,nocopyii,prhs[0]);
#pragma omp single nowait
  //printf("Njj: %d and Mjj: %d\n",Njj,Mjj );
  ok2 = getix(&jj,Mjj,Njj,&N,nocopyjj,prhs[1]);
//printf("N: %d and M: %d\n",N,M );
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
	else{
		//THIS IS NOT CORRECT AND SHOULD BE MOVED! FIX!
		plhs[0] = gnuSort(ii,jj,sr,len,M,N,Nzmax);
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
    iix = (*ix = (int*) mxMalloc(M*N*sizeof(int)));
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
    mwSize *jcS = (mwSize*) mxCalloc(N+1,sizeof(mwSize));
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
    mwSize *jcS = (mwSize*) mxCalloc(N+1,sizeof(mwSize));
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
