/* fsetop.c */
/* S. Engblom 2010-02-09 (Revision, extended union operation) */
/* S. Engblom 2006-06-13 (Major revision) */
/* S. Engblom 2005-06-02 */

#include <math.h>
// temporary fix for CC under Solaris:
#ifndef NO_STDINT
#include <stdint.h>
#endif
#include <string.h>

#include "mex.h"
#include "matrix.h"

/*------------------------------------------------------------------------*/

// forward declarations
typedef enum {CHECK = 8,INTERSECT = 4,SETDIFF = 6,SETXOR = 1,
	      UNION = 0,UNIQUE = 7,ISMEMBER = 10,UNKNOWN = -1} HASH_OP;
HASH_OP getop(const char *str);

// hash-table for two arrays A and B containing data
typedef struct {
  uint8_t **hash_tb;              // pointers to raw data
  size_t hashsz;                  // size of table
  char *hash_id;                  // identification of element
  uint32_t *hash_ixa,*hash_ixb;   // pointers from data into table
  size_t na,nb,nx;                // number of A's, B's and X's
} hashTable;

// hash-table for two arrays containing mxArrays
typedef struct {
  mxArray **hash_tb;
  size_t hashsz;
  char *hash_id;
  uint32_t *hash_ixa,*hash_ixb;
  uint32_t *hash_iixa,*hash_iixb; // pointers back into hash_ixa, hash_ixb
  size_t na,nb,nx;
} mxhashTable;

// hash-utility functions
size_t hashsize(size_t len);
uint32_t hash(const uint8_t *val,size_t nbytes,uint32_t offset);
uint32_t hash2nd(const uint8_t *val,size_t nbytes);
bool mx_IsEqual(const mxArray *A,const mxArray *B);

// Venn-topology routines
void hashSort(hashTable *T,
	      const uint8_t *A,size_t Na,
	      const uint8_t *B,size_t Nb,size_t Mbytes);
void mxhashSort(mxhashTable *T,
		const mxArray *A,size_t Na,
		const mxArray *B,size_t Nb);

// post-processing of hash-tables
typedef void (*fsetopOut)(hashTable,size_t,size_t,size_t,
			  int,mxArray **,size_t,mxClassID);

void intersectOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		  int nlhs,mxArray *plhs[],
		  size_t,mxClassID);
void setdiffOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		int nlhs,mxArray *plhs[],
		size_t,mxClassID);
void setxorOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	       int nlhs,mxArray *plhs[],
	       size_t,mxClassID);
void unionOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	      int nlhs,mxArray *plhs[],
	      size_t,mxClassID);
void uniqueOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	       int nlhs,mxArray *plhs[],
	       size_t,mxClassID);
void ismemberOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		 int nlhs,mxArray *plhs[],
		 size_t,mxClassID);

// mx-versions
typedef void (*mxfsetopOut)(mxhashTable,size_t,size_t,
			    int,mxArray **);

void mxintersectOut(mxhashTable T,size_t Na,size_t Nb,
		    int nlhs,mxArray *plhs[]);
void mxsetdiffOut(mxhashTable T,size_t Na,size_t Nb,
		  int nlhs,mxArray *plhs[]);
void mxsetxorOut(mxhashTable T,size_t Na,size_t Nb,
		 int nlhs,mxArray *plhs[]);
void mxunionOut(mxhashTable T,size_t Na,size_t Nb,
		int nlhs,mxArray *plhs[]);
void mxuniqueOut(mxhashTable T,size_t Na,size_t Nb,
		 int nlhs,mxArray *plhs[]);
void mxismemberOut(mxhashTable T,size_t Na,size_t Nb,
		   int nlhs,mxArray *plhs[]);

/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  // fetch first argument op (read at most 9 characters)
  char strop[10];
  HASH_OP op; 

  // determine the type of the set-operation
  if (nrhs < 2)
    mexErrMsgIdAndTxt("fsetop:e1","Expecting at least two inputs.");
  if (!mxIsChar(prhs[0]))
    mexErrMsgIdAndTxt("fsetop:e2","First argument must be a "
		      "character array.");
  if (mxGetString(prhs[0],strop,10) != 0 || (op = getop(strop)) == UNKNOWN)
    mexErrMsgIdAndTxt("fsetop:e3","Unknown set operation.");

  // syntax
  if (op != CHECK && op != UNIQUE) {
    if (nrhs != 3)
      mexErrMsgIdAndTxt("fsetop:e5","Expecting exactly three inputs.");
  }
  else if (nrhs != 2)
    mexErrMsgIdAndTxt("fsetop:e8","Expecting exactly two inputs.");
  if (op == UNION) {
    if (nlhs > 5)
      mexErrMsgIdAndTxt("fsetop:e13","Expecting at most five outputs.");
  }
  else if (op != CHECK && op != SETDIFF && op != ISMEMBER) {
    if (nlhs > 3)
      mexErrMsgIdAndTxt("fsetop:e9","Expecting at most three outputs.");
  }
  else if (op != CHECK) {
    if (nlhs > 2)
      mexErrMsgIdAndTxt("fsetop:e10","Expecting at most two outputs.");
  }
  else if (nlhs > 1)
    mexErrMsgIdAndTxt("fsetop:e12","Expecting at most one output."); 

  // two cases
  const bool isCell = mxIsCell(prhs[1]);

  // case 1: input is numeric array
  if (!isCell) {
    uint8_t *A,*B = NULL;
    size_t Na,Nb = 0,Mbytes;

    // fetch second input A
    if (!(mxIsNumeric(prhs[1]) || mxIsChar(prhs[1]) || mxIsLogical(prhs[1])) 
	|| mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) ||
	mxGetNumberOfDimensions(prhs[1]) > 2)
      mexErrMsgIdAndTxt("fsetop:e4","Input argument must be a 2-D real "
			"and non-sparse array.");
    A = (uint8_t *)mxGetData(prhs[1]);

    Mbytes = mxGetM(prhs[1])*mxGetElementSize(prhs[1]);
    Na = mxGetN(prhs[1]);

    // quick return
    if (op == CHECK) {
      uint32_t *C;
      plhs[0] = mxCreateNumericMatrix(1,Na,mxUINT32_CLASS,mxREAL);
      C = mxGetData(plhs[0]);
      for (size_t i = 0; i < Na; i++)
	C[i] = hash(&A[i*Mbytes],Mbytes,0);
      return;
    }

    // fetch third input B, if any
    if (op != UNIQUE) {
      if (mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]) ||
	  mxGetNumberOfDimensions(prhs[2]) > 2)
	mexErrMsgIdAndTxt("fsetop:e4","Input argument must be a 2-D real "
			  "and non-sparse array.");
      if (mxGetClassID(prhs[1]) != mxGetClassID(prhs[2]))
	mexErrMsgIdAndTxt("fsetop:e6","Arguments must be "
			  "of the same class.");
      if (mxGetM(prhs[1]) != mxGetM(prhs[2]))
	mexErrMsgIdAndTxt("fsetop:e7","The number of rows must match.");
 
      B = (uint8_t *)mxGetData(prhs[2]);
      Nb = mxGetN(prhs[2]);
    }

    // sort it out!
    hashTable T;
    if (op != ISMEMBER)
      hashSort(&T,A,Na,B,Nb,Mbytes);
    else
      hashSort(&T,B,Nb,A,Na,Mbytes); // 2nd input must be the master array

    // produce whatever output is asked for
    static const fsetopOut ftab[] = {
      &unionOut,&setxorOut,NULL,NULL,&intersectOut,NULL,
      &setdiffOut,&uniqueOut,NULL,NULL,&ismemberOut
    };

    (*ftab[op])(T,Na,Nb,Mbytes,nlhs,plhs,
		mxGetElementSize(prhs[1]),mxGetClassID(prhs[1]));

    // dealllocate the hash table
    mxFree(T.hash_ixb);
    mxFree(T.hash_ixa);
    mxFree(T.hash_id);
    mxFree(T.hash_tb);
  }
  // case 2: cell-array of arrays containing real and non-sparse data
  else {
    mxArray *A,*B = NULL;
    size_t Na,Nb = 0;

    // fetch second input A
    A = (mxArray *)prhs[1];
    Na = mxGetNumberOfElements(A);
    for (size_t i = 0; i < Na; i++) {
      const mxArray *Ai = mxGetCell(A,i);
      if (!(mxIsNumeric(Ai) || mxIsChar(Ai) || mxIsLogical(Ai)) ||
	  mxIsComplex(Ai) || mxIsSparse(Ai))
	mexErrMsgIdAndTxt("fsetop:e11","Cell-vector must contain real, "
			  "non-sparse arrays.");
    }

    // quick return
    if (op == CHECK) {
      uint32_t *C;
      plhs[0] = mxCreateNumericMatrix(1,Na,mxUINT32_CLASS,mxREAL);
      C = mxGetData(plhs[0]);
      for (size_t i = 0; i < Na; i++) {
	const mxArray *Ai = mxGetCell(A,i);
	uint32_t hashkey =  hash(mxGetData(Ai),mxGetNumberOfElements(Ai)*
				 mxGetElementSize(Ai),mxGetClassID(Ai));
	hashkey = hash((uint8_t *)mxGetDimensions(Ai),
		       mxGetNumberOfDimensions(Ai)*sizeof(int),hashkey);
	C[i] = hashkey;
      }
      return;
    }

    // fetch third input B, if any
    if (op != UNIQUE) {
      if (!mxIsCell(prhs[2]))
	mexErrMsgIdAndTxt("fsetop:e6","Arguments must be "
			  "of the same class.");

      B = (mxArray *)prhs[2];
      Nb = mxGetNumberOfElements(B);
      for (size_t i = 0; i < Nb; i++) {
	const mxArray *Bi = mxGetCell(B,i);
	if (!(mxIsNumeric(Bi) || mxIsChar(Bi) || mxIsLogical(Bi)) ||
	    mxIsComplex(Bi) || mxIsSparse(Bi))
	  mexErrMsgIdAndTxt("fsetop:e11","Cell-vector must contain real, "
			    "non-sparse arrays.");
      }
    }

    // sort it out!
    mxhashTable T;
    if (op != ISMEMBER)
      mxhashSort(&T,A,Na,B,Nb);
    else
      mxhashSort(&T,B,Nb,A,Na); // 2nd input must be the master array

    // produce whatever output is asked for
    static const mxfsetopOut ftab[] = {
      &mxunionOut,&mxsetxorOut,NULL,NULL,&mxintersectOut,NULL,
      &mxsetdiffOut,&mxuniqueOut,NULL,NULL,&mxismemberOut
    };

    (*ftab[op])(T,Na,Nb,nlhs,plhs);

    // deallocate the hash table
    mxFree(T.hash_iixb);
    mxFree(T.hash_iixa);
    mxFree(T.hash_ixb);
    mxFree(T.hash_ixa);
    mxFree(T.hash_id);
    mxFree(T.hash_tb);
  }
}
/*------------------------------------------------------------------------*/
void intersectOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		  int nlhs,mxArray *plhs[],
		  size_t siz,mxClassID id)
/* [C,IA,IB] = FSETOP('intersect',A,B) returns the columns common to
   both A and B. C = A(:,IA) = B(:,IB). */
{
  uint8_t *C;
  double *ia = NULL,*ib = NULL;

  // allocate output
  plhs[0] = mxCreateNumericMatrix(Mbytes/siz,T.nx,id,mxREAL);
  C = (uint8_t *)mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.nx,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nx,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'X') {
      memcpy(&C[j*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ia != NULL) {
	ia[j] = i+1;
	if (ib != NULL)
	  // IB follows the order determined by IA
	  T.hash_tb[hashkey] = &C[j*Mbytes];
      }

      // first unique element only
      T.hash_id[hashkey] = 'Y';
      j++;
    }
  }

  if (ib != NULL)
    for (size_t i = 0; i < Nb; i++) {
      const uint32_t hashkey = T.hash_ixb[i];
      if (T.hash_id[hashkey] == 'Y') {
	ib[(T.hash_tb[hashkey]-&C[0])/(Mbytes ? Mbytes : 1)] = i+1;

	// again, first unique element only
	T.hash_id[hashkey] = 'Z';
      }
    }
}
/*------------------------------------------------------------------------*/
void mxintersectOut(mxhashTable T,size_t Na,size_t Nb,
		    int nlhs,mxArray *plhs[])
{
  mxArray *C;
  double *ia = NULL,*ib = NULL;

  C = plhs[0] = mxCreateCellMatrix(1,T.nx);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.nx,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nx,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'X') {
      mxSetCell(C,j,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ia != NULL) {
	ia[j] = i+1;
	if (ib != NULL)
	  ib[j] = T.hash_iixb[hashkey]+1; // new construction
      }

      T.hash_id[hashkey] = 'Y';
      j++;
    }
  }
}
/*------------------------------------------------------------------------*/
void setdiffOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		int nlhs,mxArray *plhs[],
		size_t siz,mxClassID id)
/* [C,IA] = FSETOP('setdiff',A,B) returns the columns in A that are
   not in B. C = A(:,IA). */
{
  uint8_t *C;
  double *ia = NULL;

  // allocate output
  plhs[0] = mxCreateNumericMatrix(Mbytes/siz,T.na,id,mxREAL);
  C = (uint8_t *)mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
  }

  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A') {
      memcpy(&C[j*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ia != NULL) ia[j] = i+1;

      // first unique element only
      T.hash_id[hashkey] = 'Y';
      j++;
    }
  }
}
/*------------------------------------------------------------------------*/
void mxsetdiffOut(mxhashTable T,size_t Na,size_t Nb,
		  int nlhs,mxArray *plhs[])
{
  mxArray *C;
  double *ia = NULL;

  C = plhs[0] = mxCreateCellMatrix(1,T.na);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
  }

  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A') {
      mxSetCell(C,j,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ia != NULL) ia[j] = i+1;

      T.hash_id[hashkey] = 'Y';
      j++;
    }
  }
}
/*------------------------------------------------------------------------*/
void setxorOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	       int nlhs,mxArray *plhs[],
	       size_t siz,mxClassID id)
/* [C,IA,IB] = FSETOP('setxor',A,B) returns the columns that are not
   in the intersection of A and B. C = [A(:,IA) B(:,IB)]. */
{
  uint8_t *C;
  double *ia = NULL,*ib = NULL;

  // allocate output
  plhs[0] = mxCreateNumericMatrix(Mbytes/siz,T.na+T.nb,id,mxREAL);
  C = (uint8_t *)mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nb,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  size_t k = 0;
  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A') {
      memcpy(&C[k++*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ia != NULL) ia[j++] = i+1;

      // first unique element only
      T.hash_id[hashkey] = 'Y';
    }
  }

  for (size_t i = 0, j = 0; i < Nb; i++) {
    const uint32_t hashkey = T.hash_ixb[i];
    if (T.hash_id[hashkey] == 'B') {
      memcpy(&C[k++*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ib != NULL) ib[j++] = i+1;

      T.hash_id[hashkey] = 'Y';
    }
  }
}
/*------------------------------------------------------------------------*/
void mxsetxorOut(mxhashTable T,size_t Na,size_t Nb,
		 int nlhs,mxArray *plhs[])
{
  mxArray *C;
  double *ia = NULL,*ib = NULL;

  C = plhs[0] = mxCreateCellMatrix(1,T.na+T.nb);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nb,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  size_t k = 0;
  for (size_t i = 0,j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A') {
      mxSetCell(C,k++,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ia != NULL) ia[j++] = i+1;

      T.hash_id[hashkey] = 'Y';
    }
  }

  for (size_t i = 0, j = 0; i < Nb; i++) {
    const uint32_t hashkey = T.hash_ixb[i];
    if (T.hash_id[hashkey] == 'B') {
      mxSetCell(C,k++,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ib != NULL) ib[j++] = i+1;

      T.hash_id[hashkey] = 'Y';
    }
  }
}
/*------------------------------------------------------------------------*/
void unionOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	      int nlhs,mxArray *plhs[],
	      size_t siz,mxClassID id)
/* [C,IA,IB,JA,JB] = FSETOP('union',A,B) returns the combined columns
   from A and B but with no repetitions. C = [A(:,IA) B(:,IB)] and A =
   C(:,JA), B = C(:,JB). */
{
  uint8_t *C;
  double *ia = NULL,*ib = NULL,*ja = NULL,*jb = NULL;

  plhs[0] = mxCreateNumericMatrix(Mbytes/siz,T.na+T.nx+T.nb,id,mxREAL);
  C = (uint8_t *)mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na+T.nx,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nb,mxREAL);
      ib = mxGetPr(plhs[2]);
      if (nlhs > 3) {
	plhs[3] = mxCreateDoubleMatrix(1,Na,mxREAL);
	ja = mxGetPr(plhs[3]);
	if (nlhs > 4) {
	  plhs[4] = mxCreateDoubleMatrix(1,Nb,mxREAL);
	  jb = mxGetPr(plhs[4]);
	}
      }
    }
  }

  size_t k = 0;
  for (size_t i = 0, j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A' || T.hash_id[hashkey] == 'X') {
      memcpy(&C[k*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ia != NULL) {
	ia[j++] = i+1;
	if (ja != NULL) {
	  ja[i] = k+1;
	  T.hash_tb[hashkey] = &C[k*Mbytes];
	}
      }
      k++;

      // first unique element
      T.hash_id[hashkey] = 'Y';
    }
    else if (ja != NULL && T.hash_id[hashkey] == 'Y')
      ja[i] = (T.hash_tb[hashkey]-&C[0])/(Mbytes ? Mbytes : 1)+1;
  }

  for (size_t i = 0, j = 0; i < Nb; i++) {
    const uint32_t hashkey = T.hash_ixb[i];
    if (T.hash_id[hashkey] == 'B') {
      memcpy(&C[k*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ib != NULL) {
	ib[j++] = i+1;
	if (jb != NULL) {
	  jb[i] = k+1;
	  T.hash_tb[hashkey] = &C[k*Mbytes];
	}
      }
      k++;

      T.hash_id[hashkey] = 'Y';
    }
    else if (jb != NULL && T.hash_id[hashkey] == 'Y')
      jb[i] = (T.hash_tb[hashkey]-&C[0])/(Mbytes ? Mbytes : 1)+1;
  }
}
/*------------------------------------------------------------------------*/
void mxunionOut(mxhashTable T,size_t Na,size_t Nb,
		int nlhs,mxArray *plhs[])
{
  mxArray *C;
  double *ia = NULL,*ib = NULL,*ja = NULL,*jb = NULL;

  C = plhs[0] = mxCreateCellMatrix(1,T.na+T.nx+T.nb);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na+T.nx,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,T.nb,mxREAL);
      ib = mxGetPr(plhs[2]);
      if (nlhs > 3) {
	plhs[3] = mxCreateDoubleMatrix(1,Na,mxREAL);
	ja = mxGetPr(plhs[3]);
	if (nlhs > 4) {
	  plhs[4] = mxCreateDoubleMatrix(1,Nb,mxREAL);
	  jb = mxGetPr(plhs[4]);
	}
      }
    }
  }

  size_t k = 0;
  for (size_t i = 0, j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];
    if (T.hash_id[hashkey] == 'A' || T.hash_id[hashkey] == 'X') {
      mxSetCell(C,k,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ia != NULL) {
	ia[j++] = i+1;
	if (ja != NULL) {
	  ja[i] = k+1;
	  T.hash_iixa[hashkey] = k;
	}
      }
      k++;

      T.hash_id[hashkey] = 'Y';
    }
    else if (ja != NULL && T.hash_id[hashkey] == 'Y')
      ja[i] = T.hash_iixa[hashkey]+1;
  }

  for (size_t i = 0, j = 0; i < Nb; i++) {
    const uint32_t hashkey = T.hash_ixb[i];
    if (T.hash_id[hashkey] == 'B') {
      mxSetCell(C,k,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ib != NULL) {
	ib[j++] = i+1;
	if (jb != NULL) {
	  jb[i] = k+1;
	  T.hash_iixb[hashkey] = k;
	}
      }
      k++;

      T.hash_id[hashkey] = 'Z';
    }
    else if (jb != NULL) {
      if (T.hash_id[hashkey] == 'Y')
	jb[i] = T.hash_iixa[hashkey]+1;
      else if (T.hash_id[hashkey] == 'Z')
	jb[i] = T.hash_iixb[hashkey]+1;
    }
  }
}
/*------------------------------------------------------------------------*/
void uniqueOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
	       int nlhs,mxArray *plhs[],
	       size_t siz,mxClassID id)
/* [B,IA,IB] = FSETOP('unique',A) returns the same columns as in A but
   with no repetitions. B = A(:,IA) and A = B(:,IB). */
{
  uint8_t *B;
  double *ia = NULL,*ib = NULL;

  plhs[0] = mxCreateNumericMatrix(Mbytes/siz,T.na,id,mxREAL);
  B = (uint8_t *)mxGetData(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,Na,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  for (size_t i = 0, j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];

    if (T.hash_id[hashkey] == 'A') {
      memcpy(&B[j*Mbytes],T.hash_tb[hashkey],Mbytes);
      if (ia != NULL) {
	ia[j] = i+1;
	if (ib != NULL)
	  // because IB points into B
	  T.hash_tb[hashkey] = &B[j*Mbytes];
      }

      // first unique element only
      T.hash_id[hashkey] = 'Y';
      j++;
    }

    // IB does not refer to unique elements
    if (ib != NULL) ib[i] = ((T.hash_tb[hashkey]-&B[0]))/
      (Mbytes ? Mbytes : 1)+1;
  }
}
/*------------------------------------------------------------------------*/
void mxuniqueOut(mxhashTable T,size_t Na,size_t Nb,
		 int nlhs,mxArray *plhs[])
{
  mxArray *B;
  double *ia = NULL,*ib = NULL;

  B = plhs[0] = mxCreateCellMatrix(1,T.na);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,T.na,mxREAL);
    ia = mxGetPr(plhs[1]);
    if (nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix(1,Na,mxREAL);
      ib = mxGetPr(plhs[2]);
    }
  }

  for (size_t i = 0, j = 0; i < Na; i++) {
    const uint32_t hashkey = T.hash_ixa[i];

    if (T.hash_id[hashkey] == 'A') {
      mxSetCell(B,j,mxDuplicateArray(T.hash_tb[hashkey]));
      if (ia != NULL) {
	ia[j] = i+1;
	if (ib != NULL)
	  T.hash_iixa[hashkey] = j;
      }

      T.hash_id[hashkey] = 'Y';
      j++;
    }

    // new construction
    if (ib != NULL) ib[i] = T.hash_iixa[hashkey]+1;
  }
}
/*------------------------------------------------------------------------*/
void ismemberOut(hashTable T,size_t Na,size_t Nb,size_t Mbytes,
		 int nlhs,mxArray *plhs[],
		 size_t siz,mxClassID id)
/* [IA,IB] = FSETOP('ismember',A,B) returns a logical vector IA
   containing 1 where the columns of A are also columns of B and 0
   otherwise. IB contains the index in B of each column in A and zero
   if no such index exists. A(:,IA) = B(:,IB(IA). */
{
  mxLogical *ia;
  double *ib = NULL;

  plhs[0] = mxCreateLogicalMatrix(1,Na);
  ia = mxGetLogicals(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,Na,mxREAL);
    ib = mxGetPr(plhs[1]);
  }

  for (size_t j = 0; j < Na; j++) {
    // the meaning of hash_ixa and hash_ixb has now switched!
    const size_t hashkey = T.hash_ixb[j];
    if (T.hash_id[hashkey] == 'X') {
      ia[j] = 1;

      if (ib != NULL)
	ib[j] = (T.hash_tb[hashkey]-T.hash_tb[T.hash_ixa[0]])/
	  (Mbytes ? Mbytes : 1)+1;
    }
  }
}
/*------------------------------------------------------------------------*/
void mxismemberOut(mxhashTable T,size_t Na,size_t Nb,
		   int nlhs,mxArray *plhs[])
{
  mxLogical *ia;
  double *ib = NULL;

  plhs[0] = mxCreateLogicalMatrix(1,Na);
  ia = mxGetLogicals(plhs[0]);
  if (nlhs > 1) {
    plhs[1] = mxCreateDoubleMatrix(1,Na,mxREAL);
    ib = mxGetPr(plhs[1]);
  }

  for (size_t j = 0; j < Na; j++) {
    const size_t hashkey = T.hash_ixb[j];
    if (T.hash_id[hashkey] == 'X') {
      ia[j] = 1;

      // new construction
      if (ib != NULL) ib[j] = T.hash_iixa[hashkey]+1;
    }
  }
}
/*------------------------------------------------------------------------*/
HASH_OP getop(const char *str)
/* Small perfect hash for the supported set operations. Note that
   str[4] must be readable for this to work. Returns UNKNOWN on
   failure. */
{
  static const char *tab[] = {
    "union","setxor","","","intersect","",
    "setdiff","unique","check","","ismember"
  };
  const unsigned hash = ((unsigned)str[4])%11;

  if (str[0] != '\0' && strcmp(str,tab[hash]) == 0)
    return hash;
  else
    return UNKNOWN;
}
/*------------------------------------------------------------------------*/

// table is only allowed to become about 70% full
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

// should work on most platforms
#define get16bits(d) (*((const uint16_t *)(d)))

uint32_t hash(const uint8_t *val,size_t nbytes,uint32_t offset)
/* Computes a hashkey from bytes. The input offset may be used to
   generate a new hashkey from the same data. Use zero as the default.

   The original author of this function is Paul Hsieh. See
   http://www.azillionmonkeys.com/qed/hash.html for further details.
*/
{
  uint32_t hash = offset,temp,rem;

  // incorporate 4 bytes at a time
  rem = nbytes&3;
  nbytes >>= 2;

  // main loop
  for ( ; nbytes > 0; nbytes--, val += 4) {
    hash += get16bits(val);
    temp = (get16bits(val+2) << 11)^hash;
    hash = (hash << 16)^temp;
    hash += hash >> 11;
  }

  // handle end cases
  switch (rem) {
  case 3:
    hash += get16bits(val);
    hash ^= hash << 16;
    hash ^= val[2] << 18;
    hash += hash >> 11;
    break;
  case 2:
    hash += get16bits(val);
    hash ^= hash << 11;
    hash += hash >> 17;
    break;
  case 1:
    hash += val[0];
    hash ^= hash << 10;
    hash += hash >> 1;
  }

  // force "avalanching" of final 127 bits
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 2;
  hash += hash >> 15;
  hash ^= hash << 10;

  return hash;
}
/*------------------------------------------------------------------------*/
uint32_t hash2nd(const uint8_t *val,size_t nbytes)
/* Computes a secondary (incremental) hashkey from bytes. */
{
  // a fairly independent and odd hashkey
  return hash(val,nbytes,0x0f0f0f0f)|1;
}
/*------------------------------------------------------------------------*/
bool mx_IsEqual(const mxArray *A,const mxArray *B)
/* Bitwise equality for real non-sparse mxArrays containing data. This
   function will not work properly for sparse arrays, cell-arrays,
   struct-arrays or for arrays containing imaginary data. */
{
  // class
  if (mxGetClassID(A) != mxGetClassID(B)) return false;

  // real data
  const int nelems = mxGetNumberOfElements(A);
  const int sizelem = mxGetElementSize(A);
  if (nelems != mxGetNumberOfElements(B)) return false;
  if (memcmp(mxGetData(A),mxGetData(B),nelems*sizelem)) return false;

  // shape
  int ndims = mxGetNumberOfDimensions(A);
  if (ndims != mxGetNumberOfDimensions(B)) return false;
  if (memcmp(mxGetDimensions(A),mxGetDimensions(B),ndims*sizeof(int)))
    return false;

  return true;
}
/*------------------------------------------------------------------------*/
void hashSort(hashTable *T,
	      const uint8_t *A,size_t Na,
	      const uint8_t *B,size_t Nb,size_t Mbytes)
/* Determines the full Venn-topology of the arrays A and B using
   hashing.

   Input are the arrays A and B of sizes Mbytes-by-Na and
   Mbytes-by-Nb.

   On return, hash_tb is the hash-table pointing back into the columns
   of A and B, hash_id contains 'A' for pointers strictly into A, 'B'
   for pointers strictly into B and 'X' for pointers into A but with
   equivalent elements contained in B as well. The number of A's, B's
   and X's are returned in na, nb and nx respectively. The vectors
   hash_ixa and hash_ixb has the same length as A and B and points
   into hash_tb/hash_id.
*/
{
  // hash table and index into hash table
  T->hashsz = hashsize(Na+Nb);
  T->hash_tb = mxMalloc(T->hashsz*sizeof(uint8_t *));
  T->hash_id = mxCalloc(T->hashsz,sizeof(char));
  T->hash_ixa = mxMalloc(Na*sizeof(uint32_t));
  T->hash_ixb = mxMalloc(Nb*sizeof(uint32_t));

  // clear counters
  T->na = T->nb = T->nx = 0;

  // hash columns of A
  for (size_t i = 0; i < Na; i++) {
    uint32_t hashkey = hash(&A[i*Mbytes],Mbytes,0)&(T->hashsz-1);

    // collision?
    if (T->hash_id[hashkey] != '\0' && 
	memcmp(&A[i*Mbytes],T->hash_tb[hashkey],Mbytes)) {
      const uint32_t hashkey2nd = hash2nd(&A[i*Mbytes],Mbytes);

      // rehash until empty slot found
      do
	hashkey = (hashkey+hashkey2nd)&(T->hashsz-1);
      while (T->hash_id[hashkey] != '\0' && 
	     memcmp(&A[i*Mbytes],T->hash_tb[hashkey],Mbytes));
    }

    // fill in table
    if (T->hash_id[hashkey] == '\0') {
      T->hash_tb[hashkey] = &((uint8_t *)A)[i*Mbytes];
      T->hash_id[hashkey] = 'A';
      T->hash_ixa[i] = hashkey;
      T->na++;
    }
    else
      T->hash_ixa[i] = hashkey;
  }

  // hash columns of B
  for (size_t i = 0; i < Nb; i++) {
    uint32_t hashkey = hash(&B[i*Mbytes],Mbytes,0)&(T->hashsz-1);

    if (T->hash_id[hashkey] != '\0' && 
	memcmp(&B[i*Mbytes],T->hash_tb[hashkey],Mbytes)) {
      const uint32_t hashkey2nd = hash2nd(&B[i*Mbytes],Mbytes);

      do
	hashkey = (hashkey+hashkey2nd)&(T->hashsz-1);
      while (T->hash_id[hashkey] != '\0' && 
	     memcmp(&B[i*Mbytes],T->hash_tb[hashkey],Mbytes));
    }

    if (T->hash_id[hashkey] == '\0') {
      T->hash_tb[hashkey] = &((uint8_t *)B)[i*Mbytes];
      T->hash_id[hashkey] = 'B';
      T->hash_ixb[i] = hashkey;
      T->nb++;
    }
    else {
      T->hash_ixb[i] = hashkey;

      // special case here
      if (T->hash_id[hashkey] == 'A') {
	T->na--;
	T->hash_id[hashkey] = 'X';
	T->nx++;
      }
    }
  }
}
/*------------------------------------------------------------------------*/
void mxhashSort(mxhashTable *T,
		const mxArray *A,size_t Na,
		const mxArray *B,size_t Nb)
/* Corresponding routine for two cell-arrays of mxArrays. */
{
  T->hashsz = hashsize(Na+Nb);
  T->hash_tb = mxMalloc(T->hashsz*sizeof(mxArray *));
  T->hash_id = mxCalloc(T->hashsz,sizeof(char));
  T->hash_ixa = mxMalloc(Na*sizeof(uint32_t));
  T->hash_ixb = mxMalloc(Nb*sizeof(uint32_t));

  // inverse pointers (back into hash_ixa, hash_ixb)
  T->hash_iixa = mxMalloc(T->hashsz*sizeof(uint32_t));
  T->hash_iixb = mxMalloc((Nb != 0)*T->hashsz*sizeof(uint32_t));

  T->na = T->nb = T->nx = 0;

  for (size_t i = 0; i < Na; i++) {
    // hashing over data and class only...
    const mxArray *Ai = mxGetCell(A,i);
    uint32_t hashkey = hash(mxGetData(Ai),mxGetNumberOfElements(Ai)*
			    mxGetElementSize(Ai),mxGetClassID(Ai));
    hashkey = hashkey&(T->hashsz-1);

    if (T->hash_id[hashkey] != '\0' && 
	!mx_IsEqual(Ai,T->hash_tb[hashkey])) {
      // ...hashing over the shape is used as the secondary key
      const uint32_t hashkey2nd = hash2nd((const uint8_t *)
					  mxGetDimensions(Ai),
					  mxGetNumberOfDimensions(Ai)*
					  sizeof(int));

      do
	hashkey = (hashkey+hashkey2nd)&(T->hashsz-1);
      while (T->hash_id[hashkey] != '\0' && 
	     !mx_IsEqual(Ai,T->hash_tb[hashkey]));
    }

    if (T->hash_id[hashkey] == '\0') {
      T->hash_tb[hashkey] = (mxArray *)Ai;
      T->hash_id[hashkey] = 'A';
      T->hash_ixa[i] = hashkey;
      T->hash_iixa[hashkey] = i;
      T->na++;
    }
    else
      T->hash_ixa[i] = hashkey;
  }

  for (size_t i = 0; i < Nb; i++) {
    const mxArray *Bi = mxGetCell(B,i);
    uint32_t hashkey = hash(mxGetData(Bi),mxGetNumberOfElements(Bi)*
			    mxGetElementSize(Bi),mxGetClassID(Bi));
    hashkey = hashkey&(T->hashsz-1);

    if (T->hash_id[hashkey] != '\0' && 
	!mx_IsEqual(Bi,T->hash_tb[hashkey])) {
      const uint32_t hashkey2nd = hash2nd((const uint8_t *)
					  mxGetDimensions(Bi),
					  mxGetNumberOfDimensions(Bi)*
					  sizeof(int));

      do
	hashkey = (hashkey+hashkey2nd)&(T->hashsz-1);
      while (T->hash_id[hashkey] != '\0' && 
	     !mx_IsEqual(Bi,T->hash_tb[hashkey]));
    }

    if (T->hash_id[hashkey] == '\0') {
      T->hash_tb[hashkey] = (mxArray *)Bi;
      T->hash_id[hashkey] = 'B';
      T->hash_ixb[i] = hashkey;
      T->hash_iixb[hashkey] = i;
      T->nb++;
    }
    else {
      T->hash_ixb[i] = hashkey;

      if (T->hash_id[hashkey] == 'A') {
	T->na--;
	T->hash_id[hashkey] = 'X';
	T->hash_iixb[hashkey] = i;
	T->nx++;
      }
    }
  }
}
/*------------------------------------------------------------------------*/
