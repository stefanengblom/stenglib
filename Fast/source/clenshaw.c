/* clenshaw.c */
/* S. Engblom 2007-04-19 (Revision) */
/* S. Engblom 2005-07-28 */

#include <string.h>

#include "mex.h"
#include "matrix.h"

/* forward declarations */
void SizesAndType(const mxArray **A,
		  int *dimM,int *dimN,
		  mxComplexity *type);

void recurrence(double *prY,double *piY,int M,int N,
		const mxArray *A,const mxArray *B);
void clenshaw(double *prY,double *piY,int M,int N,int KC,int NC,
	      const mxArray *A,const mxArray *B,
	      const double *prC,const double *piC);

#define ISDOUBLEMATRIX(A) (mxIsDouble(A) && !mxIsSparse(A) && \
			   mxGetNumberOfDimensions(A) == 2)
#define ISDOUBLETENSOR(A) (mxIsDouble(A) && !mxIsSparse(A))

/*------------------------------------------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  /* syntax */
  if (nrhs < 2 || 3 < nrhs || nlhs > 1)
    mexErrMsgIdAndTxt("clenshaw:e1","Expecting 2 or 3 inputs "
		      "and one output.");

  /* 'recurrence'-syntax */
  if (nrhs == 2) {
    /* chech input and allocate output */
    int dimM,dimN;
    mxComplexity type;
    SizesAndType(prhs,&dimM,&dimN,&type);

    /* allocate and evaluate recurrence */
    plhs[0] = mxCreateDoubleMatrix(dimM,dimN,type);
    recurrence(mxGetPr(plhs[0]),mxGetPi(plhs[0]),
	       dimM,dimN,prhs[0],prhs[1]);
  }
  /* 'coefficient'-syntax (Clenshaw summation) */
  else {
    /* check input and allocate output */
    int dimM,dimN,KC,NC;
    mxComplexity type;
    SizesAndType(prhs,&dimM,&dimN,&type);

    /* coefficient matrix */
    if (!ISDOUBLETENSOR(prhs[2]))
      mexErrMsgIdAndTxt("clenshaw:e2",
			"Expecting a double, non-sparse matrix.");
    if (mxGetM(prhs[2]) != dimM)
      mexErrMsgIdAndTxt("clenshaw:e4",
			"The dimensions of the coefficients must "
			"match the length of each recursion.");

    if (mxGetNumberOfDimensions(prhs[2]) > 3)
      mexErrMsgIdAndTxt("clenshaw:e5","Coefficient array can have at "
			"most three dimensions.");
    else if (mxGetNumberOfDimensions(prhs[2]) == 3) {
      KC = mxGetDimensions(prhs[2])[1];
      NC = mxGetDimensions(prhs[2])[2];
    }
    else {
      KC = mxGetN(prhs[2]);
      NC = 1;
    }

    if (dimN == 1)
      dimN = NC;
    else if (dimN != NC && NC != 1)
      mexErrMsgIdAndTxt("clenshaw:e6","The number of recursions "
			"must agree.");

    if (mxIsComplex(prhs[2])) type = mxCOMPLEX;

    /* allocate and evaluate sum */
    plhs[0] = mxCreateDoubleMatrix(KC,dimN,type);
    clenshaw(mxGetPr(plhs[0]),mxGetPi(plhs[0]),
	     dimM,dimN,KC,NC,
	     prhs[0],prhs[1],mxGetPr(prhs[2]),mxGetPi(prhs[2]));
  }
}
/*------------------------------------------------------------------------*/
void SizesAndType(const mxArray **A,
		  int *dimM,int *dimN,
		  mxComplexity *type)
/* Given the arrays A[0] and A[1], this routine determines the
   dimensions dimM and dimN of the recursions defined by this
   input. The common type (mxREAL/mxCOMPLEX) is returned in type. */
{
  *dimM = *dimN = 1;
  *type = mxREAL;

  for (int j = 0; j < 2; j++) {
    const mxArray *Aj = A[j];
    if (!ISDOUBLEMATRIX(Aj))
      mexErrMsgIdAndTxt("clenshaw:e2",
			"Expecting a double, non-sparse matrix.");
    if (*dimM == 1)
      *dimM = mxGetM(Aj);
    else if (*dimM != mxGetM(Aj) && mxGetM(Aj) != 1)
      mexErrMsgIdAndTxt("clenshaw:e3",
			"Dimensions must either match or be singletons.");
    if (*dimN == 1)
      *dimN = mxGetN(Aj);
    else if (*dimN != mxGetN(Aj) && mxGetN(Aj) != 1)
      mexErrMsgIdAndTxt("clenshaw:e3",
			"Dimensions must either match or be singletons.");
    if (mxIsComplex(Aj)) *type = mxCOMPLEX;
  }
}
/*------------------------------------------------------------------------*/
void recurrence(double *prY,double *piY,int M,int N,
		const mxArray *A,const mxArray *B)
/* Evaluates the recurrence as defined by the matrices A and B and
   stores the result in (prY,piY), which must be allocated prior to
   call. [M,N] = max(size(A),size(B)). All dimensions must either
   match or be singletons. */
{
  const double *prA = mxGetPr(A),*piA = mxGetPi(A);
  const double *prB = mxGetPr(B),*piB = mxGetPi(B);

  /* increments for columns and rows */
  const int ia = mxGetM(A) == M;
  const int ja = (ia ^ mxGetN(A) == N)*((mxGetN(A) == N)-mxGetM(A));
  const int ib = mxGetM(B) == M;
  const int jb = (ib ^ mxGetN(B) == N)*((mxGetN(B) == N)-mxGetM(B));

  /* evaluate recurrence */
  if (piY == NULL)
    for (int j = 0; j < N; j++) {
      double r1 = 1.0,r2 = 1.0,r3;
      for (int i = 0; i < M; i++) {
	r3 = r2;
	r2 = r1;
	r1 = *prA*r2+*prB*r3;
	*prY++ = r1;
	prA += ia;
	prB += ib;
      }
      prA += ja;
      prB += jb;
    }
  /* the remaining 3 complex case */
  else if (piB == NULL)
    for (int j = 0; j < N; j++) {
      double r1 = 1.0,r2 = 1.0,r3;
      double i1 = 0.0,i2 = 0.0,i3;
      for (int i = 0; i < M; i++) {
	r3 = r2; i3 = i2;
	r2 = r1; i2 = i1;
	r1 = *prA*r2-*piA*i2+*prB*r3;
	i1 = *prA*i2+*piA*r2+*prB*i3;
	*prY++ = r1; *piY++ = i1;
	prA += ia; piA += ia;
	prB += ib;
      }
      prA += ja; piA += ja;
      prB += jb;
    }
  else if (piA == NULL)
    for (int j = 0; j < N; j++) {
      double r1 = 1.0,r2 = 1.0,r3;
      double i1 = 0.0,i2 = 0.0,i3;
      for (int i = 0; i < M; i++) {
	r3 = r2; i3 = i2;
	r2 = r1; i2 = i1;
	r1 = *prA*r2+*prB*r3-*piB*i3;
	i1 = *prA*i2+*prB*i3+*piB*r3;
	*prY++ = r1; *piY++ = i1;
	prA += ia;
	prB += ib; piB += ib;
      }
      prA += ja;
      prB += jb; piB += jb;
    }
  else
    for (int j = 0; j < N; j++) {
      double r1 = 1.0,r2 = 1.0,r3;
      double i1 = 0.0,i2 = 0.0,i3;
      for (int i = 0; i < M; i++) {
	r3 = r2; i3 = i2;
	r2 = r1; i2 = i1;
	r1 = *prA*r2-*piA*i2+*prB*r3-*piB*i3;
	i1 = *prA*i2+*piA*r2+*prB*i3+*piB*r3;
	*prY++ = r1; *piY++ = i1;
	prA += ia; piA += ia;
	prB += ib; piB += ib;
      }
      prA += ja; piA += ja;
      prB += jb; piB += jb;
    }
}
/*------------------------------------------------------------------------*/
void clenshaw(double *prY,double *piY,int M,int N,int KC,int NC,
	      const mxArray *A,const mxArray *B,
	      const double *prC,const double *piC)
/* Evaluates the recurrence as defined by the matrices A and B,
   multiplies it by the coefficients C and stores the sum of the
   result in (prY,piY), which must be allocated prior to call. [M,N] =
   max(size(A),size(B),size(C,[1 3])) and C is M-by-KC-by-NC (where NC
   is either 1 or N). The result has the dimensions KC-by-N. All
   dimensions must either match or be singletons, except for C's first
   dimension which must be M. */
{
  const double *prA = mxGetPr(A),*piA = mxGetPi(A);
  const double *prB = mxGetPr(B),*piB = mxGetPi(B);

  /* start at the end since Clenshaws algorithm runs backwards */
  if (M > 2) {
    prA += mxGetNumberOfElements(A)-1;
    prB += mxGetNumberOfElements(B)-1;
    prC += M*KC*NC-1;
    prY += KC*N-1;
  }

  /* increments for columns and rows */
  const int ia = mxGetM(A) == M;
  const int ja = mxGetM(A)*(mxGetN(A) == N);
  const int ka = ia-ia*M;
  const int ib = mxGetM(B) == M;
  const int jb = mxGetM(B)*(mxGetN(B) == N);
  const int kb = ib-ib*M;
  const int jc = -(NC != N)*M*KC;

  /* all real case */
  if (piY == NULL) {
    /* two special cases: sum includes boundary values only */
    if (M == 1)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  const double R1 = prA[0]+prB[0];
	  *prY++ = prC[0]*R1; prC++;
	}
	prA += ka+ja;
	prB += kb+jb;
	prC += jc;
      }
    else if (M == 2)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  const double R2 = prA[0]+prB[0];
	  const double R1 = prA[ia]*R2+prB[ib];
	  *prY++ = prC[0]*R2+prC[1]*R1; prC += 2;
	}
	prA += ia+ka+ja;
	prB += ib+kb+jb;
	prC += jc;
      }
    /* evaluate using Clenshaw summation */ 
    else if (M >= 3)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  double r1 = *prC,r2,r3;
	  prC--;

	  /* initial boundary */
	  r2 = r1;
	  r1 = *prA*r2+*prC;
	  prA -= ia;
	  prC--;

	  /* internal formula */
	  for (int i = 0; i < M-3; i++) {
	    r3 = r2;
	    r2 = r1;
	    r1 = *prA*r2+*prB*r3+*prC;
	    prA -= ia;
	    prB -= ib;
	    prC--;
	  }

	  /* final boundary */
	  r3 = r2;
	  r2 = r1;
	  r1 = *prB*r3+*prC;

	  /* initial values of the recurrence */
	  prA -= ia;
	  prB -= 2*ib;
	  const double R2 = prA[0]+prB[0];
	  const double R1 = prA[ia]*R2+prB[ib];

	  /* accumulated sum */
	  *prY-- = R2*r1+R1*r2;

	  prA -= ka;
	  prB -= kb;
	  prC--;
	}
	prA -= ja;
	prB -= jb;
	prC -= jc;
      }
  }
  else {
    /* in order to avoid a messy code the remaining 7 complex cases
       are handled in the same way as the pure complex case */
    const int ia2 = piA == NULL ? 0 : ia;
    const int ja2 = piA == NULL ? 0 : ja;
    const int ka2 = piA == NULL ? 0 : ka;
    const int ib2 = piB == NULL ? 0 : ib;
    const int jb2 = piB == NULL ? 0 : jb;
    const int kb2 = piB == NULL ? 0 : kb;
    const int ic2 = piC == NULL ? 0 : 1;
    const int jc2 = piC == NULL ? 0 : jc;

    if (M > 2) {
      if (piA != NULL) piA += mxGetNumberOfElements(A)-1;
      if (piB != NULL) piB += mxGetNumberOfElements(B)-1;
      if (piC != NULL) piC += M*KC*NC-1;
      piY += KC*N-1;
    }

    const double zero = 0.0;
    if (piA == NULL) piA = &zero;
    if (piB == NULL) piB = &zero;
    if (piC == NULL) piC = &zero;

    if (M == 1)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  const double R1 = prA[0]+prB[0],I1 = piA[0]+piB[0];
	  *prY++ = prC[0]*R1-piC[0]*I1;
	  *piY++ = prC[0]*I1+piC[0]*R1;
	  prC++; piC += ic2;
	}
	prA += ka+ja; piA += ka2+ja2;
	prB += kb+jb; piB += kb2+jb2;
	prC += jc; piC += jc2;
      }
    else if (M == 2)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  const double R2 = prA[0]+prB[0],I2 = piA[0]+piB[0];
	  const double R1 = prA[ia]*R2-piA[ia2]*I2+prB[ib],
	    I1 = prA[ia]*I2+piA[ia2]*R2+piB[ib2];
	  *prY++ = prC[0]*R2-piC[0]*I2+prC[1]*R1-piC[ic2]*I1;
	  *piY++ = prC[0]*I2+piC[0]*R2+prC[1]*I1+piC[ic2]*R1;
	  prC += 2; piC += 2*ic2;
	}
	prA += ia+ka+ja; piA += ia2+ka2+ja2;
	prB += ib+kb+jb; piB += ib2+kb2+jb2;
	prC += jc; piC += jc2;
      }
    else if (M >= 3)
      for (int j = 0; j < N; j++) {
	for (int k = 0; k < KC; k++) {
	  double r1 = *prC,r2,r3;
	  double i1 = *piC,i2,i3;
	  prC--; piC -= ic2;
	  r2 = r1; i2 = i1;
	  r1 = *prA*r2-*piA*i2+*prC;
	  i1 = *prA*i2+*piA*r2+*piC;
	  prA -= ia; piA -= ia2;
	  prC--; piC -= ic2;
	  for (int i = 0; i < M-3; i++) {
	    r3 = r2; i3 = i2;
	    r2 = r1; i2 = i1;
	    r1 = *prA*r2-*piA*i2+*prB*r3-*piB*i3+*prC;
	    i1 = *prA*i2+*piA*r2+*prB*i3+*piB*r3+*piC;
	    prA -= ia; piA -= ia2;
	    prB -= ib; piB -= ib2;
	    prC--; piC -= ic2;
	  }
	  r3 = r2; i3 = i2;
	  r2 = r1; i2 = i1;
	  r1 = *prB*r3-*piB*i3+*prC;
	  i1 = *prB*i3+*piB*r3+*piC;
	  prA -= ia; piA -= ia2;
	  prB -= 2*ib; piB -= 2*ib2;
	  const double R2 = prA[0]+prB[0],I2 = piA[0]+piB[0];
	  const double R1 = prA[ia]*R2-piA[ia2]*I2+prB[ib],
	    I1 = prA[ia]*I2+piA[ia2]*R2+piB[ib2];
	  *prY-- = R2*r1-I2*i1+R1*r2-I1*i2;
	  *piY-- = R2*i1+I2*r1+I1*r2+R1*i2;
	  prA -= ka; piA -= ka2;
	  prB -= kb; piB -= kb2;
	  prC--; piC -= ic2;
	}
	prA -= ja; piA -= ja2;
	prB -= jb; piB -= jb2;
	prC -= jc; piC -= jc2;
      }
  }
}
/*------------------------------------------------------------------------*/
