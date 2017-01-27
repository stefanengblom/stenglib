function B = frepmat(A,rep)
%FREPMAT Fast replication of array.
%   FREPMAT does the same job as the MATLAB-function REPMAT but runs
%   faster since the main part of the code is written in C.
%
%   B = FREPMAT(A,[M N]) creates a matrix B consisting of M-by-N
%   copies of A.
%
%   In general, B = FREPMAT(A,REP) repeats A in the first dimension
%   REP(1) times, in the second dimension REP(2) times and so on for
%   all of the sizes in REP.
%
%   Unlike REPMAT, missing sizes in REP are consistently defined to be
%   1. This means that if m is a scalar, then FREPMAT(A,m) does not
%   produce the same result as REPMAT(A,m) (which is defined to be
%   REPMAT(A,[m m])). Instead, FREPMAT(A,m) produces the same result
%   as would REPMAT(A,[m 1]). Also, the traditional syntax
%   FREPMAT(A,M,N,...) is not supported.
%
%   Examples:
%     a = frepmat(1,1000);            % same as ones(1000,1)
%     b = frepmat(int8(1),[100 100]); % same as ones(100,100,'int8')
%
%     s = sprand(10,10,0.5);
%     a1 = frepmat(s,[100 10]); [nnz(a1) nzmax(a1)]
%     a2 = repmat(s,[100 10]);  [nnz(a2) nzmax(a2)] % overallocation
%
%   See also REPMAT.

%   - For doubles (real or complex), FREPMAT is about 15-20% faster
%   than REPMAT.
%
%   - For other numerical types (single precision floats, integers,
%   characters and logicals), FREPMAT is about 50% faster than REPMAT.
%
%   - For sparse matrices (real, complex or logical), FREPMAT is about
%   85-95% faster than REPMAT and allocates much less memory.
%
%   - For other data-types, such as cell-, structure- and
%   function-arrays, FREPMAT is about 10% faster than REPMAT.

% S. Engblom 2004-10-28

if isnumeric(A) || ischar(A) || islogical(A)
  % runs fast
  B = mexfrepmat(A,rep);
else
  % general code using the traditional 'indexing' style
  rep = rep(:);
  lenrep = size(rep,1);
  len = max(lenrep,tndims(A));
  sizA = int32(tsize(A,1:len));
  rep = [rep; ones(len-lenrep,1)];

  % while building the indices the speed of mexfrepmat is exploited
  ix = cell(1,len);
  for i = 1:len
    ix{i} = mexfrepmat([1:sizA(i)]',rep(i));
  end
  B = A(ix{:});
end
