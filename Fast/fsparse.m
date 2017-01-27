function S = fsparse(ii,jj,ss,siz,flag,nthreads)
%FSPARSE Fast assembly of sparse matrix.
%   FSPARSE does the same job as the MATLAB-function SPARSE but runs faster
%   since the assembly is based on indirect addressing and on hashing rather
%   than on sorting. FSPARSE also in general allocates the sparse matrix
%   exactly, that is, the resulting sparse matrix S satisfies NNZ(S) =
%   NZMAX(S). In addition, FSPARSE supports an efficient and powerful
%   'assembly' syntax which makes it easy to create sparse banded matrices
%   or sparse matrices arising in finite difference computations.
%
%   S = FSPARSE(X), where X is a full or sparse matrix, constructs a sparse
%   matrix S by squeezing out any zero elements.
%
%   S = FSPARSE(II,JJ,SS) creates a sparse matrix S from triplet data
%   (II,JJ,SS). The inputs are all matrices with either matching or
%   singleton dimensions according to certain rules (see below). For
%   instance, II may be 4-by-10, JJ 1-by-10 and SS 4-by-1.
%
%   The result is a sparse matrix formally satisfying the relation
%   S(II(i,j),JJ(i,j)) = SS(i,j), except of course for singleton
%   dimensions. In the example above for instance, we have the relation
%   S(II(i,j),JJ(1,j)) = SS(i,1). Repeated indices are as usual summed
%   together.
%
%   If II is IM-by-IN, JJ JM-by-JN and SS SM-by-SN, then it is required that
%   (1) IN = JN or 1, (2) JM = IM or 1, (3) SM = IM or 1 and (4) SN = JN or
%   1.
%
%   S = FSPARSE(II,JJ,SS,[M N Nzmax]) also specifies the dimensions of the
%   output. All arguments are optional and when they are omitted M =
%   MAX(II), N = MAX(JJ) and Nzmax = NNZ(S) are used as defaults.
%
%   S = FSPARSE(II,JJ,SS,SIZ,'nosort') produces a sparse matrix S
%   where the columns are not sorted with respect to row-indices (if
%   the input is sorted with respect to II, then the output is in
%   order anyway). This runs faster and the resulting unordered matrix
%   can be used in many, but not all, situations in
%   MATLAB. Cautionary: on some platforms an unordered sparse matrix
%   causes MATLAB to crash when performing certain operations.
%
%   One or both of the index-matrices II and JJ may be integers instead of
%   doubles. This is faster still since no typecast or deep copy is
%   needed. An 'integer' is the type which corresponds to the C declaration
%   'int' on the platform; -- on most modern platforms this is int32.
%
%   Differences between FSPARSE and SPARSE:
%     - FSPARSE generally allocates the sparse matrix so that it fits
%     exactly. However, zero elements resulting from cancellation or from
%     zero data causes a slight over-allocation. Thus, both FSPARSE([1 1],[1
%     1],[1 -1]) and FSPARSE(1,1,0) allocates space for one entry.
%
%     - The call S = FSPARSE(X), where X is sparse, always makes a deep copy
%     of the sparse matrix so that NNZ(S) = NZMAX(S).
%
%     - Integer indices are not supported by SPARSE.
%
%     - According to the specification above, the call FSPARSE(1:n,1,1) is
%     not allowed. The corrected version of this example is simply
%     FSPARSE([1:n]',1,1) where the row-index is associated with the first
%     dimension of the input matrix.
%
%     - Logical sparse matrices cannot be created by FSPARSE. Use logical
%     operations for this purpose.
%
%     - SPARSE allows the input (II,JJ,SS) itself to be sparse or even
%     logical sparse. This syntax is not supported by FSPARSE.
%
%   Examples:
%     N = 10;
%     S1 = fsparse(1:N,1:N,1);    % same as speye(N)
%     S2 = fsparse([1:N]',1:N,1); % same as sparse(ones(N));
%
%     % sparse matrix from triplet
%     ii = ceil(4*rand(1,N));
%     jj = ceil(4*rand(1,N));
%     ss = rand(1,N);
%     S3 = fsparse(ii,jj,ss); [nnz(S3) nzmax(S3)]
%     s3 = sparse(ii,jj,ss);  [nnz(s3) nzmax(s3)] % over-allocation
%
%     % 3-point stencil
%     S4 = fsparse([2:N-1]',[1:N-2; 2:N-1; 3:N]',[1 -2 1],[N N]);
%
%     % 5-point stencil on an N-by-N grid
%     N2 = N*N; ix = reshape(1:N2,N,N);
%     ix([1 N],:) = []; ix(:,[1 N]) = []; ix = ix(:);
%     S5 = fsparse(ix,[ix-N ix-1 ix ix+1 ix+N],[1 1 -4 1 1],[N2 N2]);
%
%     % band-matrix
%     B = rand(3,N); B(end,1) = 0; B(1,end) = 0;
%     % same as spdiags(B',[-1 0 1],N,N):
%     S6 = fsparse([[2:N 1]; [1:N]; [N 1:N-1]],1:N,B);
%
%     % circulant matrix
%     jj = 1+mod(cumsum([0:3; ones(N-1,4)]),N);
%     S7 = fsparse([1:N]',jj,1:4);
%
%   See also SPARSE, NNZ, NZMAX, SPDIAGS, FIND.


%   Hidden argument nthreads and OpenMP-support
%
%   The full call currently supported is actually
%     S = FSPARSE(II,JJ,SS,SIZ,[] | 'sort' | 'nosort',nthreads)
%   with nthreads an integer >= 1. The default value is
%   omp_get_max_threads(). The value of nthreads is remembered until
%   the mex-file is cleared from memory:
%     S1 = fsparse(...,4); % use 4 threads
%     S2 = fsparse(...);   % continues to use 4 threads
%     clear all;
%     S3 = fsparse(...);   % uses omp_get_max_threads()
%   The threaded version is a beta-release and is only supported for
%   a subset of the syntaxes covered by FSPARSE. Also, compilation
%   must be done under OpenMP, see MAKE.

%   Hidden output timing
%
%   With #define FSPARSE_TIME an extra output time vector is supported
%   as follows:
%     [S,t] = FSPARSE(...)
%   returns an 1-by-6 vector containing timings of the internal
%   phases of FSPARSE.
%     t(1) getix()-function
%     t(2) Part 1
%     t(3) Part 2
%     t(4) Part 3
%     t(5) Part 4
%     t(6) sparse_insert()
%   The suppport for timing is a beta-release and is only supported
%   for a subset of the syntaxes covered by FSPARSE. See MAKE.

%   Early performance observations for the serial version
%
%   In general, when the input indices and values are of the same
%   sizes, FSPARSE is about 30-70% faster than SPARSE depending on the
%   number of colliding indices, the usage of integer indices and the
%   'nosort' option. When converting from full to sparse storage
%   format, FSPARSE is about 50% faster than MATLAB. Finally, for the
%   'assembly' syntax, FSPARSE is about 50-60% faster than any
%   equivalent construction in MATLAB.

% S. Engblom 2013-12-03 (Revision, OpenMP added, hidden argument nthreads).
% S. Engblom 2010-01-13 (Minor revision, caution with 'nosort' added)
% S. Engblom 2004-11-12

error('.MEX-file not found on path.');
