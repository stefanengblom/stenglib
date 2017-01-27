function [c,ia,ib] = fsetop(op,a,b)
%FSETOP Fast set operations based on hashing.
%   FSETOP(OP,A,B) performs the set operation OP on the columns of the
%   arrays A and B. The main differences from the set operations
%   implemented in MATLAB are:
%     - FSETOP runs faster because it is based on hashing instead of on
%     sorting.
%
%     - In FSETOP, the order of the elements in the input arrays are
%     retained in the output array, in contrast to the corresponding
%     MATLAB functions where the result is sorted. This means that all
%     indices for new elements are strictly increasing (except for IB
%     in the INTERSECT operation and for JA and JB in the UNION
%     operation). Also, in case of duplications, the first element
%     (lowest column number) is selected.
%
%     - FSETOP always uses the columns of its inputs as elements, like the
%     'rows' switch in the MATLAB commands, but for columns instead of
%     rows. In particular, column vectors are not special cases, they are
%     treated as single elements.
%
%     - FSETOP supports the following data types: double, single, char,
%     logicals, int8, uint8, int16, uint16, int32, uint32, int64, uint64
%     and cell-arrays containing any of the above data types.
%
%   Cell-arrays are treated as vectors and all output indices refer to a
%   linear ordering. The output C is always a row cell-vector and cells
%   containing the same data but of different types or shapes are
%   considered unequal.
%
%   C = FSETOP('check',A) computes hash-values for each column in A. The
%   output C is a row-vector containing uint32-type integers. This is not
%   a set-operation but is provided as a useful tool for computing
%   check-sums for large data-sets in a uniform way.
%
%   [C,IA,IB] = FSETOP('intersect',A,B) returns the columns common to
%   both A and B. C = A(:,IA) = B(:,IB).
%
%   [C,IA] = FSETOP('setdiff',A,B) returns the columns in A that are
%   not in B. C = A(:,IA).
%
%   [C,IA,IB] = FSETOP('setxor',A,B) returns the columns that are not
%   in the intersection of A and B. C = [A(:,IA) B(:,IB)].
%
%   [C,IA,IB,JA,JB] = FSETOP('union',A,B) returns the combined columns
%   from A and B but with no repetitions. C = [A(:,IA) B(:,IB)] and A
%   = C(:,JA), B = C(:,JB). Note that the outputs JA and JB are not
%   produced by the MATLAB-function UNION.
%
%   [B,IA,IB] = FSETOP('unique',A) returns the same columns as in A
%   but with no repetitions. B = A(:,IA) and A = B(:,IB).
%
%   [IA,IB] = FSETOP('ismember',A,B) returns a logical vector IA
%   containing 1 where the columns of A are also columns of B and 0
%   otherwise. IB contains the index in B of each column in A and zero if
%   no such index exists. A(:,IA) = B(:,IB(IA)).
%
%   Note: two elements are considered equal if and only if their
%   bitwise representations are identical. This can sometimes give
%   unexpected results. For example, with doubles, -0.0 ~= 0.0 and NaN
%   == NaN.
%
%   Examples:
%     % intersection of integers
%     a = int8(ceil(3*rand(2,10))), b = int8(ceil(3*rand(2,10)))
%     aandb = fsetop('intersect',a,b)
%
%     % unique strings in cell-array
%     strs = {'foobar' 'foo' 'foobar' 'bar' 'barfoo'}
%     strsunq = fsetop('unique',strs)
%
%     % find indices ix in b = a(ix)
%     a = randperm(6), b = ceil(6*rand(1,6))
%     [foo,ix] = fsetop('ismember',b,a)
%
%     % remove indices from set of indices
%     a = {[1 2 4] [2 3 5 6] [4 2 1] [2 3 5 6]' int32([1 2 4])};
%     b = {[4 2 1]};
%     c = fsetop('setdiff',a,b)
%
%     % union of structs
%     a = struct('foo',1,'bar',NaN,'foobar','hello')
%     b = struct('faa',2,'bar',Inf,'foobar','goodbye')
%     [cf,ia,ib] = fsetop('union',fieldnames(a),fieldnames(b));
%     ac = struct2cell(a); bc = struct2cell(b);
%     c = cell2struct([ac(ia); bc(ib)],cf)
%
%     % a single 32-bit checksum from strings
%     s = {'check','intersect','setdiff','setxor', ...
%          'union','unique','ismember'};
%     c = fsetop('check',fsetop('check',s)')
%
%   See also INTERSECT, SETDIFF, SETXOR, UNION, UNIQUE, ISMEMBER.

% S. Engblom 2010-02-09 (Revision)
% S. Engblom 2006-06-13 (Major revision)
% S. Engblom 2005-06-21
% Based on a concept by P-O Persson, COMSOL AB. The internal hash
% function used is based on a code by P. Hsieh, see
% http://www.azillionmonkeys.com/qed/hash.html.

error('.MEX-file not found on path.');
