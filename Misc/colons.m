function c = colons(a,b)
%COLONS Repeated colon operations.
%   C = COLONS(A,B) returns A:B with A and B row vectors of the same
%   length. This is just C = [A(1):B(1) A(2):B(2) ...].
%
%   No error-checking is performed.
%
%   Example:
%     a = 1;
%     b = 0:3;
%     c = colons(a,b) % returns c = [1 1 2 1 2 3]

% S. Engblom 2014-09-18

% a hack, really!
c = zeros(1,0);
for i = 1:numel(a)
  c = [c a(i):b(i)];
end
