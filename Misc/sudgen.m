function G = sudgen(S,ix)
%SUDGEN Generate Sudoku puzzle.
%   G = SUDGEN(S) generates a uniquely solvable Sudoku G from a
%   fully solved Sudoku S.
%
%   G = SUDGEN(S,IX) does the same thing, but also ensures that the
%   index of the Sudoku is <= IX. The default is IX = 1.
%
%   SUDGEN depends on the random number generator via RANDPERM.
%
%   Example:
%     S = [9     2     5     6     3     1     8     4     7; ...
%          6     1     8     5     7     4     2     9     3; ...
%          3     7     4     9     8     2     5     6     1; ...
%          7     4     9     8     2     6     1     3     5; ...
%          8     5     2     4     1     3     9     7     6; ...
%          1     6     3     7     9     5     4     8     2; ...
%          2     8     7     3     5     9     6     1     4; ...
%          4     9     1     2     6     7     3     5     8; ...
%          5     3     6     1     4     8     7     2     9];
%   G1 = sudgen(S,1)
%   G2 = sudgen(S,2) % takes more time!
%
%   S-sudsolve(G1), S-sudsolve(G2)
%
%   See also SUDSOLVE.

% S. Engblom 2016-10-26

% default
if nargin < 2, ix = 1; end

% *** check that S is correct

% attempt removal of elements in this random order
p = randperm(81);

G = S;
for rem = p
  % remove element, then attempt to solve
  G_ = G;
  G_(rem) = 0;
  try
    [~,jx] = sudsolve(G_); % should throw an error when ix > jx!
    if jx <= ix
      G = G_; % accept removal of element
    end
  catch
    ; % otherwise attempt removing the next
  end
end
