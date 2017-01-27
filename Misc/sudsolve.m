function [S,ix] = sudsolve(G)
%SUDSOLVE Solve Sudoku puzzle.
%   S = SUDSOLVE(G) solves the Sudoku puzzle with initial conditions
%   given by G. The algorithm used is a recursive greedy selection
%   procedure.
%
%   [S,IX] = SUDSOLVE(G) does the same thing but also returns the
%   index IX of the Sudoku. The index is defined as the highest branch
%   point taken by the algorithm and indicates the degree of
%   difficulty. For most human solvable Sudokus, IX <= 2.
%
%   Example:
%     % a Sudoku of medium difficulty:
%     G = zeros(9,9);
%     G([1 2 6 10 13 14 16 18 24 25 29 36 37 41 ...
%        45 46 53 57 58 64 66 68 69 72 76 80 81]) = ...
%     [7 5 8 3 4 9 2 6 1 4 6 8 5 1 3 1 2 1 6 8 6 2 9 7 7 6 1];
%     % (Published in Sydsvenska Dagbladet, Malmoe, Sweden, 2005-07-13)
%     S = sudsolve(G)
%
%     % a rather tricky one:
%     M = zeros(9,9);
%     M([2 5 8 10 18 21 25 31 33 37 41 ...
%        45 49 51 57 61 64 72 74 77 80]) = ...
%     [6 8 4 2 3 4 7 8 7 3 1 4 6 5 5 6 4 2 3 6 8];
%     [S,ix] = sudsolve(M) % takes some time
%     % (From Cleve's corner) 

% S. Engblom 2010-02-25 (Added Moler's example)
% S. Engblom 2005-07-15

% used for tracking branches
stack = cell(3,0);
kmax = 0;

% start with all possibilities
S = true(9,9,9);

% insert given values G
[i,j,v] = find(G);

% check G
if any(v ~= ceil(v)) || any(v < 1) || any(9 < v)
  error('Illegal initial conditions.');
end
if ndims(G) > 2 || any(size(G) ~= 9)
  error('Illegal initial conditions.');
end

for k = 1:size(v,1)
  ii = i(k); jj = j(k); vv = v(k);

  if ~S(ii,jj,vv), error('Inconsistent initial conditions.'); end

  % set individual number
  S(ii,jj,:) = false;

  % remove possibilities in same column...
  S(:,jj,vv) = false;

  % ...same row...
  S(ii,:,vv) = false;

  % ...and same 3-by-3 square
  I = ceil(ii/3); I = 3*(I-1)+1:3*I;
  J = ceil(jj/3); J = 3*(J-1)+1:3*J;
  S(I,J,vv) = false;
end

% remove possibilities in a greedy fashion
while true
  % no more possibilities?
  while ~any(S(:))
    % done/unsolvable/pop stack?
    if all(G(:))
      % the answer is in G...
      S = G;
      if nargout > 1, ix = kmax; end
      return;
    elseif size(stack,2) == 0
      error('Sudoku is unsolvable.');
    else
      S = stack{1,end};
      G = stack{2,end};
      kmax = stack{3,end};
      stack = stack(:,1:end-1);
    end
  end

  % find the branch with the lowest complexity
  for k = 1:9
    [ii,jj] = find(sum(S,3) == k,1);
    if ~isempty(ii), break; end
  end % will always succeed
  kmax = max(k,kmax);

  % increase stack if branch is taken
  if k > 1
    v = find(S(ii,jj,:));
    for l = 1:k
      vv = v(l);

      % S and G must be saved
      SS = S; GG = G;

      % make the move
      GG(ii,jj) = vv;
      SS(ii,jj,:) = false;
      SS(:,jj,vv) = false;
      SS(ii,:,vv) = false;
      I = ceil(ii/3); I = 3*(I-1)+1:3*I;
      J = ceil(jj/3); J = 3*(J-1)+1:3*J;
      SS(I,J,vv) = false;

      % increase stack
      stack(:,end+1) = {SS; GG; kmax};
    end

    % continue with the front of the stack
    S = stack{1,end};
    G = stack{2,end};
    kmax = stack{3,end};
    stack = stack(:,1:end-1);
  else
    % no need for stack
    vv = find(S(ii,jj,:));

    % make the move
    G(ii,jj) = vv;
    S(ii,jj,:) = false;
    S(:,jj,vv) = false;
    S(ii,:,vv) = false;
    I = ceil(ii/3); I = 3*(I-1)+1:3*I;
    J = ceil(jj/3); J = 3*(J-1)+1:3*J;
    S(I,J,vv) = false;
  end
end
