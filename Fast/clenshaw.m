function Y = clenshaw(A,B,C)
%CLENSHAW Evaluation of 3-term recurrences.
%   Y = CLENSHAW(A,B) evaluates the 3-term recurrence as defined by the
%   matrices A and B.
%
%   Let [M,N] be MAX(SIZE(A),SIZE(B)). Then the output Y is M-by-N and, for
%   1 <= j <= N, the recurrence is defined by Y(1,j) = A(1,j)+B(1,j), Y(2,j)
%   = A(2,j)*Y(1,j)+B(2,j) and, for i > 2, Y(i,j) =
%   A(i,j)Y(i-1,j)+B(i,j)Y(i-2,j).
%
%   The dimensions of the matrices A and B must either match or be
%   singletons. Singleton dimensions are as usual silently repeated.
%
%   Y = CLENSHAW(A,B,C) evaluates the 3-term recurrence defined by A and B,
%   multiplies it by the coefficients C and computes the sum of the result
%   by means of Clenshaws summation formula. C is M-by-K-by-N and the result
%   is K-by-N. The first dimension of C must be M while the second and third
%   dimensions of C can be singletons.
%
%   Examples:
%     % Fibonacci numbers
%     n = 10; clenshaw(1,[0; 0; ones(n-2,1)])'
%
%     % Legendre polynomials
%     x = linspace(-1,1); n = 4; j = [0.5 1:n]';
%     A = (2*j-1)./j*x;
%     B = (1-j)./j;
%     y = clenshaw(A,B);
%     figure, plot(x,y);
%
%     % two functions defined by Chebyshev coefficients
%     x = linspace(-1,1,50);
%     n = 4; j = [0; 1; 2*ones(n-1,1)];
%     A = j*x; B = [1; 0; -ones(n-1,1)];
%     C = [[1 0.75 0.5 0.25 0.125]' [1 -1.75 1.5 -1.25 2.125]'];
%     figure, plot(x,clenshaw(A,B,C));
%
%   See also FILTER.

% S. Engblom 2007-04-19 (Major revision)
% S. Engblom 2005-07-28

error('.MEX-file not found on path.');
