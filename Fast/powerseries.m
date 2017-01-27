y = powerseries(c,x,tol)
%POWERSERIES Sum power series.
%   Y = POWERSERIES(C,X,tol) evaluates the power series with coefficients C
%   at points X where both C and X are vectors. The output Y is the same
%   size as X and contains Y(i) = sum_(j >= 1) C(j)*X(i)^(j-1). The sum is
%   truncated when the last contribution is relatively less than tol.
%
%   Note: use
%     warning('off','powerseries:w1');
%   to turn the warning for inaccurate evaluation off. As an
%   alternative, you may explicitly set C(end) = 0.
%
%   Example:
%     % the complex exponential function
%     c = 1./cumprod([1 1:10]);
%     x = complex(linspace(-2,3),linspace(-4,6));
%     y = powerseries(c,x,0.5e-7);
%     figure, plot(real(x),real(y),'b',real(x),real(exp(x)),'r');
%     legend('Series','Exp()');

% S. Engblom 2007-08-17 (Minor revision)
% S. Engblom 2007-06-13 (Minor revision)
% S. Engblom 2007-01-22

error('.MEX-file not found on path.');
