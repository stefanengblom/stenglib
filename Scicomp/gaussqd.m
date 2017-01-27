function [x,w,y] = gaussqd(name,n,p1,p2,p3)
%GAUSSQD Quadratures for discrete measures.
%   [X,W] = GAUSSQD(NAME,N,P1,P2,...) returns abscissas X and weights W
%   for the Nth discrete quadrature associated with the polynomial NAME
%   using the parameters P1, P2, ...
%
%   The resulting summation formula is
%     sum_0^inf f(x) weight(x) = sum_i f(x(i))*w(i),
%   where the weight-function is defined in GAUSSPD.
%
%   [X,W,Y] = GAUSSQD(NAME,N,P1,P2,...) additionally returns an N-by-N
%   matrix Y containing values of the normalized polynomials of order < N
%   in all the quadrature points X. Each column Y(:,i) contains the values
%   of the polynomials evaluated at X(i).
%
%   Unlike GAUSSPD, GAUSSQD carefully checks the parameters P1, P2, ... to
%   see if the corresponding set of polynomials makes sense. The
%   conditions put on the coefficients are that enough moments exist and
%   that the measure is strictly positive. In addition, the order of the
%   polynomial must be bounded whenever the natural domain of summation is
%   limited.
%
%   Reference:
%     [1] S. Engblom: "Gaussian quadratures with respect to discrete
%     measures". Technical Report 2006-007, Dept of Information
%     Technology, Uppsala University, 2006. Available at
%     http://www.it.uu.se/research.
%
%   Example:
%     f = inline('(x.^5+2*x.^4+3*x.^3+4*x.^2+5*x)');
%     [x,w] = gaussqd('charlier',2,1); sum1 = f(x)'*w
%     [x,w] = gaussqd('charlier',3,1); sum2 = f(x)'*w
%     [x,w] = gaussqd('charlier',4,1); sum3 = f(x)'*w
%
%   See also GAUSSPD.

% S. Engblom 2006-01-12

if ~isscalar(n) || n ~= ceil(n) || n <= 0
  error('Order of quadrature must be a positive integer.');
end

wmsg = 'Parameter(s) out of range. Quadrature may not make sense.';
switch name
 case 'charlier'
  if nargin ~= 3 || ~isscalar(p1)
    error('Parameter of Charlier polynomial must be a single scalar.');
  end
  ok = p1 > 0;
  if ~ok, warning(wmsg); end

  % construct a symmetric positive definite and tridiagonal matrix J
  nn = 1:n;
  b = -sqrt(nn*p1);
  b(end) = 0;
  J = fsparse([[1 1:n-1]; nn; [2:n n]],nn, ...
              [b([end 1:end-1]); nn+p1-1; b],[],'nosort');

  % eigenvalues and -vectors
  [y,x] = eig(full(J));

  % the abscissas are the eigenvalues of J
  x = diag(x);

  % the weights can be found from the first element of each eigenvector
  w = reshape(y(1,:).^2,[],1);

  % the values of the polynomials themselves can be found from the
  % eigenvectors (normalizing the first element)
  if nargout > 2
    y = tprod(y,1./y(1,:),[1 2],[3 2]);
  end

 case 'krawtchouk'
  if nargin ~= 4 || ~isscalar(p1) || ~isscalar(p2)
    error('Parameters of Krawtchouk polynomial must be two scalars.');
  end
  [x,w,y] = gaussqd('meixner',n,-p1/(1-p1),-p2);

 case 'meixner'
  if nargin ~= 4 || ~isscalar(p1) || ~isscalar(p2)
    error('Parameters of Meixner polynomial must be two scalars.');
  end

  % Conditions on (1) convergence of the inner product, (2) the measure
  % being positive, (3) when p2 is a non-positive integer (i.e. the
  % measure has finite support), then n must be sufficiently small.
  xmax = inf;
  if p2 <= 0 && p2 == ceil(p2)
    xmax = -p2;
  end
  pp = ceil([-p2 1-p2 -1-p2]);
  pp = min(max(0,pp),xmax);
  pp = sign(p1.*(p2+pp-1)).^pp;
  ok = xmax < inf || abs(p1) < 1 || abs(p1) == 1 && -p2 >= 2*n;
  ok = ok && all(pp > 0);
  ok = ok && n <= xmax+1;
  if ~ok, warning(wmsg); end

  nn = 1:n;
  b = -sqrt(p1*nn.*(p2+nn-1));
  b(end) = 0;
  J = fsparse([[1 1:n-1]; nn; [2:n n]],nn, ...
              [b([end 1:end-1]); (p1+1)*(nn-1)+p1*p2; b],[],'nosort');

  [y,x] = eig(full(J./(1-p1)));

  x = diag(x);
  w = reshape(y(1,:).^2,[],1);

  if nargout > 2
    y = tprod(y,1./y(1,:),[1 2],[3 2]);
  end

 case 'chebyshev'
  if nargin ~= 3 || ~isscalar(p1)
    error('Parameter of Chebyshev polynomial must be a single scalar.');
  end
  [x,w,y] = gaussqd('hahn',n,1-p1,1,1-p1);

 case 'hahn'
  if nargin ~= 5 || ~isscalar(p1) || ~isscalar(p2) || ~isscalar(p3)
    error('Parameters of Hahn polynomial must be three scalars.');
  end
  % Conditions on (1) convergence of the inner product, (2) the measure
  % being positive, (3) the denominator not containing zeros and (4) when
  % p1 and/or p2 are non-positive integers (i.e. the measure has finite
  % support), then n must be sufficiently small.
  xmax = inf;
  if p1 <= 0 && p1 == ceil(p1)
    xmax = -p1;
  end
  if p2 <= 0 && p2 == ceil(p2)
    xmax = min(xmax,-p2);
  end
  pp = [p1 p2 p3];
  pp = ceil([-pp 1-pp -1-pp]);
  pp = min(max(0,pp),xmax);
  pp = sign((p1+pp-1).*(p2+pp-1).*(p3+pp-1)).^pp;
  ok =  xmax < inf || p3-p1-p2 >= 2*n;
  ok = ok && all(pp > 0);
  ok = ok && (p3 > 0 || p3 == ceil(p3) && xmax <= -p3);
  ok = ok && n <= xmax+1;
  if ~ok, warning(wmsg); end

  nn = 1:n;
  w = p3-p1-p2;
  w1 = w-2*nn+1;
  if w == -1
    % Special case for essentially the Gauss-Chebyshev formula. It's a guess
    % that this is the only singular case.
    a = [p1*p2/(w-1) ...
        polyval([1-p1-p2-p3 -(1-p1-p2-p3)*w (1+w)*p1*p2],nn(2:end)-1)./ ...
        (w1(2:end).*(w1(2:end)+2))];
  else
    a = polyval([1-p1-p2-p3 -(1-p1-p2-p3)*w (1+w)*p1*p2],nn-1)./(w1.*(w1+2));
  end
  n2 = 1:n-1;
  w2 = w-2*n2;
  b = [-sqrt(n2.*(w+1-n2).* ...
             (p1+n2-1)./(w2+1).* ...
             (p2+n2-1)./(w2+1).* ...
             (p1+w-n2)./(w2+2).* ...
             (p2+w-n2)./w2) 0];
  J = fsparse([[1 1:n-1]; nn; [2:n n]],nn, ...
              [b([end 1:end-1]); a; b],[],'nosort');

  [y,x] = eig(full(J));

  x = diag(x);
  w = reshape(y(1,:).^2,[],1);

  if nargout > 2
    y = tprod(y,1./y(1,:),[1 2],[3 2]);
  end
 otherwise
  error('Unknown type of polynomial.');
end
