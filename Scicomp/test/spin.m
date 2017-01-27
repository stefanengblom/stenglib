function spin(ix)
%SPIN Tests for SCICOMP.

% S. Engblom 2013-07-15 (GAUSSQD added)
% S. Engblom 2008-06-17

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4 @l_spin5 @l_spin6 ...
          @l_spin7 @l_spin8 @l_spin9};

n = size(ftests,2);
stests = cellstr(reshape(sprintf('Test #%d (qaussqd)',1:n-1)',[],n-1)')';
stests = [stests {'Test #9 (rtsafe)'}];
if nargin == 0, ix = 1:n; end
runtest('SPIN (scicomp)',ftests(ix),stests(ix));

%--------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Test of GAUSSPD/'charlier'.

ok = 1;

% scalar order of polynomial
x = linspace(0,5,10);
[A,B] = gausspd('charlier',x,0,1);
y = clenshaw(A,B);
ok = ok && norm(y-1) <= 1e-14;


[A,B] = gausspd('charlier',x,1,4);
y = clenshaw(A,B,[0; 1]);
ok = ok && norm(2*y-(4-x)) <= 1e-14;

a = 9;
[A,B] = gausspd('charlier',x,2,a);
y = clenshaw(A,B,[0; 0; 1]);
ok = ok && norm(y- ...
                polyval([1 -1-2*a a*a]/sqrt(2)/a,x),inf) <= 1e-14;

% multiple orders
a = 17;
[A,B] = gausspd('charlier',x,2,a);
y = clenshaw(A,B);
ok = ok && norm(y(3,:)-...
                polyval([1 -1-2*a a*a]/sqrt(2)/a,x),inf) <= 1e-14;
ok = ok && norm(y(2,:)- ...
                polyval([-1 a]/sqrt(a),x),inf) <= 1e-14;

% recurrence
a = 13;
[A,B] = gausspd('charlier',x,7,a);
y = clenshaw(A,B);
y3 = ((6+a)-x)/sqrt(7*a).*y(7,:)-sqrt(6/7)*y(6,:);
ok = ok && norm(y3-y(8,:),inf) <= 1e-12;

% difference equation
[A,B] = gausspd('charlier',[x-1 x x+1],6,a);
y = reshape(clenshaw(A,B,[0 0 0 0 0 0 1]'),[],3);
y3 = (x'+a-6)/a.*y(:,2)-x'/a.*y(:,1);
ok = ok && norm(y3-y(:,3),inf) <= 1e-12;

% value at zero
a = 5.5;
[A,B] = gausspd('charlier',0,5,a);
y = clenshaw(A,B)';
ok = ok && norm(y-a.^([0:5]/2)./sqrt(gamma(1:6)),inf) <= 1e-14;

% derivative with respect to parameter
n = 5;
[A,B] = gausspd('charlier',x,n,a);
yda = clenshaw(A,B);
yda = [sqrt(n/a) -n/(2*a)]*yda([n n+1],:);
[A,B] = gausspd('charlier',x,n,a+1e-20i);
ydan = imag(clenshaw(A,B,[zeros(n,1); 1]))*1e20;
ok = ok && norm(yda-ydan,inf) <= 1e-14;

% derivative with respect to variable
n = 7;
nn = 1:n;

coeff = cumsum(log(nn));
nn = -nn(end:-1:1);
coeff = (coeff(end)-[0 coeff(1:end-1)])+log(a)*nn;
coeff = exp(0.5*coeff)./nn;

[A,B] = gausspd('charlier',x,n,a);
ydx = clenshaw(A,B,[coeff'; 0]);

[A,B] = gausspd('charlier',x+1e-20i,n,a);
ydxn = imag(clenshaw(A,B,[zeros(n,1); 1]))*1e20;

ok = ok && norm(ydx-ydxn,inf) <= 1e-14;

%--------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Test of GAUSSQD/'charlier'.

ok = 1;

% normalization
[x,w] = gaussqd('charlier',5,10);
[A,B] = gausspd('charlier',x,4,10);
chp = clenshaw(A,B);
ok = ok && norm((chp.*chp)*w-1,inf) <= 1e-14;
ok = ok && norm((chp.*chp([5 1:4],:)*w),inf) <= 1e-14;

% third output
[x,w,y] = gaussqd('charlier',6,3);
[A,B] = gausspd('charlier',x,6,3);
y2 = clenshaw(A,B);
ok = ok && norm(y-y2(1:end-1,:),inf) <= 1e-12;

% stress
for alpha = linspace(0.1,1,20)
  [x,w] = gaussqd('charlier',5,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(1,5,20)
  [x,w] = gaussqd('charlier',10,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(5,20,30)
  [x,w] = gaussqd('charlier',20,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(20,1000,50)
  [x,w] = gaussqd('charlier',40,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end

for alpha = linspace(0.1,1,50)
  [x,w] = gaussqd('charlier',7,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(1,5,70)
  [x,w] = gaussqd('charlier',13,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(5,20,95)
  [x,w] = gaussqd('charlier',24,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end
for alpha = linspace(20,1000,153)
  [x,w] = gaussqd('charlier',40,alpha);
  ok = ok && norm(sum(w)-1,inf) <= 1e-10;
end

%--------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 Test of GAUSSPD/'meixner'.

ok = 1;

%--------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 Test of GAUSSQD/'meixner'.

ok = 1;

%--------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 Test of GAUSSPD/'hahn'.

ok = 1;

%--------------------------------------------------------------------------
function ok = l_spin6
%L_SPIN6 Test of GAUSSQD/'hahn'.

ok = 1;

% warm-ups
x = gaussqd('hahn',2,2,3,15);
xx = [31/14-1/14*709^(1/2); 31/14+1/14*709^(1/2)];
ok = ok && norm(x-xx,inf) <= 1e-14;

x = gaussqd('hahn',5,2,3,15);
xx = [0.1370075211, 2.454204755, 8.583365284, 29.15895577, 179.6664667]';
ok = ok && norm((x-xx)./xx,inf) <= 1e-9;

% special cases for measures of finite support
x = gaussqd('hahn',4,-3,-4.5,1.3);
ok = ok && norm(x-[0:3]',inf) <= 1e-14;

x = gaussqd('hahn',4,-3,-3.01,1.1);
ok = ok && norm(x-[0:3]',inf) <= 2e-14;

%--------------------------------------------------------------------------
function ok = l_spin7
%L_SPIN7 Test of Krawtchouk and Chebyshev polynomials.

ok = 1;

%--------------------------------------------------------------------------
function ok = l_spin8
%L_SPIN8 Test of limiting cases of parameters.

ok = 1;

% This is really difficult to test! Here, no warnings should be generated.
[x,w] = gaussqd('meixner',5,-1,-4);
ok = ok && norm(x-[0:4]',inf) <= 1e-12;

[x,w] = gaussqd('meixner',4,-1,-9);

[x,w] = gaussqd('meixner',10,-1,-9);
ok = ok && norm(x-[0:9]',inf) <= 1e-12;

[x,w] = gaussqd('krawtchouk',5,3/4,4);
ok = ok && norm(x-[0:4]',inf) <= 1e-12;

[x,w] = gaussqd('krawtchouk',7,2/3,6);
ok = ok && norm(x-[0:6]',inf) <= 1e-12;

[x,w] = gaussqd('hahn',3,-2,3,-2);
ok = ok && norm(x-[0:2]',inf) <= 1e-12;

[x,w] = gaussqd('hahn',4,-3,-3,2);
ok = ok && norm(x-[0:3]',inf) <= 1e-12;

[x,w] = gaussqd('hahn',2,2,2,8);
ok = ok && norm(x-[1 18]',inf) <= 1e-12;

%--------------------------------------------------------------------------
function [f,df] = l_tfun(x)
%L_TFUN Small testfunction for L_SPIN1.

f = (0.5 <= x & x <= 1.5).*(1-x)+(x < 0.5)*0.5-(1.5 < x)*0.5;
if nargout > 1, df = (0.5 <= x & x <= 1.5).*(-1); end

%--------------------------------------------------------------------------
function ok = l_spin9
%L_SPIN9 Test of RTSAFE.

ok = 1;

x1 = [linspace(0,1,10) 1 linspace(1,2,10)]';
x2 = [linspace(1,2,10) 1 linspace(0,1,10)]';

x = rtsafe(x1,l_tfun(x1),x2,l_tfun(x2),@l_tfun);
ok = ok & all(x == 1);

%--------------------------------------------------------------------------
