clear; close all;
flag = 3;

tf = 10;
n = 1000;

u0 = zeros(n,1);
e = ones(n,1);

for i = 1:n
    u0(i) = i-1;
end

if flag == 1
    M = 5*speye(n, n);
elseif flag == 2
    M = spdiags([e 5*e e], -1:1, n, n);
    M(1,1) = 1; M(n,n) = 1;
elseif flag == 3
    M = spdiags([3*e 5*e e], -1:1, n, n);
    M(1,1) = 1; M(n,n) = 1;
elseif flag == 4
    M = spdiags([3*e 5*e e], -1:1, n, n);
    M(1,1) = 1; M(n,n) = 1;
    
    M = 0.1*M;
else
    M = eye( n,n );
end
rtol = 1e-4;

S = spdiags([e e e], -1:1, n, n);
% computing solutions
%options = [];
options = odeset('Mass',M,'RelTol', rtol ); % options for ode23s

sol2 = ode23t(@funMass,[0 tf],u0,options);

figure(2)
plot( sol2.x,sol2.y );


sol1 =  ode1s(@funMass,[0 tf],u0,{'mass',M,'jpattern',S,'rtol',rtol,'atol',rtol/100});
%sol1 =  ode1s(@funMass,[0 tf],u0,{'mass',M,'rtol',rtol,'atol',rtol/100});
%sol1 =  ode1s(@funMass,[0 tf],u0,{'rtol',rtol,'atol',rtol/100});
x = sol1.x;
y1 = sol1.y;

y = zeros(size(u0),size(x) );
y(:,1) = u0;



figure(1)
plot( sol1.x,sol1.y );



e(1)  = 0; % same initial condition
ey(1) = 0; % same initial condition

% for i = 2:size(sol1.y,2)
%     sol = ode23t(@funMass,[0 x(i)],u0,options);
%     y(:,i) = sol.y(:,size(sol.y,2));
%     
%     % comparison calculations
%     ey(i)   = y(:,i)'*y(:,i);
%     dy(:,i) = y(:,i) - y1(:,i);
%     e(i)    = dy(:,i)'*dy(:,i); % squared difference (each time)
% end
ey = sqrt( ey );

e  = sqrt( e  );
error = sqrt(e'*e)/norm(ey); % shall I take error/100? 

