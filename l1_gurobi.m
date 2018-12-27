function [x,out] = l1__gurobi(x0, A, b, mu, opts)
clear model;

[m,n]=size(A);
Q=zeros(2*n+m,2*n+m);
Q(2*n+1:2*n+m,2*n+1:2*n+m)=0.5*eye(m);
model.Q=sparse(Q);
model.obj=[zeros(n,1);mu*ones(n,1);zeros(m,1)];

model.A=sparse([eye(n) -eye(n) zeros(n,m); -eye(n) -eye(n) zeros(n,m); A zeros(m,n) -eye(m)]);
model.rhs=[zeros(2*n,1);b];
model.lb=repmat(-inf,2*n+m,1); % Important! By default the lower bound is zero! 
model.sense=[repmat('<',2*n,1);repmat('=',m,1)];

results=gurobi(model);
x=results.x(1:n);
out=results;

end