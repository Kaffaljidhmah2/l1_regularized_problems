function [x,out] = l1__mosek(x0, A, b, mu, opts)
[m,n]=size(A);
Q=zeros(2*n+m,2*n+m);
Q(2*n+1:2*n+m,2*n+1:2*n+m)=eye(m);
Q=sparse(Q);
c=[zeros(n,1);mu*ones(n,1);zeros(m,1)];

L=[eye(n) -eye(n) zeros(n,m); -eye(n) -eye(n) zeros(n,m); A zeros(m,n) -eye(m)];
L=sparse(L);
blc=[-inf(2*n,1);b];
buc=[zeros(2*n,1);b];


res=mskqpopt(Q,c,L,blc,buc,[],[]);

x=res.sol.itr.xx(1:n);
out=res.sol.itr;
end