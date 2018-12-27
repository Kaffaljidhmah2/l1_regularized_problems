% opts.maxIter: list of maximum iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.beta_ sets the step size.
% opts.gaptol : tolerance of duality gap
% opts.tol : tolerance of absolute progress
% opts.reltol : tolerance of relative progress

function [x,out]=l1_dual_admm(x0, A , b, mu, opts)
out.name="Dual ADMM";
[m,n]=size(A);

%initialize
x=x0; u=A*x-b; z=-A'*u;

U=(eye(m)+opts.beta_.*(A*A'))^-1;

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
out.str=[];

for cont_id=1:opts.cont_scheme
	outstr="";
	tmp_mu=cont_m(cont_id);
	for t=1:opts.maxIter(cont_id)
		u=U*(opts.beta_.* (A*(x-z)) - b);
		A_u=A'*u;
		z=x-A_u;
		z(z>tmp_mu)=tmp_mu;
		z(z<-tmp_mu)=-tmp_mu;
		xstep=A_u+z;
		x=x-xstep;

		%Stopping conditions
		if norm(A'*u,'inf')<=tmp_mu
			fprintf("ball constraint met\n");
			primal=0.5*norm(A*opts.beta_*x-b,2)^2+tmp_mu*norm(x*opts.beta_,1);
			dual=-0.5*norm(u,2)^2-b'*u;
			if primal<dual - opts.gaptol
				fprintf("BUG !\n");
			end
			if primal-dual < opts.gaptol
				outstr="duality gap drops below gaptol";
				break
			end
		end
		prog=norm(xstep);
		if prog<opts.tol
			outstr="absolute progress drops below tol";
			break
		end
		if prog<opts.reltol*norm(x)
			outstr="relative progress drops below reltol";
			break
		end
	end
	if outstr==""
		outstr="maximum iterations exceeds";
	end
	out.str=[out.str ; outstr];
end

x=x*opts.beta_;
end