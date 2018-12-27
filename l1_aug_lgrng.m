% opts.maxIter: list of maximum iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.beta_ sets the ADMM step size.
% opts.tol : tolerance of absolute progress
% opts.reltol : tolerance of relative progress
% opts.subopt.maxIter : Inner loop
% opts.subopt.a : Inner step size

function [x,out]=l1_aug_lgrng(x0, A , b, mu, opts)
out.name="augmented lagrangian";
[m,n]=size(A);

%initialize
x=x0; u=A*x-b; z=-A'*u;
U=(eye(m)+opts.beta_.*(A*A'))^-1;

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
out.str=[];

for cont_id=1:opts.cont_scheme
	outstr="";
	mu_=cont_m(cont_id);
	for t=1:opts.maxIter(cont_id)
		
		% Solve Subproblem
		z_prev=z;
		residual=A*x-b;

		for sub_t=1:opts.subopt.maxIter
			ratio=(sub_t-2)/(sub_t+1);
			z_ = z + ratio*(z-z_prev);
			z_prev=z;
			u=U*(residual-opts.beta_*(A*z_));
			grad=opts.beta_*(A'*u+z_)-x;
			z=z_ - opts.subopt.a * grad;
			z(z>mu_)=mu_;z(z<-mu_)=-mu_;
		end
		u=U*(residual-opts.beta_*(A*z_));

		% ------

		xstep=opts.beta_.*(A'*u+z);
		x=x-xstep;
		%Stopping conditions
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

end