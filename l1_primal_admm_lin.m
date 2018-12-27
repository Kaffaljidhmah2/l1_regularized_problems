% opts.maxIter: list of maximum iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.beta_ sets the ADMM step size.
% opts.gamma_ sets the linearization step size.
% opts.gaptol : tolerance of duality gap
% opts.tol : tolerance of absolute progress
% opts.reltol : tolerance of relative progress

function [x,out]=l1_primal_admm_lin(x0, A , b, mu, opts)
out.name="Primal ADMM lin";
[m,n]=size(A);

%initialize
x=x0;y=x0;u=A'*(A*x-b);

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
out.str=[];

for cont_id=1:opts.cont_scheme
	outstr="";
	mu_beta=cont_m(cont_id)/opts.beta_; % save computation

	for t=1:opts.maxIter(cont_id)
		residual=A*x-b;
		grad=A'*residual + opts.beta_*(x-y+u);
		grad= opts.gamma_ * grad;

		x=x- grad;

		y_tmp=x+u;
		y=zeros(n,1);
		mask=y_tmp>mu_beta;
		y(mask)=y_tmp(mask)-mu_beta;
		mask=y_tmp<-mu_beta;
		y(mask)=y_tmp(mask)+mu_beta;

		u=u+x-y;

		%Stopping conditions
		prog=norm(grad);
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