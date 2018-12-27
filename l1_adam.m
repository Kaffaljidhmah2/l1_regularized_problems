% Usage: opts.step_size_list is the step size list.
% opts.maxIter: maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.gamma: list of smooth coefficient.
% opts.rho_1 controls how strong the history affects the current first moment
% opts.rho_2 controls how strong the history affects the current second moment
% opts.delta_ is a small number.

function [x,out]=l1_adam(x0, A, b, mu, opts)
out.name="Adam Method";

n=length(x0);

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
x=x0;

out.str=[];
uni_t=1;
for cont_id=1:opts.cont_scheme
	a=opts.step_size_list(cont_id);
	outstr="";
	r=zeros(n,1);
	s=zeros(n,1);
	for t=1:opts.maxIter
		resid=A*x-b;
		nabla_H=x./opts.gamma(cont_id);
		mask=abs(nabla_H)>1;
		nabla_H(mask)=sign(nabla_H(mask));
		grad=A'*resid + cont_m(cont_id)*nabla_H;
		s=opts.rho_1.*s+(1-opts.rho_1).*grad;
		r=opts.rho_2.*r+(1-opts.rho_2).*(grad.*grad);
		s_=s./(1-opts.rho_1^uni_t);
		r_=r./(1-opts.rho_2^uni_t);
		v=a./(sqrt(r_)+opts.delta_);
		v=v.*s_;
		x=x-v;
		uni_t=uni_t+1;

		abserr=norm(v);

		if abserr<opts.tol
			outstr="absolute error drops below tol";
			break;
		end
		if abserr<opts.reltol*norm(x)
			outstr="relative error drops below reltol";
			break;
		end
	end
	if outstr==""
		outstr="maximum iterations exceeds";
	end
	out.str=[out.str ; outstr];
end

end