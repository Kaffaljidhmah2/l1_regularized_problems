% Usage: opts.step_size is the step size.
% opts.maxIter: maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.delta_ is a small number.
% opts.gamma: list of smooth coefficient.

function [x,out]=l1_adagrad(x0, A, b, mu, opts)
out.name="Adagrad Method";

n=length(x0);
a=opts.step_size;
cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
x=x0;
r=zeros(n,1);
out.str=[];
for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter
		resid=A*x-b;
		nabla_H=x./opts.gamma(cont_id);
		mask=abs(nabla_H)>1;
		nabla_H(mask)=sign(nabla_H(mask));
		grad=A'*resid + cont_m(cont_id)*nabla_H;
		r=r+grad.*grad;
		v=a./(opts.delta_+sqrt(r));
		v=v.*grad;
		x=x-v;

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