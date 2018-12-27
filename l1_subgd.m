% Usage: opts.step_size_scheme is a function handle than returns the step size at t-th iteration.
% opts.maxIter: maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.

function [x,out]=l1_subgd(x0, A, b, mu, opts)
out.name="Subgradient Method";

n=length(x0);
a=opts.step_size_scheme;
cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
x=x0;

out.str=[];
for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter
		residual=A*x-b;
		g=zeros(n,1);
		g(x>0)=1;
		g(x<0)=-1;
		subgradient=A'*residual+cont_m(cont_id)*g;
		x_new=x-a(t)*subgradient;

		abserr=norm(x_new-x);
		x=x_new;

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