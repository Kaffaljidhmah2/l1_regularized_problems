% Usage: opts.step_size_scheme is a function handle that returns the step size at t-th iteration.
% opts.maxIter: maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.

function [x,out]=l1_pgd(x0, A, b, mu, opts)
out.name="Projected Gradient Method";

a=opts.step_size_scheme;
cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;


x1=max(x0,0);
x2=-min(x0,0);
out.str=[];
for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter
		residual=A*(x1-x2)-b;
		derivative=A'*residual;
		x1_new=x1-a(t)*(derivative+cont_m(cont_id));
		x2_new=x2-a(t)*(-derivative+cont_m(cont_id));
		x1_new=max(x1_new,0);
		x2_new=max(x2_new,0);

		abserr=norm([x1_new-x1;x2_new-x2]);
		x1=x1_new;
		x2=x2_new;
		if abserr<opts.tol
			outstr="absolute error drops below tol";
			break;
		end
		if abserr<opts.reltol*norm([x1,x2])
			outstr="relative error drops below reltol";
			break;
		end
	end
	if outstr==""
		outstr="maximum iterations exceeds";
	end
	out.str=[out.str;outstr];
end


x=x1-x2;

end