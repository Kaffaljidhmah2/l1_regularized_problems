% Usage: opts.step_size_scheme is a function handle that returns the step size at t-th iteration.
% opts.maxIter: a list of maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.
% opts.gamma: list of smooth coefficient.

function [x,out]=l1_smooth_fgd(x0, A, b, mu, opts)
out.name="Fast Smooth Gradient";
out.str=[];

a=opts.step_size_scheme;
x_prev=x0;
x=x0;

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;

for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter(cont_id)
		y = x + ((t-2)/(t+1))*(x-x_prev);

		resid=A*y-b;
		nabla_H=y./opts.gamma(cont_id);
		mask=abs(nabla_H)>1;
		nabla_H(mask)=sign(nabla_H(mask));
		grad=A'*resid + cont_m(cont_id)*nabla_H;

		x_prev=x;
		x=y-a(t)*grad;

		prog=a(t)*norm(grad,2);
		if prog<opts.tol
			outstr="absolute progress drops below tol.";
			break;
		end
		if prog<opts.reltol*norm(x,2)
			outstr="relative progress drops below reltol.";
			break;
		end

	end
	if outstr==""
		outstr="maximum iterations exceeds";
	end
	out.str=[out.str;outstr];
end

end