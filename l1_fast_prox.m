% Usage: opts.step_size_scheme is a function handle that returns the step size at t-th iteration.
% opts.maxIter: a list of maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.

function [x,out]=l1_fast_prox(x0, A, b, mu, opts)
out.name="FISTA";
out.str=[];
n=length(x0);
a=opts.step_size_scheme;
x_prev=x0;
x=x0;

cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;

for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter(cont_id)
		y=x+((t-2)/(t+1))*(x-x_prev);
		x_prev=x;

		resid=A*y-b;
		y=y-a(t)*(A'*resid);
		x=zeros(n,1);
		thres=a(t)*cont_m(cont_id);
		mask=y>thres;
		x(mask)=y(mask)-thres;
		mask=y<-thres;
		x(mask)=y(mask)+thres;

		prog=norm(x-x_prev,2);

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