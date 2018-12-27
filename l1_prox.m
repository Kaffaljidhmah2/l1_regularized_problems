% Usage: opts.step_size_scheme is a function handle that returns the step size at t-th iteration.
% opts.maxIter: a list of maximum iterations.
% opts.reltol: tolerance of relative error of x's of two succesive iterations.
% opts.tol: tolerance of absolute error of x's of two succesive iterations.
% opts.cont_scheme sets the continuation scheme.

function [x,out]=l1_prox(x0, A, b, mu, opts)
out.name="Proximal Gradient Method";
out.str=[];

n=length(x0);
a=opts.step_size_scheme;
cont_m=10.^(opts.cont_scheme-1:-1:0)*mu;
x=x0;

for cont_id=1:opts.cont_scheme
	outstr="";
	for t=1:opts.maxIter(cont_id)
		residual=A*x-b;
		y=x-a(t)*(A'*residual);

		thres=a(t)*cont_m(cont_id);
		x_new=zeros(n,1);
		mask=y>thres;
		x_new(mask)=y(mask)-thres;
		mask=y<-thres;
		x_new(mask)=y(mask)+thres;

		prog=norm(x_new-x);
		x=x_new;

		if prog<opts.tol
			outstr="absolute progress drops below tol.";
			break;
		end
		if prog<opts.reltol*norm(x)
			outstr="relative progress drops below reltol.";
			break;
		end

	end
	if outstr==""
		outstr="maximum iterations exceeds";
	end
	out.str=[out.str; outstr];
end

end