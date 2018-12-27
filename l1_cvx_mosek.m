function [x,out] = l1_cvx_mosek(x0, A, b, mu, opts)
n=length(x0);
cvx_solver mosek
cvx_begin
	variable x(n)
	minimize( 0.5*sum_square(A*x-b) + mu*norm(x,1))
cvx_end
out.optval=cvx_optval;
out.status=cvx_status;
end