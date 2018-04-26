xtrue = randn(10,1);
A = randn(20,10);
y = A*xtrue + 0.01*randn(20,1);

gradx = @(x) 2*(A'*A*x - A'*y);
lin_srch = @(pk,x) (pk'*A'*y - pk'*A'*A*x)/(pk'*A'*A*pk);

x = cg(lin_srch,gradx,zeros(10,1),20);