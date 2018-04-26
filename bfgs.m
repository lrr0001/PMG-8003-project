function [x,Binv] = bfgs(line_search,grad,Binv,x,nosteps)
gradx = grad(x);
if norm(gradx) <1e-10
    return;
end
for ii = 1:nosteps
pk = -Binv*gradx;
alph = line_search(pk,x);
sk = alph*pk;
x = x + sk;
newgradx = grad(x);
if norm(newgradx) < 1e-10
    return;
end
yk = newgradx - gradx;
if norm(yk) < 1e-10
    continue;
end

skTyk = sk'*yk;
if skTyk == 0
    warning('Algorithm intends to divide by zero. Terminating now');
    return;
end

ykTBinv = yk'*Binv;
ykTBinvyk = ykTBinv*yk;
skskT = sk*sk';
skykTBinv = sk*ykTBinv;

Binv = Binv + (skTyk + ykTBinvyk)/(skTyk)^2*skskT - (skykTBinv + skykTBinv')/skTyk;
gradx = newgradx;
end

end
