function x = cg(line_search,grad,x,nosteps)
step = 0;
newgradx = grad(x);
graddiff = newgradx;
for ii = 1:nosteps
    denom = step'*graddiff;
    if denom ~= 0
        bet = newgradx'*graddiff/denom;
    else
        bet = 0;
    end
    pk = -(newgradx - step*bet);
    alph = line_search(pk,x);
    step = alph*pk;
    x = x + step;
    gradx = newgradx;
    newgradx = grad(x);
    graddiff = newgradx - gradx;
end

end
