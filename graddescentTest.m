% computeDelLDelthetFull = @(Xlessr,prodMean,lambda,D,currThet) Xlessr'*(1 - sigmf(Xlessr*currThet,[2,0])) - prodMean' + lambda*D*currThet;


computeDelLDelthetFull = @(Xlessr,prodMean,lambda,D,currThet,N) 1/N*Xlessr'*tanh(Xlessr*currThet) - prodMean' + lambda*D*currThet;


computeLFull = @(Xlessr,meanProd,lambda,D,currThet) mean(log(cosh(Xlessr*currThet)) + log(2)) - meanProd*currThet + lambda/2*currThet'*D*currThet;
rmax = 12;
N = 2000;
r = 1;
nosperupdate = 200;
lambda = 0.1;

X = randn(N,rmax);
X = sign(X);

thet = randn(rmax - 1,1);
Xlessr = X(:,[1:r - 1,r + 1:rmax]);

Xprod = (X(:,r)*ones(1,rmax - 1)).*Xlessr;
meanProd = mean(Xprod,1);

Q = orth(randn(rmax - 1,rmax - 1));
eigenvals = abs(randn(rmax - 1,1) + 0.3);
D = Q*diag(eigenvals)*Q';

ss = 0.1;

iterloss = zeros(1,nosperupdate + 1);
iterloss(1) = computeLFull(Xlessr,meanProd,lambda,D,thet);
for ii = 1:nosperupdate
    thet = thet - ss*computeDelLDelthetFull(Xlessr,meanProd,lambda,D,thet,N);
    
    iterloss(ii + 1) = computeLFull(Xlessr,meanProd,lambda,D,thet);
end
plot(iterloss)