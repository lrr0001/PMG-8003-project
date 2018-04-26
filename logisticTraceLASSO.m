function thet = logisticTraceLASSO(X,thet,lamb,mu,noiperupdate)
tic
noupdate = numel(mu);
[N,M] = size(X);
%Xlessr = X(:,[1:r - 1,r + 1:size(X,2)]);

%prodMat = (X(:,r).*ones(1,size(X,2) - 1)).*Xlessr;

%prodMean = mean(prodMat,1);
XMean = mean(X,1);
%normalizedProdMat = prodMat./sqrt(size(X,1));

computeDelLDelthet = @(d,currThet) 1/N*X'*tanh(X*currThet) - XMean' + lamb*d.*currThet;

%computeS = @(currThet,currMu) (X*diag(currThet.^2)*X' + currMu*eye(N))^(1/2);

%computeD = @(currS) diag(diag(X'*(currS\X)));
%Binv = eye(M);
for ii = 1:noupdate
    t = toc;
    if t > 300
        fprintf('ii = %d\n',ii);
        tic;
    end
        
%    S = computeS(thet,mu(ii));
    [U,sigmas,~] = svd(X*diag(thet.^2)*X');
    temp = X'*U;
    d = diag(temp*diag(1./sqrt(diag(sigmas) + mu(ii)))*temp');
    %D = computeD(S);
    %{
    for jj = noiperupdate
        normalizationComputation = X*thet;
        prevThet = thet;
        thet = thet - 0.1*computeDelLDelthet(normalizationComputation,D,thet);
        if any(isnan(thet))
            error('We have a NaN.');
        end
    end
    %}
    lin_srch = @(pk,thet) line_search(XMean,d,X,pk,thet,1,1,lamb);
    gradthet = @(thet) computeDelLDelthet(d,thet);
    %[thet,Binv] = bfgs(lin_srch,gradthet,Binv,thet,noiperupdate);
    thet = cg(lin_srch,gradthet,thet,noiperupdate);
end