clear
load('lattice4by4.mat');

X = samples';
[nos,nov] = size(X);
thetT = zeros(nov,nov);
thetL = zeros(nov,nov);
for vnos = 1024
    thetT = zeros(nov,nov);
    for ii = 1:16
        Xii = X(1:vnos,ii)*ones(1,nov);
        Y = X(1:vnos,:);
        Y(:,ii) = ones(vnos,1);
        Y = Y.*Xii;
        init = randn(16,1);
        temp = logisticTraceLASSO(Y,init,0.05,2.^(1:-1.5:-20),32);
        thetT(:,ii) = temp;
        temp = logisticLASSO(Y,init,0.2,32*16);
        thetL(:,ii) = temp;
    end
    save(sprintf('lattice4by4TLresult_%d',vnos),'thetT','thetL','w');
end