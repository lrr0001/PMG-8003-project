
clear
noc = 4;
submap = reshape(1:noc^2,noc,noc);
count = 1;
factors(2*(noc -2)^2 + 6*(noc - 2) + 4) = PGM_factor([],1);
h = zeros(noc,noc);
w = zeros(noc^2,noc^2);

for ii = 1:noc
    for jj = 1:noc
        temp = randn;
        w(submap(ii,jj),submap(ii,jj)) = temp;
        factors(count) = PGM_factor(submap(ii,jj),exp(temp).^([1,-1]));
        count = count + 1;
        if ii + 1 <= noc
            temp = randn;
            w(submap(ii,jj),submap(ii + 1,jj)) = temp;
            a11 = exp(temp);
            a12 = 1/a11;
            factors(count) = PGM_factor([submap(ii,jj),submap(ii + 1,jj)],[a11,a12;a12,a11]);
            count = count + 1;
        end
        if jj + 1 <= noc
            temp = randn;
            w(submap(ii,jj),submap(ii,jj + 1)) = temp;
            a11 = exp(temp);
            a12 = 1/a11;
            factors(count) = PGM_factor([submap(ii,jj),submap(ii,jj + 1)],[a11,a12;a12,a11]);
            count = count + 1;
        end
    end
end
w = w + w' - diag(diag(w));

%blks = cell(2,1);
%blk1set1 = submap(2:2:end,1:2:end);
%blk1set2 = submap(1:2:end,2:2:end);
%blk2set1 = submap(1:2:end,1:2:end);
%blk2set2 = submap(2:2:end,2:2:end);

%blks{1} = [blk1set1(:)',blk1set2(:)'];
%blks{2} = [blk2set1(:)',blk2set2(:)'];

gSampler = Gibbs_sampling(factors);
samples = gSampler.getGSamples(1024,128,16);
save('lattice4by4.mat','samples','factors','w');