fileID = fopen('mushroom.txt');
formatSpec = '%c,';

nof = 23;
nos = 8124;
sizeA = [nof,1];

%scanning in data
A = zeros(nof,nos);
for ii = 1:nos
A(:,ii) = fscanf(fileID,formatSpec,sizeA);
b = fscanf(fileID,'\n');
end
A = char(A);

keys = cell(nof,1);

% legend for data
keys{1} = 'ep';
keys{2} = 'bcxfks';
keys{3} = 'fgys';
keys{4} = 'nbcgrpuewy';
keys{5} = 'tf';
keys{6} = 'alcyfmnps';
keys{7} = 'adfn';
keys{8} = 'cwd';
keys{9} = 'bn';
keys{10} = 'knbhgropuewy';
keys{11} = 'et';
keys{12} = 'bcuezr?';
keys{13} = 'fyks';
keys{14} = 'fyks';
keys{15} = 'nbcgopewy';
keys{16} = 'nbcgopewy';
keys{17} = 'pu';
keys{18} = 'nowy';
keys{19} = 'not';
keys{20} = 'ceflpszn';
keys{21} = 'knbhrouwy';
keys{22} = 'acsynv';
keys{23} = 'glmpuwd';


%combining abundant, several, and numerous categories
A(22,A(22,:) == 'a'|A(22,:) == 'n'|A(22,:) == 'v') = 'a';


%converting characters into numbers
B = zeros(nof,nos);
catsize = zeros(nof,1);
for ii = 1:nof
    catsize(ii) = numel(keys{ii});
%    for jj = 1:numel(keys{ii})
%        B(ii,A(ii,:) == keys{ii}(jj)) = jj;
%    end
end
%reduce category size for feature with combined categories.
catsize(22) = catsize(22) - 2;
%keep = [1,2,3,5,7,8,9,11,13,14,17,19,20,22];
keep = [1,2,3,5];

nov = sum(catsize(keep) - 1);
binVars = zeros(nov + 1,nos);
count = 1;
varCorres = zeros(nov,1);
for ii = keep
    for jj = 1:catsize(ii) - 1
        trueVars = A(ii,:) == keys{ii}(jj);
        binVars(count,trueVars) = 1;
        binVars(count,~trueVars) = -1;
        varCorres(count) = ii;
        count = count + 1;
        
    end
end
binVars(count,:) = 1;

Y = ones(nov,1)*binVars(1,:);

nots = 4028;
thet = logisticTraceLASSO(Y(:,1:nots)'.*binVars(2:end,1:nots)',randn(nov,1),0.05,2.^(1:-3:-20),8);

% given the random initialization, future results may not be consistent.

%save('X_1_trace_lassoMuchroom.mat','thet');

%%
load('X_1_trace_lassoMuchroom.mat');

log_regress_accuracy = mean(sign(thet([2,5,7,10])'*binVars([3,6,8,11],nots + 1:end)) == sign(binVars(1,nots + 1:end)))

% logistic regression: 21.73%