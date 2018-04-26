function alph = line_search(AMean,d,A,pk,thet,inc,nos,lamb)

Apk = A*pk;
Athet = A*thet;

thetDsminsmu = lamb*thet'.*d'*pk - AMean*pk;
sDs = pk'.*d'*pk;

alph = 1;
for ii = 1:nos
f2 = mean(Apk.^2.*sech(alph*Apk + Athet).^2) + lamb*sDs;
f1 = mean(Apk.*tanh(alph*Apk + Athet)) + alph*lamb*sDs + thetDsminsmu;
if ii ~= nos
    alph = alph - inc*f1/f2;
else
    alph = alph - f1/f2;
end
end
%{
alph = 1;
for ii = 1:nos
    gradalph = mean(Apk.*tanh(Athet + alph*Apk)) + thetDsminsmu + alph*sDs;
    currSign = sign(gradalph);
    if abs(gradalph) < 1e-4
        return;
    end
    if ii ~= 1
        if currSign == prevSign
            return;
        end
    end
    prevSign = currSign;
    alph = alph - currSign*inc;
end
if alph > 1
    return;
else
    for ii = 1:nos
    gradalph = mean(Apk.*tanh(Athet + alph*Apk)) + thetDsminsmu + alph*sDs;
    currSign = sign(gradalph);
    
    if currSign == prevSign
        return;
    end
    prevSign = currSign;
    alph = alph*(1 - inc);
    end
end
%}