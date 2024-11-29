function [C,Ck] = StatisticalComplexity_Entropic(A)

n = size(A,1);
d = sum(sum(A))/n/(n-1);
K = sum(A);

pij = (A./repmat(K,n,1)).*A;

Si = nansum(pij.*log(pij));

S = sum(Si)/n/log(n-1);

for k = 1:10
    ER = threshold_proportional(squareform(rand(1,n*(n-1)/2)),d)>0;
    
    Ke = sum(ER);

    peij = (ER./repmat(Ke,n,1)).*ER;

    Sei = nansum(peij.*log(peij));   
    Se(k) = sum(Sei)/n/log(n-1);

    Qi = (nansum(0.5*(pij+peij).*log(0.5*(pij+peij)))-0.5*(Si+Sei))/log(2);

    Q(k) = nanmean(Qi);

    Ck(k) = Q(k)*S;
end

C = nanmean(Ck);


