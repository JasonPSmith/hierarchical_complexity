n = 1000;
for i = 1:99
    for j = 1:100
    d(i,j) = 0.01*i;
    sigma(i,j) = 0.02*j;    
    [~,Rn(:,i,j)] = NormalisedHierarchicalComplexity(threshold_proportional(WeightedNetworkModel(n,'go3','log',sigma(i,j)),d(i,j))>0);
    end
end