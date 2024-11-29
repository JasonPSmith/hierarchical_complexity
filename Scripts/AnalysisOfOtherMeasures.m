%% Entropic
for i = 1:100
n(i) = 50 + randi(2450);
d(i) = rand();
Aer = threshold_proportional(squareform(rand(1,n(i)*(n(i)-1)/2)),d(i))>0;
Arg = threshold_proportional(RandomGeometricGraph(n(i),3),d(i))>0;
Ahg = threshold_proportional(WeightedNetworkModel(n(i),'go3','log',0.2),d(i))>0;
Arh = randmio_und(Ahg,2);
Cer(i) = StatisticalComplexity_Entropic(Aer);
Crg(i) = StatisticalComplexity_Entropic(Arg);
Chg(i) = StatisticalComplexity_Entropic(Ahg);
Crh(i) = StatisticalComplexity_Entropic(Arh);
end

%% Structural 
n = 100;
for i = 1:100
    d(i) = rand();

    Aer = threshold_proportional(squareform(rand(1,n*(n-1)/2)),d(i))>0;
    Arg = threshold_proportional(RandomGeometricGraph(n,3),d(i))>0;
    Ahg = threshold_proportional(WeightedNetworkModel(n,'go3','log',0.2),d(i))>0;
    Arh = randmio_und(Ahg,2);
    
    dGer_d(i) = StatisticalComplexity_Structural(Aer); 
    dGrg_d(i) = StatisticalComplexity_Structural(Arg);
    dGhg_d(i) = StatisticalComplexity_Structural(Ahg);
    dGrh_d(i) = StatisticalComplexity_Structural(Arh);
end

%% Structural 
n = 100;
for i = 1:100
    d(i) = rand();

    Aer = threshold_proportional(squareform(rand(1,n*(n-1)/2)),d(i))>0;
    Arg = threshold_proportional(RandomGeometricGraph(n,3),d(i))>0;
    Ahg = threshold_proportional(WeightedNetworkModel(n,'go3','log',0.2),d(i))>0;
    Arh = randmio_und(Ahg,2);
    
    dGer_d(i) = StatisticalComplexity_Structural(Aer); 
    dGrg_d(i) = StatisticalComplexity_Structural(Arg);
    dGhg_d(i) = StatisticalComplexity_Structural(Ahg);
    dGrh_d(i) = StatisticalComplexity_Structural(Arh);
end

nn = 100:100:1000;
dn = 0.01;
for i = 1:length(nn)    
    Aer = threshold_proportional(squareform(rand(1,nn(i)*(nn(i)-1)/2)),dn)>0;
    Arg = threshold_proportional(RandomGeometricGraph(nn(i),3),dn)>0;
    Ahg = threshold_proportional(WeightedNetworkModel(nn(i),'go3','log',0.2),dn)>0;
    Arh = randmio_und(Ahg,2);
    
    dGer_n(i) = StatisticalComplexity_Structural(Aer); 
    dGrg_n(i) = StatisticalComplexity_Structural(Arg);
    dGhg_n(i) = StatisticalComplexity_Structural(Ahg);
    dGrh_n(i) = StatisticalComplexity_Structural(Arh);
end