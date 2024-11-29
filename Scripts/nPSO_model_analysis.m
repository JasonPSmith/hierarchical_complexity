%% ALL PARAMETER TEST
N = 1000;                                           % Number of nodes 
m = [5,10,15,20,25];                                % Half of average degree            
T = [0,0.33,0.66,1];                                % Clustering
gamma = [1:0.25:3];                                 % Power-law exponent
distr = [0,1,3,5,7,9];                              % Number of communities?

for i = 1:length(m)
    d(i) = sum(sum(A))/N/(N-1);
    for j = 1:length(T)
        for k = 1:length(gamma)
            for l = 1:length(distr)
                for p = 1:100
                    A = full(nPSO_model(N,m(i),T(j),gamma(k),distr(l)));                   
                    NHC(i,j,k,l,p) = NormalisedHierarchicalComplexity(A);
                end
            end
        end
    end
end


%% GAMMA TEST
N = 1000;                                           % Number of nodes 
m = [5,10,15,20];                                   % Half of average degree            
T = 0;                                            % Clustering
gamma = [1:0.25:3];                                 % Power-law exponent
distr = 0;                                          % Number of communities?

for i = 1:length(m)    
    for k = 1:length(gamma)
        for p = 1:10
            A = full(nPSO_model(N,m(i),T,gamma(k),distr));       
            NHC_mg(i,k,p) = NormalisedHierarchicalComplexity(A);
            NDV_mg(i,k,p) = NormalisedDegreeVariance(A);
            if k == 1 & p == 1
                d_mg(i) = sum(sum(A))/N/(N-1);
            end
        end
    end
end

%% CLUSTERING TEST
N = 1000;                                           % Number of nodes 
m = 5;                                              % Half of average degree            
T = [0.1:0.1:1];                                    % Clustering
gamma = 2;                                          % Power-law exponent
distr = 0;                                          % Number of communities?

for i = 1:length(T)
    for p = 1:100
        A = full(nPSO_model(N,m,T(i),gamma,distr)); 
        NHC_T(i,p) = NormalisedHierarchicalComplexity(A);
        NDV_T(i,p) = NormalisedDegreeVariance(A);
    end
end

%% COMMUNITY TEST
N = 1000;                                           % Number of nodes 
m = 5;                                              % Half of average degree            
T = 0;                                              % Clustering
gamma = 2;                                          % Power-law exponent
distr = 1:20;                                       % Number of communities?

for i = 1:length(distr)
    for p = 1:100
        A = full(nPSO_model(N,m,T,gamma,distr(i))); 
        NHC_distr(i,p) = NormalisedHierarchicalComplexity(A);
        NDV_distr(i,p) = NormalisedDegreeVariance(A);
    end
end