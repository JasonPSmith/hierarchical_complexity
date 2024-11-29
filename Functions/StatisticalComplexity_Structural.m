function dG = StatisticalComplexity_Structural(A)

n = size(A,1);

a = modularity_louvain_und(A);
amod = length(unique(a))/n;
vmod = var(groupcounts(a'))/mean(groupcounts(a'));

L = eig(diag(diag(sum(A)))-A);
vlam = var(L)/mean(L);

rmot = sum(motif3struct_bin(A))/sum(motif4struct_bin(A));

dG = amod*rmot/vmod/vlam;

