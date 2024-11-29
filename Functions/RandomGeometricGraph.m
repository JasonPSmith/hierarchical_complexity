function [W,coords] = RandomGeometricGraph(n,q)

coords = rand(n,q);
W = squareform(pdist(coords));
W = ~eye(n).*((1-W/sqrt(q)));