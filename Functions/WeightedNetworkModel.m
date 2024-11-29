function A = WeightedNetworkModel(n,foundation,hierarchy,sigma,mu)

% Produces random modular hierarchy network with n nodes and k levels. 
%
%   INPUT:      n-              Number of nodes in the network
%               foundation-     enter: 'uni' for random graph, 'go2' for
%                               random geometric graph in 2D or 'go3' in 3D, 'sph' for
%                               random spherical manifold graph, 'hed' for head shape
%               hierarchy-      enter: 'log' for log-normal distribution,
%                               'exp' for exponential distribution and 'nrm' for normal
%                               distribution, 'noh' for no hierarchy
%               fitting-        'rnd for random fitting or 'ctr'
%                               for centrality fitting
%               sigma-          parameter sigma for normal, exponential and log-normal                               
%               mu-             parameter mu for log-normal distribution
%   OUTPUT:     A-              Complete weighted adjacency matrix for the weighted network model.
%
% Citation: Smith, K., Explaining the emergence of complex networks through log-normal fitness in a Euclidean Node Similarity Space, Scientific Reports, 11: 1976 (2021)

%
%
%  23-10-2018

rng('shuffle');

if ~exist('mu','var') || isempty(mu)
    mu = 0;
end

if all(foundation(1:3) == 'uni')
    R = squareform(rand(1,n*(n-1)/2));
elseif all(foundation(1:2) == 'go')
    R = RandomGeometricGraph(n,str2num(foundation(3:end)));
elseif all(foundation(1:2) == 'sp')
    m = str2num(foundation(3:end));
    for i = 1:n
        x = randn(m,1);
        coord(:,i) = x/norm(x);
    end

    [X,Y] = find(tril(ones(n),-1));
    R = ~eye(n).*(1-squareform(acos(sum(coord(:,X).*coord(:,Y)))));
elseif all(foundation(1:2) == 'hb')
    m = str2num(foundation(3:end));
    for i = 1:n
        x = randn(m,1);
        coord(:,i) = x/norm(x);
    end

    [X,Y] = find(tril(ones(n),-1));
    R = ~eye(n).*(1-squareform(acos(sum(coord(:,X).*coord(:,Y)))));
    
    %X(:,1) = randi(180,n,1)-90;
    %X(:,2) = randi(180,n,1);
    %Y = tril(ones(n),-1);
    %[J2,J1] = find(Y == 1);
    %R = ~eye(n).*(1-squareform(distance(X(J1,1),X(J1,2),X(J2,1),X(J2,2)))/180);
    %X = randn(n,3);
	%for i = 1:n
    %    X(i,:) = X(i,:)/norm(X(i,:));
    %end
    %R = ~eye(n).*(2- squareform(pdist(X)));
elseif all(foundation(1:3) == 'hed')
    X(:,1) = randi(90,2*n,1);
    X(:,2) = randi(180,2*n,1);
    A = unique(X,'rows');
    X = A(randperm(length(A),n),:);
    S = oblateSpheroid;
    S.SemimajorAxis = 9;
    S.SemiminorAxis = 6.5;
    Y = tril(ones(n),-1);
    [J2,J1] = find(Y == 1);
    R = ~eye(n).*(1-squareform(distance(X(J1,1),X(J1,2),X(J2,1),X(J2,2),S,'degrees')/9/pi));
    %[x,y,z]=sph2cart(pi*rand(1,n),2*pi*rand(1,n),ones(1,n));
    %R = squareform(pdist(0.5*[x;y;abs(z)]'));
end

% if all(fitting == 'ctr')
%     [~,index] = sort(sum(R),'descend');
%     R = R(index,index);
%     
%     if all(hierarchy == 'exp')
%         X = repmat(sort(exprnd(sigma,1,n),'descend'),n,1);
%         Hierval = ~eye(n).*(X + X');
%         A = R.*Hierval;
%         A = ~eye(n).*(A);
%     elseif all(hierarchy == 'nrm')
%         X = repmat(sort(normrnd(0,sigma,1,n),'descend'),n,1);
%         Hierval = ~eye(n).*(X + X');
%         A = R.*Hierval;
%         A = ~eye(n).*(A + abs(min(A(:))));
%     elseif all(hierarchy == 'log')
%         X = repmat(sort(lognrnd(mu,sigma,1,n),'descend'),n,1);
%         Hierval = ~eye(n).*(X + X');
%         A = R.*Hierval;
%         A = ~eye(n).*(A + abs(min(A(:))));
%     elseif all(hierarchy == 'noh')
%         A = R;
%    end
%elseif all(fitting == 'rnd')

if all(hierarchy == 'exp')
    X = repmat(exprnd(sigma,1,n),n,1);
    Hierval = ~eye(n).*(X + X');
    A = R.*Hierval;
    A = ~eye(n).*(A);
elseif all(hierarchy == 'nrm')
    X = repmat(normrnd(0,sigma,1,n),n,1);
    Hierval = ~eye(n).*(X + X');
    A = R.*Hierval;
    A = ~eye(n).*(A + abs(min(A(:))));
elseif all(hierarchy == 'log')
    X = repmat(lognrnd(mu,sigma,1,n),n,1);
    Hierval = ~eye(n).*(X + X');
    A = R.*Hierval;
    A = ~eye(n).*(A + abs(min(A(:))));
elseif all(hierarchy == 'pwr')
    X = repmat(rand(1,n).^(sigma),n,1);
    Hierval = ~eye(n).*(X + X');
    A = R.*Hierval;
    A = ~eye(n).*(A);
elseif all(hierarchy == 'noh')
    A = R;
end
%end

[~,index] = sort(sum(A),'descend');
A = A(index,index);
%A = A/max(A(:));