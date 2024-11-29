function [Rm,R] = NormalisedHierarchicalComplexity(A,directed)

% This function allows you to compute different hierarchical complexity
% normalisations. 
% INPUT
%   A - a symmetric adjacency matrix for an undirected graph
%   directed      - 1 provides directed HC measures: in-in, in-out, out-in, out-out  
%
% citation: Smith & Smith, Statistical Complexity of Heterogeneous
% Geometric Networks, PLOS Complex Systems, in press (2024). ArXiv: 
% https://doi.org/10.48550/arXiv.2310.20354




if ~exist('directed','var') || isempty(directed)
    directed = 0;
end

n = size(A,1);

if directed == 0
    
    K = sum(A);
    Kmin = min(K(K>0));
    Kmax = max(K);
    
    % Create neighbourhood degree sequences
    D = repmat(K',1,n);
    
    if n > 15000
        D = reshape(A,1,[]).*reshape(D,1,[]);
        D = reshape(D,n,n);
        D = sort(D);
        D(D==0) = nan;
    else
        D = sort(A.*D);  
        D(D==0) = nan;
    end
    
    R = nan(n,1);
    
    for i = Kmin:Kmax    % runs for each degree in the network
        if sum(K==i)>1   % makes sure to only work on degrees shared by more than 1 node
        R(i) = sum(std(D(:,K==i),0,2,'omitnan'),'omitnan');
        else
            R(i) = nan;
        end
    end
    
    R = R/(1-sum(K)/n/(n-1))/sum(K)*2;
    
    Rm = mean(R,'omitnan');

elseif directed == 1
   
    Kin = sum(A);
    Kout = sum(A');

    Kinmin = min(Kin(Kin>0));
    Kinmax = max(Kin);
    Koutmin = min(Kout(Kout>0));
    Koutmax = max(Kout);
    
    % Create neighbourhood degree sequences
    Din = repmat(Kin',1,n);
    Anan = double(A); 
    Anan(Anan==0) = nan;
    Din = sort(Anan.*Din);
        
    Dout = repmat(Kout',1,n);
    Dout = sort((Anan').*Dout);
    
    Rinin = nan(n,1);
    Routout = nan(n,1);
    Rinout = nan(n,1);
    Routin = nan(n,1);

    for i = Kinmin:Kinmax    % runs for each degree in the network
        if sum(Kin==i)>1   % makes sure to only work on degrees shared by more than 1 node
        Rinin(i) = sum(std(Din(:,Kin==i),0,2,'omitnan'),'omitnan');
        Rinout(i) = sum(std(Dout(:,Kin==i),0,2,'omitnan'),'omitnan');
        else
            Rinin(i) = nan;
            Rinout(i) = nan;
        end
    end

    for i = Koutmin:Koutmax    % runs for each degree in the network
        if sum(Kout==i)>1   % makes sure to only work on degrees shared by more than 1 node
        Routin(i) = sum(std(Din(:,Kout==i),0,2,'omitnan'),'omitnan');
        Routout(i) = sum(std(Dout(:,Kout==i),0,2,'omitnan'),'omitnan');
        else
            Routin(i) = nan;
            Routout(i) = nan;
        end
    end
    
    Rinin = Rinin/(1-sum(Kin)/n/(n-1))/sum(Kin)*2;
    Rinout = Rinout/(1-sum(Kin)/n/(n-1))/sum(Kin)*2;
    Routin = Routin/(1-sum(Kout)/n/(n-1))/sum(Kout)*2;
    Routout = Routout/(1-sum(Kout)/n/(n-1))/sum(Kout)*2;

    R = [Rinin;Rinout;Routin;Routout];
    
    Rm = [mean(Rinin,'omitnan'); mean(Rinout,'omitnan'); mean(Routin,'omitnan'); mean(Routout,'omitnan')];
end

