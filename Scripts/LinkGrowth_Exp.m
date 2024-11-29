for k = 1:100
    n = 1000;
    d(k) = rand()*0.1;
    A = threshold_proportional(WeightedNetworkModel(n,'go3','log',0.2),d(k))>0;

    m = 100; %round(0.002*N*(N-1)/2/20);

%    for p = 1:10
        Ahier = A;
        Ahs = Ahier;
        Asim = Ahs;
        %Arnd = Asim;

        for j = 1:50
            % Pref
            Krep = repmat(sum(Ahier),size(Ahier,1),1);
            Krep = (Krep + Krep').*(~Ahier);
            Krep = Krep/(n-1);
            Krep(Krep>0) = exp(Krep(Krep>0));
            Phier = squareform(Krep-diag(diag(Krep)))./sum(squareform(Krep-diag(diag(Krep))));
            Phier = Phier./sum(Phier);

            Ahier = squareform(Ahier);
            id = randsample(1:length(Phier),m,true,Phier);
            Ahier(id) = 1;
            Ahier = squareform(Ahier);
            Rexphg_hier(k,j) = NormalisedHierarchicalComplexity(Ahier);

            % Random
            %Arnd = squareform(Arnd);
            %[~,id] = find(Arnd);
            %Arnd(id(randi(length(id),m,1))) = 0;
            %Arnd = squareform(Arnd);
            %Rhg_rnd(k,j)  = NormalisedHierarchicalComplexity(Arnd);
        end

        % Sim
        for j = 1:50
            Isim = zeros(size(Asim,1));
            for i = 1:size(Asim,1)
                Isim(:,i) = Asim(i,:)*Asim;
            end
            Isim = Isim - diag(diag(Isim));
            Isim(isnan(Isim)) = 0;
            Isim = Isim.*(~Asim);

            if sum(Isim(:)>0)>2*m
                Krep = repmat(sum(Asim),size(Asim,1),1);
                Krep = (Krep + Krep');%.*Asim;
                J = (Isim./Krep);%.*Asim;
                J = J - diag(diag(J));
                J(isnan(J)) = 0;
                J(isinf(J)) = 0;
                Psim = squareform(J);
                Psim(Psim>0) = exp(Psim(Psim>0));
                Psim = Psim./sum(Psim);
                Asim = squareform(Asim);

                id = randsample(1:length(Psim),m,true,Psim);
                Asim(id) = 1;

                Asim = squareform(Asim);
                Asim = Asim(sum(Asim)>0,sum(Asim)>0);
                Rexphg_sim(k,j)  = NormalisedHierarchicalComplexity(Asim);
            end
        end

        % Pref + Sim
        for j = 1:50
            Ihs = zeros(size(Ahs,1));
            for i = 1:size(Ahs,1)
                Ihs(:,i) = Ahs(i,:)*Ahs;
            end
            Ihs = Ihs - diag(diag(Ihs));
            Ihs(isnan(Ihs)) = 0;
            Ihs = Ihs.*(~Ahs);
            if sum(Ihs(:)>0)>2*m
                Phs = squareform(Ihs);
                Phs(Phs>0) = exp(Phs(Phs>0));
                Phs = Phs./sum(Phs);
                Ahs = squareform(Ahs);
                id = randsample(1:length(Phs),m,true,Phs);
                Ahs(id) = 1;
                Ahs = squareform(Ahs);

                Rexphg_hs(k,j)   = NormalisedHierarchicalComplexity(Ahs);
            end
        end
end