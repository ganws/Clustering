function [MNDLScore] = calculate_MNDL(X, clusters, centroids)
        % Calculates splitting criterion for input clusters using minimum noiseless description length criterion.
         
        %INPUT: 
        % clusters: Clusters for which splitting criterion should be calculated.
        % centroids: Centers of the clusters.
         
        % Returns splitting criterion in line with bayesian information criterion. 
        % Low value of splitting cretion means that current structure is much better.
         
         MNDLScore = Inf;
         W = 0;
         K = length(clusters);
         N = 0;
         
         sigma_sqrt = 0;
         
         alpha = 0.9;
         beta = 0.9;
         
         for k = 1:length(clusters)
             Ni = length(clusters{k});
             if Ni == 0 
                 MNDLScore = Inf;
                 return
             end
             
             % calculate eucledian distance
             V = X(clusters{k}, : ) - centroids(k,:);
             euclidean_dist = sum(sqrt(sum(V.^2, 2)));
             
             sigma_sqrt = sigma_sqrt + euclidean_dist;
             W = W + euclidean_dist / Ni;
             N = N + Ni;
         end
         
         if N - K > 0
            sigma_sqrt = sigma_sqrt / (N - K);
            sigma = sqrt(sigma_sqrt);
            Kw = (1.0 - K / N) * sigma_sqrt;
            Ks = ( 2.0 * alpha * sigma / (N^0.5) ) * ( (alpha^2) * sigma_sqrt / N + W - Kw / 2 ) ^0.5;
            MNDLScore = sigma_sqrt * (2 * K)^0.5 * ((2 * K)^0.5 + beta) / N + W - sigma_sqrt + Ks + 2 * alpha^0.5 * sigma_sqrt / N ;
         end
            
             
    end