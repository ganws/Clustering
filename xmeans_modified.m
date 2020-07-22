function [idx, C, total_wce] = xmeans_modified(X, k_max, varargin)

%     Modified version of xmeans
%     20200721 - Written
%     Written by Gan Wei Sheng
%     
%     INPUT:
%     X: [observations, dimensions] data matrix
%     k_max: Maximum number of clusters that can be allocated
%     SPLITTING_CRITERION (Optional): type of splitting criterion
%                                     -'bic' bayesian information criterion (default)
%                                     -'mndl' minimum noise description length
%                                     -'gap' gap statistic
%     'visualize_split' (Optional): visualize splitting process specified by name-value pair,
%                         value is 'on' or 'off'(default)
%                                     
%           
%     
%     OUTPUT:
%     idx: cluster index
%     C: [k,dimensions] matrix containing centroids coordinates, with rows corresponding to
%        calculated clusters
%     total_wce: Total euclidean squared distance


%    ////////////////____TODO_____///////////////
%    - Initial centroids as input argumenmt (current default is empty) 
%    - TOLERANCE as input argument
%    - Implement degree of class separation J as splitting criterion
%    ////////////////////////////////////////////

    % ============= Arguments validation ============
   par = inputParser;
   
   defaultSplitCrit = 'bic';
   validSplitCrit = {'bic', 'mndl', 'gap'};
   checkSplitCrit = @(x) any(validatestring(x,validSplitCrit));
   
   defaultSplitVisual = 'off';
   validSplitVisual = {'on', 'off'};
   checkSPlitVisual = @(x) any(validatestring(x, validSplitVisual));
   
   addRequired(par, 'X', @ismatrix);
   addRequired(par, 'k_max', @isnumeric);
   addOptional(par, 'SPLITTING_TYPE', defaultSplitCrit, checkSplitCrit);
   addParameter(par, 'visualize_split', defaultSplitVisual, checkSPlitVisual);
   
   parse(par, X, k_max, varargin{:});
   
    % ============== Initialization ============
  	[~, C] = kmeans(X,1); %initial 1 centroid
    TOLERANCE = 0.025;
    SPLITTING_TYPE = par.Results.SPLITTING_TYPE;
    repeat = 1; % repeat (unit): How many times k-means should be run to improve parameters (by default is 1).
                % With larger 'repeat' values suggesting higher probability of finding global optimum.
	currentSplittingIteration = 1;
    
   % =============== Main algorithm steps =============
   while size(C, 1) <= k_max
       current_cluster_num = size(C, 1);
       [clusters, C, total_wce] = improve_param(C, null(1,1), false);
       allocated_C = improve_structure(clusters, C);
       
       if current_cluster_num == size(allocated_C,1)
           break
       else
           C = allocated_C;
       end
       
       [clusters, C, total_wce] = improve_param(C, null(1,1), false);
   end
   
   % NOTE: 
   % Return idx instead clusters cell index as standard MATLAB practice as in kmeans.
   % Throughout the algorithm kmeans result index is converted to cell
   % index for convenience purpose because the code is based directly on
   % pyclustering implementation.
   % Refactoring needed to avoid unneceesary conversion
   
   idx = zeros(size(X,1), 1);
   for m = 1:length(clusters)
       idx(clusters{m}) = m;
   end
   
  % display output
  disp(['Splitting criterion = ', SPLITTING_TYPE]);
  disp(['xmeans result, k = ', num2str(length(clusters))]);
  
   % ================Nested function definitions==============
   % nested functions share the same workspace with main function so that 
   % copies of X do not need to be created
   
    function[clusters, local_centroids, local_wce] = improve_param(centroids, available_index, ONE_CLUSTER)
        % Performs k-means clustering in specific region
        % INPUT:
        % centroid
        % available_index
        % FIRST_CLUSTER: true(split 2), false(split3)

        local_X = X;
        if ~isempty(available_index)
            local_X = X(available_index,:);
        end
            
        local_centroids = centroids;
        
        % Perform local kmeans
        if isempty(local_centroids)
            if ONE_CLUSTER
                [clusterindx, local_centroids, local_sumD] = kmeans(local_X, 2, 'Replicates', 10);
            else
                [clusterindx, local_centroids, local_sumD] = kmeans(local_X, 3, 'Replicates', 10);
            end
        else
            [clusterindx, local_centroids, local_sumD] = kmeans(local_X, size(local_centroids,1), 'Start', local_centroids);
        end
        clusters = {};
        clusters = convertIdx2CellIdx(clusterindx);
        
        if ~isempty(available_index)
            clusters = local2globalclusters(clusters, available_index);
        end
        
        local_wce = sum(local_sumD);
        
    end


    function[allocated_c]  = improve_structure(clusters, centers)
        % Check for best structure: divides each cluster into two and checks for best results using splitting criterion.

        % INPUT:
        % clusters: Clusters that have been allocated.
        % centers: Centroid coordinates.

        % return Allocated_c for clustering.
        
        allocated_c = [];
        amount_free_c = k_max - size(centers,1);
        
        % Cluster pairing
        cluster_num =  length(clusters);
        unpaired_index = 1:cluster_num; %
        pair_index = {};
        pair_cluster_list = {};
        currentIter = 1;
        
        % assign a one-cluster 'pair ' if there is only 1 cluster or cluster
        % number is odd
        if cluster_num == 1 || mod(cluster_num, 2) == 1
            % choose the cluster furthest from the mean
            % pdist
            dist2 = sum((centers - mean(centers,1)).^2,2); %euclidean distance square
            [~, dist_idx] = sort(dist2);
            c1 = dist_idx(end);
            
            %c1 = randi([1, cluster_num]); % choose one from list (at this point all clusters are still available)
            pair_cluster_list{1} = {clusters{c1}}; % add chosen cluster into paired list
            pair_index{1} = c1;
            unpaired_index(c1) = []; % remove chosen index from list
            currentIter = currentIter + 1;
        end

        while ~isempty(unpaired_index)

            %select a random cluster from available index
             randIdx = randi(length(unpaired_index), 1);
             select_idx = unpaired_index(randIdx) ;
             unpaired_index(randIdx) = [];
             
             % find closest centroid from the list of unpaired centroids
             closest_cluster_idx = knnsearch(centers(unpaired_index, :), centers(select_idx,:), 'K', 1);
             neighbor_idx = unpaired_index(closest_cluster_idx);
             unpaired_index(closest_cluster_idx) = [];
             
             pair_cluster_list{currentIter} = {clusters{select_idx}, clusters{neighbor_idx}};
             pair_index{currentIter} = [select_idx, neighbor_idx];
             
             currentIter  = currentIter + 1;
        end
        
        for k = 1:length(pair_cluster_list)
            
            % search for the closest centroid
           querry_centers = centers(pair_index{k},:);
           querry_clusters = pair_cluster_list{k};
           ONE_CLUSTER = false;
           
           if length(pair_cluster_list{k}) == 1
               ONE_CLUSTER = true;
               querry_idx = querry_clusters{1}; % for one-cluster 'pair'
           else
               querry_idx = [querry_clusters{1};querry_clusters{2}];
           end
           
            % solve k-means problem for children where data of parent are used.
            [parent_child_clusters, parent_child_centers, ~] = improve_param(null(1,1), querry_idx, ONE_CLUSTER);
            
            % if it is possible to split current cluster
            if length(parent_child_clusters) > 1 
                % Calculate splitting criterion
                if ~strcmp(SPLITTING_TYPE, 'gap')
                    parent_scores = splitting_criterion(querry_clusters, querry_centers);
                    child_scores = splitting_criterion(parent_child_clusters, parent_child_centers);
                end
                                
                split_require = false;
                
                % Reallocate number of cluters in line with scores
                switch SPLITTING_TYPE
                    case 'bic' % BIC
                        if parent_scores < child_scores
                            split_require = true;
                        end
                        
                    case 'mndl' % MNDL
                        if parent_scores > child_scores
                            split_require = true;
                        end
                        
                    case 'gap' % Gap statistiac
                        eval = evalclusters(X(clusters{k},:), 'kmeans', 'gap', 'KList', 1:2, 'B', 10);
                        parent_scores = eval.CriterionValues(1);
                        child_scores = eval.CriterionValues(2);
                        if parent_scores < child_scores
                            split_require = true;
                        end
                end
                
                if strcmp(par.Results.visualize_split, 'on')
                % visualize splitting process
                    figure   
                    scatter(X(:,1), X(:,2), '.');
                    hold on
                    ccolor = {'g.','r.','b.'};
                    for pr = 1:length(parent_child_clusters)
                        scatter(X(parent_child_clusters{pr},1), X(parent_child_clusters{pr},2), ccolor{pr});
                        plot(parent_child_centers(pr,1), parent_child_centers(pr,2), 'k.', 'MarkerSize', 10);
                    end

                    plot(centers(:,1), centers(:,2), 'kx', 'MarkerSize', 10);
                    plot(querry_centers(:,1), querry_centers(:,2), 'kp', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
                    
                    anotStr = {['parent score=',num2str(parent_scores)], ['child score=', num2str(child_scores)]};
                    annotation('textbox', [.2 0 .3 .3], 'String', anotStr, 'FitBoxToText', 'on');
                    title(['iter=', num2str(currentSplittingIteration), '  split=', num2str(split_require)]);
                    print(['iter-',num2str(currentSplittingIteration)], '-dpng');
                end
                
                currentSplittingIteration = currentSplittingIteration + 1;
        
                if split_require && amount_free_c > 0
                    allocated_c = [allocated_c; parent_child_centers];
                    amount_free_c = amount_free_c - size(parent_child_centers,1);
                else
                    allocated_c = [allocated_c; centers(pair_index{k},:)];
                end 
            end
            
        end
        
        
    end


    function [score] = splitting_criterion(clusters, centers)
        % return splitting criterion score
        switch SPLITTING_TYPE
            case 'bic' % BIC
                score = bayesian_information_criterion(clusters, centers);
            case 'mndl' % MNDL
                score = minimum_noiseless_description_length(clusters, centers);
                
            otherwise % Default = Beyesian
                score = bayesian_information_criterion(clusters, centers);
                
        end     
    end

    
    function [clusters] = local2globalclusters(local_clusters, available_index)
        % Converts clusters in local region defined by 'available_index' to
        % global index

        clusters = cell(1,1);
        for i = 1:length(local_clusters)
            clusters{i} = available_index(local_clusters{i});
        end
    end


    function [cluster_idx] = convertIdx2CellIdx(idx)
        % convert index into cell index
        %cluster_idx = cell(1,length(unique(idx)));
        cluster_idx = {};
        for i = unique(idx)'
            cluster_idx{i} = find(idx == i);
        end
    end


    function [BICscore] = bayesian_information_criterion(indx_class, centroids)

        % Calculate bayesian information criterion (BIC) with Kass's formula:
        % BIC = L - 0.5 * p * ln(N)
        % L = n * log(n) - n * log(N) - n * 0.5 * log(2.0 * pi) - n * sigma_multiplier - (n - K) * 0.5

        % INPUT
        % X: dataset [observations, dimensions]
        % indx_class: index of each clusters in cell array
        % centroids: [K, dimensions] matrix contaiing centroid coordinates 

        % OUTPUT:
        % BICscore: sum of all bic for all clusters 
        
        K = length(indx_class); % number of cluster
        dimension = size(X,2);
        sigma_sqrt = 0;
        N = 0;

        % Maximum Likelihood Estimation for variance
        for i = 1:K
            Xtmp = X(indx_class{i},:);
            diff = Xtmp - centroids(i,:);
            dist_sqrt = sum(sum(diff.^2));
            sigma_sqrt = sigma_sqrt + dist_sqrt;
            N = N + size(Xtmp,1);
        end

        if N-K > 0

            sigma_sqrt = sigma_sqrt / (N-K);
            p = (K-1) + dimension * K + 1;
            
            % in case of the same points , sigma_sqrt can be zero
            if sigma_sqrt <= 0
                sigma_multiplier = -Inf;
            else
                sigma_multiplier = dimension * 0.5 * log(sigma_sqrt);
            end

            % calculate and sum BIC in all clusters
            score = zeros(K,1);
            for i = 1:K
                Xi = X(indx_class{i},:);
                n = size(Xi,1);
                L = n * log(n) ...
                    - n * log(N) ...
                    - n * 0.5 * log(2.0 * pi) ...
                    - n * sigma_multiplier ...
                    - (n - K) * 0.5 ;
                score(i) = L - p * 0.5 * log(N);
            end
        end

        BICscore = sum(score);
    end


    function [MNDLScore] = minimum_noiseless_description_length(clusters, centroids)
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
             euclidean_dist2 = sum(sum(V.^2, 2));

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

end

