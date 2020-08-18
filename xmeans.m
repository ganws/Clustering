function [idx, C, total_wce] = xmeans(X, k_max, varargin)

%     MATLAB implmentation of xmeans algorithm based on Pyclustering version 0.9.3.1
%     20200721 - Improved split visualization
%     20200716 - Added mndl and gap as criterion
%              - Splitting criterion as input argument
%              - Added option to visualize splitting process
%     20200710 - Written
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
    iter = 1;
    
    
   % =============== Main algorithm steps =============
   while size(C, 1) <= k_max
       current_cluster_num = size(C, 1);
       is_split = [];
       [clusters, C, total_wce] = improve_param(C, null(1,1));
       [allocated_C, parentchild_relation, unborn_c] = improve_structure(clusters, C);
       C_old = C;
       
       if current_cluster_num == size(allocated_C,1)
           break
       else
           C = allocated_C;
       end
       
       [clusters, C, total_wce] = improve_param(C, null(1,1));
       idx = ConvertCellIdx2Idx(clusters);
       
       if strcmp(par.Results.visualize_split, 'on')
            visualize_splitting_process;
       end
       
       iter = iter+1;
   end
   
   % NOTE: 
   % Return idx instead clusters cell index as standard MATLAB practice as in kmeans.
   % Throughout the algorithm kmeans result index is converted to cell
   % index for convenience purpose because the code is based directly on
   % pyclustering implementation.
   % Refactoring needed to avoid unneceesary conversion
   
  % display output
  disp(['Splitting criterion = ', SPLITTING_TYPE]);
  disp(['xmeans result, k = ', num2str(length(clusters))]);
   
  
  
  
   % ================Nested function definitions==============
   % nested functions share the same workspace with main function so that 
   % copies of X do not need to be created
   
    function[clusters, local_centroids, local_wce] = improve_param(centroids, available_index)
        % Performs k-means clustering in specific region

        local_X = X;
        if ~isempty(available_index)
            local_X = X(available_index,:);
        end
            
        local_centroids = centroids;
        
        % Perform local kmeans
        if isempty(local_centroids)
            [clusterindx, local_centroids, local_sumD] = kmeans(local_X, 2);
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


    function[allocated_c, parent_child_relation, unborn_c]  = improve_structure(clusters, centers)
        % Check for best structure: divides each cluster into two and checks for best results using splitting criterion.

        % INPUT:
        % clusters: Clusters that have been allocated.
        % centers: Centroid coordinates.

        % return Allocated_c for clustering.
        
        allocated_c = [];
        unborn_c = [];
        parent_child_relation = []; % index showing if a parent has child or not
        amount_free_c = k_max - size(centers,1);
        
        for k = 1:length(clusters)
            % solve k-means problem for children where data of parent are used.
            [parent_child_clusters, parent_child_centers, ~] = improve_param(null(1,1), clusters{k});
            
            % if it is possible to split current cluster
            if length(parent_child_clusters) > 1 
                % Calculate splitting criterion
                parent_scores = splitting_criterion(clusters(k), centers(k,:));
                child_scores = splitting_criterion(parent_child_clusters, parent_child_centers);
                                
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
        
                if split_require && amount_free_c > 0
                    allocated_c = [allocated_c; parent_child_centers];
                    amount_free_c = amount_free_c - 1;
                    parent_child_relation = [parent_child_relation, k, k];
                    is_split = [is_split, true];
                else
                    unborn_c = [unborn_c; parent_child_centers];
                    allocated_c = [allocated_c; centers(k,:)];
                    parent_child_relation = [parent_child_relation,k];
                    is_split = [is_split, false];
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
             %euclidean_dist2 = sum(sum(V.^2, 2));

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

    function [] = visualize_splitting_process()
        figure
        gscatter(X(:,1), X(:,2), idx); % plot every points
        hold on
        
        % plot dead parents and their corresponding children
        dead_idx = find(is_split);
        for i = 1:length(dead_idx)
            
            p = C_old(dead_idx(i), :);  % parent
            c = C(parentchild_relation == dead_idx(i),:); %children
            
            h(1) = plot(p(:, 1), p(:, 2), 'ko', 'MarkerSize', 6);
            h(2) = plot(c(:, 1), c(:, 2), 'k.', 'MarkerSize', 12);
            
            %draw line between parent-children
            h(3) = plot([p(:,1), c(1,1)], [p(:,2), c(1,2)], 'k-');
            plot([p(:,1), c(2,1)], [p(:,2), c(2,2)], 'k-');
        end
        
        % plot alive parents and their corresponding unborn children
        alive_idx = find(~is_split);
        new_alive_idx = false(size(parentchild_relation));
        for a = 1:length(alive_idx)
            p = C(parentchild_relation==alive_idx(a),:); % parent
            c1 = unborn_c(2*a-1,:); % unborn 1
            c2 = unborn_c(2*a,:);   % unborn 2
            
            plot(p(:, 1), p(:, 2), 'k.', 'MarkerSize', 12);
            h(4) = plot(c1(:,1), c1(:,2), 'ko', 'MarkerSize', 3);
            plot(c2(:,1), c2(:,2), 'ko', 'MarkerSize', 3);
            
            plot([p(:,1), c1(1,1)], [p(:,2), c1(1,2)], 'k:');
            h(5) = plot([p(:,1), c2(1,1)], [p(:,2), c2(1,2)], 'k:');
        end
        hold off
        
        if(isempty(alive_idx))
            legend([h(2), h(1), h(3)],  'Centroids', 'Dead Centroid','Parent-Children Relation')
        else
            legend([ h(2), h(1), h(4), h(3), h(5)],  ...
                'Centroids', 'Dead Centroid','Unborn Children', 'Parent-Children Relation', 'Unborn Parent-Children Relation')
        end
        title(['Iter = ', num2str(iter), '  k = ', num2str(length(C))]);
        print(['iter-',num2str(iter)], '-dpng');
        
    end

    function [idx] = ConvertCellIdx2Idx(clusters)
        idx = zeros(size(X,1), 1);
        for m = 1:length(clusters)
            idx(clusters{m}) = m;
        end
    end
        

end

