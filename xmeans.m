function [clusters, C, total_wce] = xmeans(X, k_max)

%     MATLAB implmentation of xmeans algorithm based on Pyclustering version 0.9.3.1
%     20200710 - Written
%     Written by Gan Wei Sheng
%     
%     INPUT:
%     X: [observations, dimensions] data matrix
%     k_max: Maximum number of clusters that can be allocated
%     
%     OUTPUT:
%     clusters: Cell array containing index for each cluster
%     C: [k,dimensions] matrix containing centroids coordinates, with rows corresponding to
%        calculated clusters
%     total_wce: Total euclidean squared distance

    

%    ////////////////____TODO_____///////////////
%    - Initial centroids as input argumenmt (current default is empty) 
%    - spltting type as input argument (current default = 0 (BIC) ))
%    - TOLERANCE as input argument
%    - Implement degree of class separation J as splitting criterion
%    ////////////////////////////////////////////

    

    % ============== Initialization ============
	C = [];
    SPLITTING_TYPE = 0; % 0 = bayesian information criterion, 1 = degree of class separation(?)
    TOLERANCE = 0.025;
    repeat = 1; % repeat (unit): How many times k-means should be run to improve parameters (by default is 1).
                % With larger 'repeat' values suggesting higher probability of finding global optimum.
                
                
   % =============== Main algorithm steps =============
   while size(C, 1) <= k_max
       current_cluster_num = size(C, 1);
       [clusters, C, ~] = improve_param(C, null(1,1));
       allocated_C = improve_structure(clusters, C);
       
       if current_cluster_num == size(allocated_C,1)
           break
       else
           C = allocated_C;
       end
       
       [clusters, C, total_wce] = improve_param(C, null(1,1));
   end
   
   
   % ================Nested function definitions==============
   
    function[clusters, local_centroids, local_wce] = improve_param(centroids, available_index)
        % Performs k-means clustering in specific region
        % param[in] s: [k x dimension] Cluster s matrix , if None then automatically generated two s using center initialization method.

        local_X = X;
        if ~isempty(available_index)
            local_X = X(available_index,:);
        end
            
        local_centroids = centroids;
        
        % Perform local kmeans
        if isempty(local_centroids)
            [clusterindx, local_centroids, local_sumD] = kmeans(local_X, 2);
        else
            [clusterindx, local_centroids, local_sumD] = kmeans(local_X, size(centroids,1), 'Start', centroids);
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
        amount_free_s = k_max - size(centers,1);
        
        for k = 1:length(clusters)
            % solve k-means problem for children where data of parent are used.
            [parent_child_clusters, parent_child_s, ~] = improve_param(null(1,1), clusters{k});
            
            % if it is possible to split current cluster
            if length(parent_child_clusters) > 1 
                % Calculate splitting criterion
                parent_scores = splitting_criterion(X, clusters(k), centers(k,:));
                child_scores = splitting_criterion(X, parent_child_clusters, parent_child_s);
                
                split_require = false;
                
                % Reallocate number of s (cluters) in line with
                % scores
                if SPLITTING_TYPE == 0 % Bayesian
                    if parent_scores < child_scores
                        split_require = true;
                    end
                end
                
                if split_require && amount_free_s > 0
                    allocated_c = [allocated_c; parent_child_s];
                    amount_free_s = amount_free_s - 1;
                else
                    allocated_c = [allocated_c; centers(k,:)];
                end 
            end
            
        end
        
        
    end


    function [score] = splitting_criterion(X, clusters, s)
        % return splitting criterion score
        switch SPLITTING_TYPE
            case 0 %BIC
                score = bayesian_information_criterion(X, clusters, s);
            case 1 %J
                error("J as splitting criterion is not yet implemented")
        end     
    end


    function [clusters] = local2globalclusters(local_clusters, available_index)
        % Converts clusters in local region defined by 'available index' to
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
end

function [BICscore] = bayesian_information_criterion(X, indx_class, centroids)

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
    
    % Maximum Likelihood Estimation for variance
    sigma_sqr = 0;
    N = 0;
    for i = 1:K
        Xtmp = X(indx_class{i},:);
        diff = Xtmp - centroids(i,:);
        dist_square = sum(sum(diff.^2));
        sigma_sqr = sigma_sqr + dist_square;
        N = N + size(Xtmp,1);
    end
    sigma_sqr = sigma_sqr / (N-K);
    
    p = (K-1) + dimension * K + 1;
    sigma_multiplier = dimension * 0.5 * log(sigma_sqr);
    
    score = zeros(K,1);
    
    % calculate and sum BIC in all clusters
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
    
BICscore = sum(score);

end