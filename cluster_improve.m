function [idx_new, C_new] = cluster_improve(X, idx_old, C_old)

%     Algorithm to improve clustering results
%     20200727 - Written
%     Written by Gan Wei Sheng
%     
%     INPUT:
%     X: [observations, dimensions] data matrix
%     idx_old: old cluster index
%     C_old: old centroid coordinates

%     OUTPUT:
%     idx: new cluster index
%     C_new: new centroid coordinates

%    ////////////////____TODO_____///////////////
%    - In BIC calculation, calculate from the closest neighbour first 
%      doesnt make any difference if pairwise BIC is calculated
%    ////////////////////////////////////////////

    % ============== Initialization ==================
    k = size(C_old,1);
    C_Remaining_List = [1:k];
    Pairwise_DecisionMatrix = false(k,k);
    Merge_List = cell(k,2); % row1: cluster index, row2: to-merge-neighbour index
    for i = 1:k
        Merge_List{i,1} = i;
    end
    
    idx_cluster = {};
    for i = unique(idx_old)'
        idx_cluster{i} = find(idx_old == i);
    end

    % ============== Pairwise BIC Calculation ==================
    while length(C_Remaining_List) > 1

        current_centroid_idx = C_Remaining_List(1); % 
        C_Remaining_List(find(C_Remaining_List == current_centroid_idx)) = []; % remove current centroid from the list

        remaining_k = length(C_Remaining_List);

        % find rank neighbours by distance

        Closest_Neighbor_idx = knnsearch( C_old(C_Remaining_List,:), C_old(current_centroid_idx,:),'k', remaining_k);

        % test every centroid-neighbor pair bic, determine from BIC if they
        % should merge or just leave it

    %     toJoin = false(remaining_k,1); % array to store the result if centroid-neighbor pair should join or not
        iter = 1;

        for id = C_Remaining_List(Closest_Neighbor_idx)

            X1 = X(idx_cluster{current_centroid_idx}, :);
            X2 = X(idx_cluster{id}, :);
            Xtmp = [X1;X2];

            sz1 = size(X1,1);
            sz2 = size(X2,1);

            idx_separated = { 1:sz1 , sz1+1:sz1+sz2 };
            idx_merged = {1:sz1+sz2};

            C_separated = [C_old(current_centroid_idx,:);C_old(id,:)];
            C_merged = mean(C_separated);

            bic_separated = calculateBIC(Xtmp, idx_separated, C_separated);
            bic_merged = calculateBIC(Xtmp, idx_merged, C_merged);

            if bic_merged > bic_separated
                toMerge = true;
                Merge_List{current_centroid_idx,2} = [Merge_List{current_centroid_idx,2}, id];
                Merge_List{id,2} = [Merge_List{id,2}, current_centroid_idx];

            else
                toMerge = false;
            end

            Pairwise_DecisionMatrix(current_centroid_idx, id) = toMerge;
            Pairwise_DecisionMatrix(id, current_centroid_idx) = toMerge;

            fprintf("Current Centroid: %d | Current Neighbour: %d | bic_separated: %.2f | bic_joined: %.2f | Merge => %d\n", ...
                    current_centroid_idx, id, bic_separated, bic_merged, toMerge);

            iter = iter + 1;
        end

    end

    % ============= Merge Centroids ================
    % Each centroid is a node. Merged pair forms a chain link.
    % All connected nodes form a group.
    
    group_num = 0;

    while ~isempty(Merge_List)
        group_num = group_num  + 1;

        Q = Merge_List{1,1}; % nodes to be testd
        Q_tested = []; % nodes that are already tested
        
        % Q will be empty if all connected nodes are tested
        while ~isempty(Q)

            currentNode = Q(1);
            neighbors = ListNeighbors(Merge_List, currentNode);

            % add neighbor to test group if it is not in Q and Q_tested
            for nb = neighbors
                if isempty(find(Q == nb, 1)) && isempty(find(Q_tested == nb, 1))
                    Q = [Q, nb];
                end
            end

            Q(1) = []; % remove currentNode from Q
            Q_tested = [Q_tested, currentNode]; % add currentNode to Q_tested
            rmv_idx = FindIndexInList(Merge_List, currentNode); 
            Merge_List(rmv_idx,:) = []; % remove currentNode from Merge_List
        end

        Group{group_num} = Q_tested; % group all connected nodes
    end

    % New result
    k_new = length(Group);
    idx_new = ones(size(X,1),1);
    C_new = zeros(k_new,2);

    for g = 1:k_new
        for c = Group{g}
            idx_new(idx_cluster{c}) = g;
            [~, C_new(g,:)] = kmeans(X(idx_new==g, :), 1);
        end
    end
    
    % perform kmeans to improve structures
	[idx_new, C_new] = kmeans(X, k_new, 'Start', C_new);


end 

%% Local functions

function [neighbors] = ListNeighbors(List, CentroidIdx)

    for i=1:size(List,1)
        if List{i,1} == CentroidIdx
            neighbors = List{i,2};
        end
    end
end

function [index] = FindIndexInList(List, CentroidIdx)
    for i=1:size(List,1)
        if List{i,1} == CentroidIdx
            index = i;
        end
    end
end
