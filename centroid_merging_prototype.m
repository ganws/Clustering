%% Xmeans

% load data
X1 = readmatrix("6class.txt");
X2 = readmatrix("11class.txt");
X = [X1;X2];

% perfrom xmeans
k_max = 30;
[idx, centroids] = xmeans(X, 20);
k = size(centroids,1);

% convert index into cell index
%cluster_idx = cell(1,length(unique(idx)));
idx_cluster = {};
for i = unique(idx)'
    idx_cluster{i} = find(idx == i);
end
% plot result
figure
gscatter(X(:,1), X(:,2), idx);
hold on
plot(centroids(:,1), centroids(:,2), "kx");
title("Before")

%% Plot sigma eclipse 
% 
% sigma = zeros(k,2);
% for i = 1:size(centroids,1)
%     sigma(i,:) = std(X(idx==i,:));
% 
%     % plot sigma eclipse
%     rx=sigma(i,1); % horizontal radius
%     ry=sigma(i,2); % vertical radius
%     
%     x0=centroids(i,1); % x0,y0 ellipse centre coordinates
%     y0=centroids(i,2);
%     
%     t=-pi:0.01:pi;
%     x=x0+rx*cos(t);
%     y=y0+ry*sin(t);
%     plot(x,y)
%     
% end

%% Isolate centroid to study
% 
% id1 = 2;
% id2 = 3;
% id3 = 7;
% 
% sz1 = length(find(idx==id1));
% sz2 = length(find(idx==id2));
% sz3 = length(find(idx==id3));
% 
% idx2 = ones(sz1+sz2+sz3, 1);
% idx2(sz1+1 : sz1 + sz2) = 2;
% idx2(sz1 + sz2 +1 : sz1 + sz2 + sz3) = 3;
% 
% idx_cell{1} = 1:sz1;
% idx_cell{2} = sz1+1:sz1+sz2;
% idx_cell{3} = sz1+sz2+1:sz1+sz2+sz3;
% 
% X2 = [X(idx==id1,:);X(idx==id2,:);X(idx==id3,:)];
% c2 = [centroids(id1,:);centroids(id2,:);centroids(id3,:)]
% 
% cluster_num = 3;


%% Pairwise BIC Calculation

% k = 3;
C = centroids;

% Initialize
Xnew = X;
% C_new = C;
% idx_cluster = idx_cell;

C_Remaining_List = [1:k];
Pairwise_DecisionMatrix = false(k,k);
Merge_List = cell(k,2);

% row1: cluster index, row2: to-merge-neighbour index
for i = 1:k
    Merge_List{i,1} = i;
end

while length(C_Remaining_List) > 1
    
    current_centroid_idx = C_Remaining_List(1); % 
    C_Remaining_List(find(C_Remaining_List == current_centroid_idx)) = []; % remove current centroid from the list
    
    remaining_k = length(C_Remaining_List);
    
    % find rank neighbours by distance
    
    Closest_Neighbor_idx = knnsearch( C(C_Remaining_List,:), C(current_centroid_idx,:),'k', remaining_k);
    
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
        
        C_separated = [C(current_centroid_idx,:);C(id,:)];
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

%% Merge centroids

group_num = 0;

while ~isempty(Merge_List)
    group_num = group_num  + 1;
   
    Q = Merge_List{1,1};
    Q_tested = [];
    
    while ~isempty(Q)
        
        currentNode = Q(1);
        neighbors = ListNeighbors(Merge_List, currentNode);
        
        % add neighbor to test group if it is not in Q and Q_tested
        for nb = neighbors
            if isempty(find(Q == nb)) && isempty(find(Q_tested == nb))
                Q = [Q, nb];
            end
        end
        
        Q(1) = [];
        Q_tested = [Q_tested, currentNode];
        rmv_idx = FindIndexInList(Merge_List, currentNode);
        Merge_List(rmv_idx,:) = [];
    end
    
    Group{group_num} = Q_tested;
end

k_new = length(Group);

% New result
idx_new = ones(size(X,1),1);
centroids_new = zeros(k_new,2);

for g = 1:k_new
    for c = Group{g}
        idx_new(idx_cluster{c}) = g;
        centroids_new(g,:) = mean(C(c, :),1);
    end
end
%% Plot merge results

figure;
gscatter(X(:,1), X(:,2), idx_new);
hold on
plot(centroids_new(:,1), centroids_new(:,2), 'kx');
title("After")

% improve struct
[idx_new, centroids_new] = kmeans(X, k_new, 'Start', centroids_new);
figure;
gscatter(X(:,1), X(:,2), idx_new);
hold on
plot(centroids_new(:,1), centroids_new(:,2), 'kx');
title("After Structure-improve")

%% function

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


