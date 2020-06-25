%% Program to evaluate kmeans using J as evaluation index

% Gan Wei Sheng
% v20200625 - created

clear;clc;
load toy_cluster_case1.mat;
%% SETTINGS

max_k = 20; % maximum k (num of cluster) to be tested 

%% Main experiment loop

for k = 2:max_k
     [idx, centroid] = kmeans(X, k); %clustering
     
    %create new index based on clustering results
    for j = 1:k
        new_indx_class{j} = find(idx == j);
    end
    
J(k) = calculateJ(X', new_indx_class); % calculate degree of separation between clusters

end

%% Plot results

plot(1:max_k, J)
xlabel('k'); ylabel('J');