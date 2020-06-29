%% Program to evaluate kmeans using J as evaluation index

% Gan Wei Sheng
% v20200625 - created

clear;clc;
load toy_cluster_case1.mat;
%% SETTINGS

max_k = 20; % maximum k (num of cluster) to be tested 

%% TEST
% 
% k = 5;
% [idx, centroid, sumD, D] = kmeans(X, k); %clustering
%     for j = 1:k
%         new_indx_class{j} = find(idx == j);
%     end
% Jtmp = calculateJ(X', new_indx_class); % calculate degree of separation between clusters
% 
% %plot


%% MAIN EXPERIMENT LOOP

for k = 2:max_k
     [idx, centroid, sumD, D] = kmeans(X, k); %clustering
     
    %create new index based on clustering results
    for j = 1:k
        new_indx_class{j} = find(idx == j);
    end
    
[J(k), bSigma(k), wSigma(k)] = calculateJ(X', new_indx_class); % calculate degree of separation between clusters
distortion(k) = sum(sumD)/ size (X,1);
end

%% PLOT

% degree of class separation J
figure;
plot(1:max_k, J)
xlabel('k'); ylabel('J');
hold on
plot(1:max_k, bSigma)
plot(1:max_k, wSigma)
legend('J', 'bSigma', 'wSigma');

% distortion
figure;
plot(1:max_k, distortion)
xlabel('k');
ylabel('distortion per point');