clear;clc;
% 
X1 = readmatrix('11class.txt'); % read data from file
X2 = readmatrix('5class.txt');
% X =  readmatrix('sample_data.csv');
X = [X1;X2];
% % load toy_1cluster.mat
% load fisheriris.mat
% X = meas(:,3:4);
%%
% 
% ITERATION = 20;
% k_max = 30;
% % 
% % for iter = 1:ITERATION
% %     k_max = 30; % maximum number of clusters allowed
% %     [idx, centers, wce]  = xmeans_modified(X, k_max, 'bic', 'visualize_split', 'off'); % perform      
% %     k(iter) = length(centers);
% % end
% 

[idx, centers, wce]  = xmeans(X, 30, 'bic', 'visualize_split', 'off'); % perform      
% [idx, centers] = kmeans(X,13);
% visualize_cluster(X, cluster_indx,  centers);

% idx = dbscan(X, 0.5, 5);

% visualize
figure
gscatter(X(:,1), X(:,2), idx);
hold on
plot(centers(:,1), centers(:,2), 'kx');
hold off
legend("off")
title(['after xmeans, k=', num2str(length(unique(idx)))])


% 
% % improve clustering
idx_old = idx;
C_old = centers;
ITERATION = 5;
k = zeros(1, ITERATION+1);
k(1) = size(C_old,1);

for i = 1:ITERATION

    [idx_new, C_new] = cluster_improve(X, idx_old, C_old);
    k(i+1) = size(C_new, 1);
    
	% improve struct
    figure;
    gscatter(X(:,1), X(:,2), idx_new);
    hold on
    plot(C_new(:,1), C_new(:,2), 'kx');
    hold off
    title(['Improve Iter = ' ,num2str(i) , ', k = ', num2str(k(i+1))])
    
    idx_old = idx_new;
    C_old = C_new;

end

figure
plot(0:ITERATION, k)

