% Demo code to perform xmeans clustering on sample data

% 20200714 - Written
% Written by Gan Wei Sheng

clear;clc;

X = readmatrix('5class.txt'); % read data from file
k_max = 10; % maximum allocation of cluster number

% Perform x-means on sample data
[cluster_indx, centers, wce]  = xmeans(X, k_max); 

% Plot original data
figure
scatter(X(:,1),X(:,2), '.')
title("Before clustering")

% Visualize clustering results
dim = size(centers,2);
if (dim ==2)
    figure
    for i = 1:length(cluster_indx)
        scatter(X(cluster_indx{i},1), X(cluster_indx{i},2), '.');
        hold on
        scatter(centers(i,1), centers(i,2), 'ko', 'filled');
        hold on
    end
end
title(['After xmeans, k = ', num2str(length(cluster_indx))]);