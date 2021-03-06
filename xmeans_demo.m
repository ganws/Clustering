% Demo code to perform xmeans clustering on sample data

% 20200716 - Use gscatter
% 20200714 - Written
% Written by Gan Wei Sheng

clear;clc;

X = readmatrix('5class.txt');
k_max = 10; % maximum allocation of cluster number

% Perform x-means on sample data
[idx, centers, wce]  = xmeans(X, k_max, 'bic', 'visualize_split', 'on'); 
result_k = length(unique(idx));

% Plot original data
figure
scatter(X(:,1),X(:,2), '.')
title("Before clustering")

% Visualize final clustering results
figure
h1 = gscatter(X(:,1), X(:,2), idx);
hold on
h2 = plot(centers(:,1), centers(:,2), 'kx');
title(['After xmeans, k = ', num2str(result_k)]);
legend('off')