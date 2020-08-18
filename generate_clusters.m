%% Script to manually generate clusters (2d) with normal distribution

% Gan Wei Sheng
% v20200625 - created
% v20200630 - change to multivariate normal distribution function

clear,clc
%% CLUSTER GENERATION
% x and y are set to be independent variables (Cov(x,y) = 0)

% Generation settings
cluster_num = 1;
plot_style = {'rx', 'bo', 'g^'};

sz = [1000, 1000, 1000]; 

% cluster1
sigma = [1,0 ; 0,1]; % covariance matrix
mu = [2,3]; 
X1 = mvnrnd(mu, sigma, sz(1));

% cluster2
sigma = [1,0; 0, 1]; % covariance matrix
mu = [4,4]; 
X2 = mvnrnd(mu, sigma, sz(2));

% cluster3
sigma = [1,0; 0, 1]; % covariance matrix
mu = [10,1]; 
X3 = mvnrnd(mu, sigma, sz(3));

X = [X1;X2;X3]; %combine all samples into 1 matrix
%% TAGGING

C = ones(sum(sz),1); %cluster label
C(sz(1)+1:sz(1)+sz(2),:) = 2;
C(sz(1)+sz(2)+1:end,:) = 3;

%generate class_index
for i  = 1:cluster_num
    indx_class{i} = find(C==i);
end

%% Plot

for i = 1:cluster_num
    plot(X(C==i,1), X(C==i,2), plot_style{i});
    hold on
end

hold off

% save('toy_1cluster.mat', 'X', 'C', 'sz', 'cluster_num', 'plot_style', 'indx_class');