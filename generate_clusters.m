%% Script to manually generate clusters (2d) with normal distribution

% Gan Wei Sheng
% v20200625 - created


%% CLUSTER GENERATION

cluster_num = 3;
plot_style = {'rx', 'bo', 'g^'};

% variable init
sz = zeros(cluster_num,1 ); 

% cluster1
sz(1) = 50; %num of samples
sigma = [1,1]; %standard deviation
mu = [0,0]; 
x1 = normrnd(mu(1), sigma(1), [sz(1),1]);
y1 = normrnd(mu(2), sigma(2), [sz(1),1]);

% cluster2
sz(2) = 50; %num of samples
sigma = [1,1]; %standard deviation
mu = [4,4]; 
x2 = normrnd(mu(1), sigma(1), [sz(2),1]);
y2 = normrnd(mu(2), sigma(2), [sz(2),1]);

% cluster3
sz(3) = 50; %num of samples
sigma = [1,1]; %standard deviation
mu = [-2,6]; 
x3 = normrnd(mu(1), sigma(1), [sz(3),1]);
y3 = normrnd(mu(2), sigma(2), [sz(3),1]);

X = [x1,y1;x2,y2; x3,y3]; %matrix representation
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

save('toy_cluster_case1.mat', 'X', 'C', 'sz', 'cluster_num', 'plot_style', 'indx_class');