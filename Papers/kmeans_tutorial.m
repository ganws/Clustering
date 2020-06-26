%% load data

load fisheriris;
X = meas(:, 3:4);

%% cluster data

rng(1); % for reproducibility
[idx,C] = kmeans(X,3);

%% visualizing Voronoi

x1 = min(X(:,1)):0.01:max(X(:,1));
x2 = min(X(:,2)):0.01:max(X(:,2));
[x1G, x2G] = meshgrid(x1,x2);
XGrid = [x1G(:), x2G(:)];

idx2Region = kmeans(XGrid, 3, 'MaxIter', 1, 'Start', C);

% plot the cluster regions
figure;
gscatter(XGrid(:,1), XGrid(:,2), idx2Region, ...
    [0,0.75, 0.75; 0.75, 0, 0.75; 0, 0.75, 0]);

hold on;

plot(X(:,1), X(:,2), 'k*', 'MarkerSize', 5);
title('Fisher''s Iris Data');
xlabel('Petal Lengths (cm)');
ylabel('Petal Widths (cm)');    



