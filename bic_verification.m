% read pyclustering xmeans sample result data saved in csv and cross check
% with matlab implementation

clear;clc

cluster_indx_csv = readmatrix("cluster_indx.csv", "Range", 1);
centers = readmatrix("centroid_info.csv");
X = csvread("sample_data.csv");
finalBIC = -143.9659; %result taken directly from pyclustering results

cluster_num = size(centers, 1);

%convert cluster_indx into matlab cell format
cluster_indx = cell(cluster_num,1);
for cluster = 1:cluster_num
    for element = 1: size(cluster_indx_csv,2)
        tmp = cluster_indx_csv(cluster, element);
        if ~isnan(tmp)
            cluster_indx{cluster}(:,end+1) = tmp+1 % matlab index starts from 1
        end
    end
end

%plot clusters and centers
figure
ccolor = {'g', 'r', 'b', 'c'};
for i = 1:cluster_num
    scatter(X(cluster_indx{i},1), X(cluster_indx{i},2), ccolor{i});
    hold on
    scatter(centers(i,1), centers(i,2), 100, 'filled',  ccolor{i});
    hold on
end

% check BIC
matlab_bic = calculateBIC(X, cluster_indx, centers)
mndl = calculateMNDL(X, cluster_indx, centers)