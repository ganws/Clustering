function [fig] = visualize_cluster(X, cluster_idx, centers)
    % visualize clusters up to 3d spaces
    
    % 20200715
    % Written by Gan Wei Sheng
    
    
    dim = size(centers,2);
    cluster_num = size(centers,1);
    if (dim ==2)
        figure
        for i = 1:length(cluster_idx)
            scatter(X(cluster_idx{i},1), X(cluster_idx{i},2), '.');
            hold on
            scatter(centers(i,1), centers(i,2), 'ko', 'filled');
            hold on
        end
    end
    
end