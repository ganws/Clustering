function [BICscore] = calculateBIC(X, indx_class, centroids)

    % Calculate bayesian information criterion (BIC) with Kass's formula:
    % BIC = L - 0.5 * p * ln(N)
    % L = n * log(n) - n * log(N) - n * 0.5 * log(2.0 * pi) - n * sigma_multiplier - (n - K) * 0.5

    % 20200715 - Bug-fixes
    % 20200706 - Written
    % Written by Gan Wei Sheng
    
    % INPUT
    % X: dataset [observations, dimensions]
    % indx_class: index of each clusters in cell array
    % centroids: [K, dimensions] matrix contaiing centroid coordinates 
    
    % OUTPUT:
    % BICscore: sum of all bic for all clusters 

    
    % ========== Check parameters =============
    
    if length(indx_class) ~= size(centroids,1)
        error("number of k must be the same for both indx_class and centroids")
    end
    
    % =========== BIC Calculation ==============
    K = length(indx_class); % number of cluster
    dimension = size(X,2);
    
    % Maximum Likelihood Estimation for variance
    sigma_sqr = 0;
    N = 0;
    for i = 1:K
        Xtmp = X(indx_class{i},:);
        diff = Xtmp - centroids(i,:);
        dist_square = sum(sum(diff.^2));
        sigma_sqr = sigma_sqr + dist_square;
        N = N + size(Xtmp,1);
    end
    sigma_sqr = sigma_sqr / (N-K);
    
    p = (K-1) + dimension * K + 1;
    sigma_multiplier = dimension * 0.5 * log(sigma_sqr);
    
    score = zeros(K,1);
    
    % calculate and sum BIC in all clusters
    for i = 1:K
        Xi = X(indx_class{i},:);
        n = size(Xi,1);
        L = n * log(n) ...
            - n * log(N) ...
            - n * 0.5 * log(2.0 * pi) ...
            - n * sigma_multiplier ...
            - (n - K) * 0.5 ;
        score(i) = L - p * 0.5 * log(N);
    end
    
BICscore = sum(score);

end
