function [BICscore] = calculateBIC(X, indx_class, centroids)

% Calculate bayesian information criterion (BIC) with Kass's formula:
    % BIC = L - 0.5 * p * ln(N)
    % L = n * log(n) - n * log(N) - n * 0.5 * log(2.0 * pi) - n * sigma_multiplier - (n - K) * 0.5

    % INPUT
    % X: dataset [observations, dimensions]
    % indx_class: index of each clusters in cell array
    % centroids: [K, dimensions] matrix contaiing centroid coordinates 

    % OUTPUT:
    % BICscore: sum of all bic for all clusters 

    K = length(indx_class); % number of cluster
    dimension = size(X,2);
    sigma_sqrt = 0;
    N = 0;

    % Maximum Likelihood Estimation for variance
    for i = 1:K
        Xtmp = X(indx_class{i},:);
        diff = Xtmp - centroids(i,:);
        dist_sqrt = sum(sum(diff.^2));
        sigma_sqrt = sigma_sqrt + dist_sqrt;
        N = N + size(Xtmp,1);
    end

    if N-K > 0

        sigma_sqrt = sigma_sqrt / (N-K);
        p = (K-1) + dimension * K + 1;

        % in case of the same points , sigma_sqrt can be zero
        if sigma_sqrt <= 0
            sigma_multiplier = -Inf;
        else
            sigma_multiplier = dimension * 0.5 * log(sigma_sqrt);
        end

        % calculate and sum BIC in all clusters
        score = zeros(K,1);
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
    end

    BICscore = sum(score);

end
