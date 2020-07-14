%% Program to evaluate kmeans using J as evaluation index

% Gan Wei Sheng
% v20200625 - created

clear;clc;
X1 = readmatrix('5class.txt');
X2 = readmatrix('11class.txt');
X = [X1;X2];
%% SETTINGS

max_k = 30; % maximum k (num of cluster) to be tested 
%% TEST
% 
% k = 5;
% [idx, centroid, sumD, D] = kmeans(X, k); %clustering
%     for j = 1:k
%         new_indx_class{j} = find(idx == j);
%     end
% Jtmp = calculateJ(X', new_indx_class); % calculate degree of separation between clusters
% 
% %plot


%% MAIN EXPERIMENT LOOP
ITERATION = 20;

ktest = length(2:max_k);
J_final= zeros(ITERATION, ktest);
distort_final = zeros(ITERATION, ktest);
BIC_final = zeros(ITERATION, ktest);

for iter = 1:ITERATION
    
    
    IK = 1;
    for k = 2:max_k
        
        [idx, centroid{k}, sumD, D] = kmeans(X, k); %clustering

        %create new index based on clustering results
        new_indx_class = {};
        for j = 1:k
            new_indx_class{j} = find(idx == j);
        end

        [J(IK), bSigma(IK), wSigma(IK)] = calculateJ(X', new_indx_class); % calculate degree of separation between clusters
        distortion(IK) = sum(sumD)/ size (X,1);
        bic(IK) = calculateBIC(X, new_indx_class, centroid{k});

        IK = IK +1;

        %save the last kmeans result from the last iteration for vizualization
        if iter == ITERATION
            cluster_indx{k} = new_indx_class;
        end

    end
        
    J_final(iter,:) = J;
    distort_final(iter,:) = distortion;
    BIC_final(iter,:) = bic;
   
end

J_mean = mean(J_final);
distort_mean = mean(distort_final);
bic_mean = mean(BIC_final);

J_std = std(J_final);
distort_std = std(distort_final);
bic_std = std(BIC_final);

%% Visualize clusters

plot_k = 16; % choose which result to plot

figure
for k = 1:plot_k
    scatter(X(cluster_indx{plot_k}{k},1), X(cluster_indx{plot_k}{k},2), '.')
    hold on
    scatter(centroid{plot_k}(k,1), centroid{plot_k}(k,2), 50, 'ko', 'filled')
end
textannot = ['BIC = ', num2str(bic(plot_k-1))]
annotation('textbox', [.15,.7,.2,.2], 'String', textannot, 'FitBoxToText', 'on');
hold off

%% 
% % degree of class separation J
% figure;
% yyaxis left
% plot(2:max_k, J)
% xlabel('k'); ylabel('J');
% hold on
% % plot(1:max_k, bSigma)
% % plot(1:max_k, wSigma)
% 
% % distortion
% yyaxis right
% plot(2:max_k, distortion)
% xlabel('k');
% ylabel('distortion per point');
% legend('J', 'distortion');
% 
% figure
% plot(2:max_k, bic);
% xlabel('k');
% ylabel('bic');

%% Plot experiment results

% degree of class separation
figure
yyaxis left
errorbar(2:max_k, J_mean, J_std);
xlabel("k")
ylabel("J")

hold on

% distortion
% figure
yyaxis right
errorbar(2:max_k, distort_mean, distort_std);
xlabel("k")
ylabel("distortion")
legend("J", "distortion");

% BIC
figure
errorbar(2:max_k, bic_mean, bic_std);
xlabel("k")
ylabel("BIC")


