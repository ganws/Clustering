% determine the effect of class size has towards performance of kmeans and
% kmeans

clear;clc

X1 = load('2classlabelled_case1.mat'); % NG 1%
X2 = load('2classlabelled_case2.mat'); % NG 2%
X3 = load('2classlabelled_case3.mat'); % NG 3%
X4 = load('2classlabelled_case4.mat'); % NG 4%
X5 = load('2classlabelled_case5.mat'); % NG 5%
X6 = load('2classlabelled_case6.mat'); % NG 6%
X7 = load('2classlabelled_case7.mat'); % NG 7%
X8 = load('2classlabelled_case8.mat'); % NG 8%
X9 = load('2classlabelled_case9.mat'); % NG 9%

X10 = load('2classlabelled_case10.mat'); % 10%
X30 = load('2classlabelled_case30.mat'); % 30%
X50 = load('2classlabelled_case50.mat'); % 50%

%% 

% [idx, C] = xmeans(X1.X, 5, 'bic');
X = X10.X;
P = 100;
N = 1000-P;

[idx, C] = xmeans(X, 5);

figure
plot(X(idx==1, 1), X(idx==1, 2), 'r.')
hold on
plot(X(idx==2, 1), X(idx==2, 2), 'b.')
plot(C(:,1), C(:,2), 'kx')

% ng accuracy
TNR = sum(idx(1:N) == 1)/ N
TPR = sum(idx( (N+1) : end) == 2) / P

