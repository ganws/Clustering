% compute and plot the pdf of a bivariate normal distribution
% with parameter mu = [0,0] and sigma = [0.25, 0.3, 0.3, 1];

% Define parameters mu and sigma
mu = [0,0];
sigma = [1, 0; 0, 1];

% Create grid of evenly spaced points in two-dimensional space
x1 = -3:0.2:3;
x2 = -3:0.2:3;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

% Evaluate the pdf of the normnal distribution at the grid points
y = mvnpdf(X, mu, sigma);
y = reshape(y,length(x2), length(x1));

% Plot the pdf values

surf(x1,x2,y)
caxis([min(y(:))-0.5*range(y(:)), max(y(:))]);
axis([-3 3 -3 3 0 0.4]);
xlabel('x1');
xlabel('x2');
zlabel('Probability Density');