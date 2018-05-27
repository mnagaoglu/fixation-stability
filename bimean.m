function mu = bimean(X,Y,density)
% function mu = bimean(X,Y,density)
% Compute the mean of a bivariate distribution, given the density
% X and Y meshgrid coordinates for the density array
% mu is [E[X], E[Y]]

% 4/2011 bst wrote it

t = sum(density(:));

mu(1,1) = sum(sum(X.*density))/t;
mu(2,1) = sum(sum(Y.*density))/t;
