function [pv, pd, covar] = bivar(X,Y,density)
% function [pv, pd, covar] = bivar(X,Y,density)
% Compute the covariance matrix (covar) of a bivariate distribution of a 
% given density.  Also output the principal directions (columns of pd) 
% and variance along the principal directions (pv(1) is the larger one).
% X and Y are meshgrid coordinates of density

% 4/2011 bst wrote it

mu = bimean(X,Y,density);
X = X-mu(1);
Y = Y-mu(2);

t = sum(density(:));

XX = sum(sum(X.*X.*density))/t;
YY = sum(sum(Y.*Y.*density))/t;
XY = sum(sum(X.*Y.*density))/t;

covar = [XX, XY; XY, YY];
[pd, pv] = eig(covar);
pv = diag(pv);
pv = pv([2,1]);
pd = pd(:,[2,1]);