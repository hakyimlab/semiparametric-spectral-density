function res=simulachol(sigma)
%------------------------------------
% simulates a gaussian random field using chol decomposition
% mean 0
n=size(sigma,1);
e=randn(n,1);
C=chol(sigma);
res=C'*e;
%------------------------------------

