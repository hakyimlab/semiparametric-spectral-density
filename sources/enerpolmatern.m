function y=ener(theta)
%% y=ener(theta)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,bcoef,wt] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt 

smoo = exp(theta(1));
sig2 = exp(theta(2));
rango = exp(theta(3));
a = 2*sqrt(smoo)/rango;
u = exp(theta(4));
v = exp(theta(5));

cova = polmatern(dista,u,v,a,smoo)*sig2;
mineig = min(eig(cova));

if( mineig > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  xx = invcc' * Za;
  nobs = size(Za,1);
  xx=reshape(xx,nt*nobs,1);
  y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  %%==========================================
else
  disp('WARNING: NEGATIVE EIGENVALUE')
  disp(mineig)
  disp('condition number = ')
  disp(cond(cova));
  disp('theta was ')
  disp(theta)
  y = 1e100;
end


