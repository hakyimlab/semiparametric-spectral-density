function y=ener(theta)
%% y=ener(theta)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,bcoef,wt] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt reml

cova = calcCova(theta);
mineig = min(eig(cova));
nobs = size(Za,1);


if( mineig > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  nobs = length(Za(:,1));
  if reml
    invcova = invcc * invcc';
    W = M' * invcova * M;
    ww = chol(W);
    invW = inv(W);
    y = nt * (sum(log(diag(cc))) + sum(log(diag(ww))));
    SS = invcova * (eye(nobs) - M * invW * M' * invcova);
    for ii = 1:nt
      z = Za(:,ii);
      y = y + 0.5 * z' * SS * z;
    end
  else
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  end
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


