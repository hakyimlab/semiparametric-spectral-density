function y=ener(theta)
%% y=ener(theta,data,dista,M)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,dfstart,bcoef] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt knots Hpol t

  %------------------------------------
  %% values of param
  %------------------------------------
  smoo = theta(1);
%% FIX this with approx for integer smoothness
%  if(smoo-floor(smoo)<0.01)
%    smoo = floor(smoo)+.01;
%  end
%  if(ceil(smoo)-smoo<0.01)
%    smoo = ceil(smoo)-.01;
%  end
   
  gm = 2*(smoo+1);
%  irango = 2*sqrt(smoo)/rango;
  sig2 = theta(2);
  bcoef = theta(3:end);
  ll = length(bcoef) - 3;

  covacnots = bhankelx(cnots,smoo,bcoef); % bhankel forces covacnots(0)=1

  % define matrix cova
  cova = zeros(size(dista));

  % get spline coeff
  cc = spline(cnots,covacnots);
  % interpolate values of cova(rr) using spline coeff
  cova(indupper) = ppval(cc,rr);

  % complete diag and lower triangular part of cova matrix
  cova = cova + cova' + eye(size(cova));
  cova = cova * sig2;
  
  %%DEBUG
%  disp('cond number');disp(cond(cova))
%  plot(cnots,covacnots)

  
if( min(eig(cova)) > 0 )
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
  y = 1e300;
end
