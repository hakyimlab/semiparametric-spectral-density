function cova=calcCova(theta)
%% y=calcCova(theta)

global Za dista M cnots indupper rr nt 

  %------------------------------------
  %% parse values of param
  %------------------------------------
  smoo = theta(1);
  sig2 = theta(2);
  bcoef = theta(3:end-1);
  wt = theta(end);

  covacnots = bhankelx(cnots,smoo,bcoef,wt); % bhankel forces covacnots(0)=1

  % define matrix cova
  cova = zeros(size(dista));

  % get spline coeff
  cc = spline(cnots,covacnots);
  % interpolate values of cova(rr) using spline coeff
  cova(indupper) = ppval(cc,rr);

  % complete diag and lower triangular part of cova matrix
  cova = cova + cova' + eye(size(cova));
  cova = cova * sig2;
