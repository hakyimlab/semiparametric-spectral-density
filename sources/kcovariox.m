function res = kcovariox(x)
%% calculates nadwat(x,covario,rsort,50).*x
  global covario rsort rmax;
  nx = size(x,1);ny=size(x,2);
  T1 = 1500;
  T2 = rmax;
  res=x(:);
  indsmall = x<=T1;
  indlarge = x>T1 & x<=T2;
  indout = x>T2;
  xsmall = x(indsmall);
  xlarge = x(indlarge);
  res(indsmall) = nadwat(xsmall,covario,rsort,50);
  resT1 = nadwat(T1,covario,rsort,50);
  res(indlarge) = resT1 * (rmax - xlarge) / (rmax-T1);
  res(indout) = 0;
  res = res .* x(:);
  res = reshape(res,nx,ny);