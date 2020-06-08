function res = kcovario(x)
%% calculates nadwat(x,covario,rsort,50)

  global covario rsort rmax;
  nx = size(x,1);ny=size(x,2);
  T1 = 2000;
  T2 = rmax;
  res=x(:);
  indsmall = x<=T1;
  indlarge = x>T1 & x<=T2;
  indout = x>T2;
  xsmall = x(indsmall);
  xlarge = x(indlarge);
  res(indsmall) = nadwatloop(xsmall,covario,rsort,50);
  resT1 = nadwatloop(T1,covario,rsort,50);
  res(indlarge) = resT1 * (rmax - xlarge) / (rmax-T1);
  res(indout) = 0;
  res = reshape(res,nx,ny);