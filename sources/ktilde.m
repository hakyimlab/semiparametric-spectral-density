%% calc positive definite kernel covariogram at r
%% by hankel tr, taking positive part and hankel 
%% transforming back

function [res,pp] = ktilde(r,R)
  [w,fw]=sieghankel(@kcovario,R);
  pp = spline(w,fw);
  [x,kx] = sieghankel(@fdens,1/R,pp);
  pp = spline(x,kx);
  res = ppval(pp,r);
  
  
function res = fdens(x,pp) 
  %% returns truncated spectral density
  %% corresp to kcovario
  res = ppval(pp,x);
  res(find(res<0))=0;
    