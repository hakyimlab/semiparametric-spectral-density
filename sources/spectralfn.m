function y = spectralfn(w,smoo,t,bcoef)
% arbitrary spectral density function
  y = w;
%  bcoef = [.2 1 .2 .8 .4 .2];
  sp = spmak(t,bcoef);
  wt = t(end-3);
  
  indismall = find(w<=wt);
  if length(indismall)>0
    y(indismall) = spval(sp,w(indismall));
  end
  
  indilarge = find(w>wt);
  if length(indilarge)>0  
    ft = spval(sp,wt);
    gm = smoo*2 + 2;
    consta = ft*wt^(gm);
    y(indilarge) = consta./w(indilarge).^(gm);
  end
  
  y = y/2/pi;