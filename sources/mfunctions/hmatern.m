function y = hmatern(u,range,smoo)
%------------------------------------
% matern function 
  y=ones(size(u));
  unonzero=find(u);
  if smoo<100
    urange = u*2*sqrt(smoo)/range;
    y(unonzero) = (urange(unonzero).^smoo) .* besselk(smoo,urange(unonzero))/ ...
        (2^(smoo - 1))/ gamma(smoo);
  else
    y(unonzero) = exp(-u(unonzero).^2);
  end
  
%------------------------------------
