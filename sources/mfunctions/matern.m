function y = matern(u,range,smoo)
%------------------------------------
% matern function 
urange = u/range;
y=ones(size(u));
unonzero=find(u);
y(unonzero) = (urange(unonzero).^smoo) .* besselk(smoo,urange(unonzero))/ ...
        (2^(smoo - 1))/ gamma(smoo);
%------------------------------------
