function res=natwat(r,covaij,rij,h)
%  res=natwat(r,covaij,rij,h)
  nr = numel(r);
  npairs = numel(covaij);
  res = r; res(:)=0;
  
  for ii = 1:nr
    res(ii) = sum( covaij .* exp(-(rij-r(ii)).^2 / 2 / h^2 ) ) / ...
              sum( exp(-(rij-r(ii)).^2 / 2 / h^2 ) );
  end
  
  