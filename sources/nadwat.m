function res=natwat(r,covaij,rij,h)
%  res=natwat(r,covaij,rij,h)
  nr = numel(r);
  npairs = numel(covaij);
  covaij = repmat(covaij,1,nr);
  rij = repmat(rij,1,nr);
  r=reshape(r,1,nr);
  r = repmat(r,npairs,1);
  res=sum(covaij.*exp(-(rij-r).^2/2/h^2)) ./ sum(exp(-(rij-r).^2/2/h^2));
  