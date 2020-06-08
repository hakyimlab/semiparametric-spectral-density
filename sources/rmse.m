function res=rmse(x,y)
  if nargin==1
    y=x;
  end
  n=length(x);
  res=sqrt( sum( (x(:)-y(:)).^2 )/n );
