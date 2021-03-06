function sh0=struveh0(x,sh0);
% [x,sh0]=struveh0(x,sh0)
% computes struve function of order 0
%
%% HKI corrected bug where if struveh0(0:22)
%% all computation was done as if the whole
%% vector was > 20
  
pi=3.141592653589793d0;


%% saving format informatin
[nrow ncol] = size(x);
x = reshape(x, nrow*ncol, 1);
sh0small=[];
sh0large=[];

s=1.0d0;
r=1.0d0;
xsmall = x(x<=20.0d0);
xlarge = x(x>20.0d0);
if (length(xsmall)>0) 
  x = xsmall;
  a0=2.0.*x./pi;
  for  k=1:60;
    r=-r.*x./(2.0d0.*k+1.0d0).*x./(2.0d0.*k+1.0d0);
    s=s+r;
    if (abs(r) < abs(s).*1.0d-12) break; end;
  end;
  sh0small=a0.*s;
end
s=1.0d0;
r=1.0d0;
if(length(xlarge>0))
  x=xlarge;
  km=fix(.5.*(x+1.0));
  if (x >= 50.0) km=25; end;
  for  k=1:km;
    r=-r.*((2.0d0.*k-1.0d0)./x).^2;
    s=s+r;
    if (abs(r) < abs(s).*1.0d-12) break; end;
  end;
  t=4.0d0./x;
  t2=t.*t;
  p0=((((-.37043d-5.*t2+.173565d-4).*t2-.487613d-4).*t2+.17343d-3).*t2-.1753062d-2).*t2+.3989422793d0;
  q0=t.*(((((.32312d-5.*t2-.142078d-4).*t2+.342468d-4).*t2-.869791d-4).*t2+.4564324d-3).*t2-.0124669441d0);
  ta0=x-.25d0.*pi;
  by0=2.0d0./sqrt(x).*(p0.*sin(ta0)+q0.*cos(ta0));
  sh0large=2.0d0./(pi.*x).*s+by0;
end;
sh0=reshape([sh0small; sh0large],nrow,ncol);

return;

