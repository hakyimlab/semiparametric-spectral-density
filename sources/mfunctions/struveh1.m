function sh1=struveh1(x,sh1);
% sh1=struveh1(x,sh1)
% calculates struve function of order 1

%% HKI corrected bug where if struveh0(0:22)
%% all computation was done as if the whole
%% vector was > 20
  
pi=3.141592653589793d0;

% saving format information
[nrow ncol] = size(x);
x = reshape(x, nrow*ncol, 1);

r=1.0d0;
xsmall = x(x<=20.0d0);
xlarge = x(x>20.0d0);
sh1small=[];
sh1large=[];

if (length(xsmall)>0) ;
  x = xsmall;
  s=0.0d0;
  a0=-2.0d0./pi;
  for  k=1:60;
    r=-r.*x.*x./(4.0d0.*k.*k-1.0d0);
    s=s+r;
    if (abs(r) < abs(s).*1.0d-12) break; end;
  end;
  sh1small=a0.*s;
end

r=1.0d0;
if (length(xlarge)>0)
  x = xlarge;
  s=1.0d0;
  km=fix(.5.*x);
  if (x > 50.d0) km=25; end;
  for  k=1:km;
    r=-r.*(4.0d0.*k.*k-1.0d0)./(x.*x);
    s=s+r;
    if (abs(r) < abs(s).*1.0d-12) break; end;
  end;
  t=4.0d0./x;
  t2=t.*t;
  p1=((((.42414d-5.*t2-.20092d-4).*t2+.580759d-4).*t2-.223203d-3).*t2+.29218256d-2).*t2+.3989422819d0;
  q1=t.*(((((-.36594d-5.*t2+.1622d-4).*t2-.398708d-4).*t2+.1064741d-3).*t2-.63904d-3).*t2+.0374008364d0);
  ta1=x-.75d0.*pi;
  by1=2.0d0./sqrt(x).*(p1.*sin(ta1)+q1.*cos(ta1));
  sh1large=2.0./pi.*(1.0d0+s./(x.*x))+by1;
end;

sh1 = reshape([sh1small; sh1large],nrow,ncol);

return;

