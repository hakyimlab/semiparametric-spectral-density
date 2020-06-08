function sl0=struvel0(x,sl0);
% sl0=stvl0(x,sl0)
% calculates struve function of order 1
pi=3.141592653589793d0;
s=1.0d0;
r=1.0d0;
if (x <= 20.0d0) ;
a0=2.0d0.*x./pi;
for  k=1:60;
r=r.*(x./(2.0d0.*k+1.0d0)).^2;
s=s+r;
if (abs(r./s) < 1.0d-12) break; end;
end;
sl0=a0.*s;
else;
km=fix(.5.*(x+1.0));
if (x >= 50.0) km=25; end;
for  k=1:km;
r=r.*((2.0d0.*k-1.0d0)./x).^2;
s=s+r;
if (abs(r./s) < 1.0d-12) break; end;
end;
a1=exp(x)./sqrt(2.0d0.*pi.*x);
r=1.0d0;
bi0=1.0d0;
for  k=1:16;
r=0.125d0.*r.*(2.0d0.*k-1.0d0).^2./(k.*x);
bi0=bi0+r;
if (abs(r./bi0) < 1.0d-12) break; end;
end;
bi0=a1.*bi0;
sl0=-2.0d0./(pi.*x).*s+bi0;
end;
return;

