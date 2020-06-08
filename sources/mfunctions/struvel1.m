function sl1=struvel1(x,sl1);
% sl1=stvl1(x,sl1)
% calculates struve function 
pi=3.141592653589793d0;
r=1.0d0;
if (x <= 20.0d0) ;
s=0.0d0;
for  k=1:60;
r=r.*x.*x./(4.0d0.*k.*k-1.0d0);
s=s+r;
if (abs(r) < abs(s).*1.0d-12) break; end;
end;
sl1=2.0d0./pi.*s;
else;
s=1.0d0;
km=fix(.50.*x);
if (x > 50) km=25; end;
for  k=1:km;
r=r.*(2.0d0.*k+3.0d0).*(2.0d0.*k+1.0d0)./(x.*x);
s=s+r;
if (abs(r./s) < 1.0d-12) break; end;
end;
sl1=2.0d0./pi.*(-1.0d0+1.0d0./(x.*x)+3.0d0.*s./x.^4);
a1=exp(x)./sqrt(2.0d0.*pi.*x);
r=1.0d0;
bi1=1.0d0;
for  k=1:16;
r=-0.125d0.*r.*(4.0d0-(2.0d0.*k-1.0d0).^2)./(k.*x);
bi1=bi1+r;
if (abs(r./bi1) < 1.0d-12) break; end;
end;
sl1=sl1+a1.*bi1;
end;
return;

