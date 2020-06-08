
n0=256; %n0=power(2,15);
m=2;
pp=2;
t=2;
del=t/(n0+0.0);
del1=ceil(log(pi/del)/log(2.0));
n= 4^del1; %(int) power(4, (int) del1);
if (n<n0) n=n0; end
mn=m*n;
rt = 5;
phi = rt/(mn-1)/del; 
phi=1;

xx = (0:(mn-1))*del*phi;  
%aa = kcovario(xx);

aa = exp(-xx);
aa(n0+1:end)=0;

[fw,del2]=knhankelc(n0, n,  mn, pp, del, aa);

ww = del2*(0:(mn-1))/phi;
fw = fw*phi^2;

plot(ww(1:n0),fw(1:n0))

