%% test bhankelx(cnots,gm,bcoef)~hmatern(cnots,rango,smoo)
clear all
global Za dista M cnots indupper rr nt knots Hpol t tailmat smoovect
ngm = 50;
smoovect = linspace(.05,5,ngm)+.0001;
smoo=2.0001;
gm=2*smoo+2;
rango=50;
irango=2*sqrt(smoo)/rango;
wt=2*irango;
cnots = (.1:.1:4)*rango;
%% hankel transform of polynomials
Hpol = calcHpol(knots,cnots);
%% hankel tranform of tail
tailmat = calcTailmat(cnots*wt);

ll=100;
t = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*wt;
knots = (0:ll)/ll*wt;
bcoef = fmatern(abs(t(3:end-2)),irango,smoo);
sp = spmak(t,bcoef)
xx = linspace(0,wt,20);
hold off
plot(xx,spval(sp,xx),'r.')
hold on
plot(xx,fmatern(xx,irango,smoo))
hold off

Hpol = calcHpol(knots,cnots);
trucova = hmatern(cnots,rango,smoo);
mycova = bhankelx(cnots,smoo,bcoef);

plot(cnots(1:end),trucova(1:end))
hold on
plot(cnots(1:end),mycova(1:end),'r.')
hold off

rmse(trucova,mycova)

[thetaff,fval]=fminunc(ener,thetaf)

hold off
sp = spmak(t,thetaf(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmf./xx.^gmf,'r.')
grid
print -depsc plotf.eps

hold off
sp = spmak(t,thetaff(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmff./xx.^gmff,'r.')
grid
print -depsc plotff.eps

hold off
sp = spmak(t,thetafff(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmfff./xx.^gmfff,'r.')
grid
print -depsc plotfff.eps



kksmoof = 2.4251;
kksig2f = 1.0129;
kkbcoeff=[ 
    1.5273
    8.1425
    1.5273
    0.2079
    0.0135
    0.2214
];

hold off
sp = spmak(t,bcoeff);
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gm./xx.^gm,'r.')




ener([smoof sig2f [bcoeff(1:1) 100 bcoeff(3:end)]])



smoof = 2.5216;
sig2f = 1.1226;
kkbcoeff=[    0.4566
    0.3851
    0.4566
    0.0581
  100.4991
    0.6671
    0.7671
   20.4019
    7.9151
    6.6429
    0.3399
    0.8463
    3.6148
    0.8002
    0.7406
    0.1371
    0.6588
    0.2072];





xx =linspace(1,100,10) ;
yy = linspace(1,100,10);
[xxx,yyy]=meshgrid(xx,yy);
for ii=1:length(xxx);
   for jj=1:length(yyy);
     zlik(ii,jj)=-ener([smoof sig2f [bcoeff(1:3) xxx(ii,jj) yyy(ii,jj) bcoeff(6:end)]]);
   end;
 end
meshc(xxx,yyy,zlik)
xlabel('f(4/3 1e-3'),ylabel('f(2e-3)'),zlabel('likeli')


rango=2000;smoo=6.01;gm=2*(smoo+1);sig2=3;
irango = 2*sqrt(smoo)/rango;
knots = (0:ll)/ll*irango*alpha; %% interval of interest (0:wt)
t= [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*irango*alpha; %% deBoor knots
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); %% coef of bsplines approx
bcoef(end) = contderiv(bcoef,gm);  %% continuity of deriv at wt
cnots = 



%% TEST: bhankel ~= true covariance
rango=500;smoo=2.01;gm=2*(smoo+1);
ll=200;cnots = rango*(.1:.1:2);
irango = 2*sqrt(smoo)/rango;
wt=3*irango;
knots = (0:ll)/ll*wt; %% interval of interest (0:wt)
t= [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*wt; %% deBoor knots
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); %% coef of bsplines approx
bcoef(end) = contderiv(bcoef,gm);  
figure(1)
xx = linspace(0,wt,20);
sp=spmak(t,bcoef); 
hold off
plot(xx,spval(sp,xx),'r.')
hold on
plot(xx,fmatern(xx,irango,smoo))
hold off

Hpol = calcHpol(knots,cnots);

trucova=hmatern(cnots,rango,smoo);
mycova=bhankelx(cnots,gm,bcoef);
disp(rmse(trucova,mycova))
hold off
plot(cnots,trucova,'b-')
hold on
plot(cnots,mycova,'.r')

print -depsc plot.eps



%% TEST: profile likelihood
simrango=2000;simsmoo=6;simsig2=1;simgm=2*(simsmoo+1);
z = simulachol(hmatern(dista,simrango,simsmoo)*simsig2);

xx=.2:.4:6;
pll = zeros(size(xx));
for ii=1:length(xx)
  pll(ii)=enersmoo(xx(ii),simrango,sig2,z,dista)
end
plot(xx,pll,'.')
pll
hold on; plot([0 6],(min(pll)+2)*[1 1]);hold off





indupper = find(triu(dista,1));
rr=dista(indupper);
kkcova = zeros(size(dista));kksimcova=kkcova;

rmin=min(rr); rmax=max(rr);
kkcnots = linspace(rmin,rmax,100);
alpha = 5;

simrango = 1000;
simsmoo = 1.1; simgm = (simsmoo+1)*2;
simirango = 2*sqrt(simsmoo)/simrango;
simknots = (0:ll)*simirango/alpha;
simt = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]*simirango/alpha;
simbcoef = fmatern(abs(simt(3:(end-2))),simirango,simsmoo); 
simbcoef(end)=contderiv(simbcoef,simgm);
kksimcovacnots = bhankel(kkcnots,simgm,simknots,simbcoef);

rango = 1000;
smoo = 5.1; gm = (smoo+1)*2;
irango = 2*sqrt(smoo)/rango;
knots = (0:ll)*irango/alpha; 
t = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]*irango/alpha;
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); 
bcoef(end)=contderiv(bcoef,gm);
kkcovacnots = bhankel(kkcnots,gm,knots,bcoef);

cc = spline(kkcnots,kksimcovacnots);
kksimcova(indupper)=ppval(cc,rr);
cc = spline(kkcnots,kkcovacnots);
kkcova(indupper)=ppval(cc,rr);

kksimcova = kksimcova + kksimcova' + eye(size(kkcova));
kkcova = kkcova + kkcova' + eye(size(kkcova));

cond(kksimcova)
cond(kkcova)

plot(kkcnots,kksimcovacnots,'b.');hold on
plot(kkcnots,kkcovacnots,'r.');hold off




xx=t(3:end-2);
plot(xx,simbcoef,'b.');hold on
plot(xx,bcoef,'r.');hold off








%%here
rango = 500;
gm=3.1;smoo =gm/2 - 1;
cnots = (.1:.1:5)*rango;
protknots=0:20;
knots = protknots/rango/5;
lend = protknots(end)
t = [-3 -2 -1 protknots lend+1 lend+2 lend+3]/rango/5;
xint = t(3:(end-2));
bcoef = fmatern(xint,1/rango,smoo)*2*pi;
covacnots = bhankel(cnots,gm,knots,bcoef);
notailcovacnots = kkbhankel(cnots,gm,knots,bcoef)
plot(cnots,matern(cnots,rango,smoo));hold on
plot(cnots,covacnots,'r.');
plot(cnots,notailcovacnots,'g.');hold off















tic
for ii=1:NN
  st = cnots(ii)/rango*2;
  covatail(ii)=hg1f2(1-gm/2,1,2-gm/2,-st^2/4)*st^(2-gm)/(gm-2) ...
               - gm*gamma(-gm/2)/gamma(gm/2)/2^gm ;
  covatailapp(ii)=kktail(st,gm);
end
toc

plot(cnots(100:end)*2/rango,covatail(100:end),'.');hold on;grid
plot(cnots(100:end)*2/rango,covatailapp(100:end),'r.');hold off


gm=15.1;
knots=nodos(1/30);
fknots=repmat(1,1,length(cnots));
kkcova=bhankel(cnots,gm,knots,fknots);
plot(cnots*2/rango,kkcova)

%=============================================
l=10;
knots = 0:10;
xx= [-3 -2 -1 knots 11 12 13];
xint = xx(3:(end-2));
n = length(xint)
%yy=(xint-3).*(xint-6).*(xint-9);
%yy = repmat(0,1,n);yy(6)=1;
yy= xint.^3;
kk=mexbsplpp(xx,yy)

ll = n - 3;
kk = kk./repmat(cumprod([1 1:3]'),1,ll)
kk=kk(4:-1:1,:)';
pp=mkpp(xx(4:n+1),kk);


plot(xint,yy,'r');
xtest=xint(1)+(1:100)/100*(xint(end)-xint(1));
hold on;
plot(xtest,ppval(pp,xtest),'b-')
hold off
