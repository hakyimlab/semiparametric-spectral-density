%%function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values

%------------------------------------
%% read sites locations
%% simul.dat: "lat"   "wnh4"  "lon"   "cwnh4" "site"
%------------------------------------
[lat,wnh4,lon,cwnh4]=textread('simul.dat','%f%f%f%f%*[^\n]');
M = load('M.txt');
nobs = length(wnh4);
dista = cordist([lon,lat]);
indupper = find(triu(dista,1)); 
rr = dista(indupper); rmin = min(rr); rmax = max(rr);
nnots=100; cnots = linspace(rmin,rmax,nnots);

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

ll = 6; %number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

%% calculate or read tabulated values
%% hankel transform of polynomials
%tic
%Hpol = calcHpol(ll,cnots,wtvect);
%disp('hPol calc time')
%toc
%save strcat('tablahpol_',num2str(ll)) Hpol

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


%------------------------------------
%% read observations or simulate
%------------------------------------
simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end)
simt = nodos(simwt,ll);
simknots = simt(4:end-3);
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simsig2=1;
simbcoef = [0.2 1.0 0.2 0.8 0.4 0.2];
simbcoef(end)=contderiv(simbcoef,simgm)
simtheta = [simsmoo simsig2 simbcoef simwt];
randn('state',0);
%simcova = calcCova(simtheta);

%% test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
simrango = 500; simirango = 2*sqrt(simsmoo)/simrango;
u = 0.5*simirango; v=.1*simirango;
simcova = polmatern(dista,u,v,simirango,simsmoo);
simwt = 1.5*simirango;

%% test matern cova
%simrango = 500; simirango = 2*sqrt(simsmoo)/simrango;
%simcova = hmatern(dista,simrango,simsmoo);

for ii=1:nt
  Za(:,ii) = simulachol(simcova*simsig2);
end

cc = chol(simcova);
invcc = inv(cc);
xx = invcc' * Za;
nobs = size(Za,1);
xx=reshape(xx,nt*nobs,1);
truene = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
disp(truene)

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=0.2); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=3;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= simwt*1); wt = wt(end)

t = nodos(wt,ll);
theta = [smoo,sig2,bcoef,wt];thetaini=theta;
Ei = ener(theta); Eini = Ei;


%------------------------------------
%% improve initial guess fitting matern
%------------------------------------
matheta = log([smoo sig2 500]);
%[paramf fval exitflag] = fminunc(@enermatern,matheta)
[paramf fval exitflag] = fminsearch(@enermatern,matheta);
matheta = exp(paramf);
matsmoo = smoovect(smoovect<=matheta(1));matsmoo=matsmoo(end);
matsig2 = matheta(2);
matirango = 2*sqrt(matsmoo)/matheta(3);
matwt = wtvect(wtvect <= matirango*.75); matwt = matwt(end);

theta = [matsmoo matsig2 bcoef matwt ];
matEi = ener(theta);

%------------------------------------
%% loop for simulated annealing
%------------------------------------
%% minimize (-l(f)+penalty)
L = 20000; %% max number of iterations
Ti = 1; %% scale temperature
lam = 0.995;
tsca = Ti * lam.^(1:L);

fixwt = false;
conta=1;
llvect = nan(L,1);
for T=tsca
  % generate new state %logf m-vector
  thetaf = bproposal(theta,fixwt);  %% bproposal.m generates f and proposal.m log(f)
  %% evaluate energy 
  Ef = ener(thetaf);
  %disp(Ef)
  % accept or reject
  dE = Ef-Ei;
  if dE < 0
    disp(conta)
    disp('descent')
    theta = thetaf;Ei = Ef;%disp(theta)
    llvect(conta) = - Ef;
  else
    p = exp(-dE/T)/(1+exp(-dE/T));
    pu = rand(1);
    if pu<p
      theta = thetaf; Ei = Ef; %disp(theta)
      llvect(conta)=-Ef;
      disp('ascent')
    end
  end
  conta=conta+1;
end
thetaf = theta;
Ef = Ei;

%save matlaboutput simtheta truene thetaf Ef;

smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll)

disp('number of indep replicates')
disp(nt)
disp('sim  ini  after-matern  final')
disp([simsmoo smoo matsmoo smoof ])
disp([simsig2 sig2 matsig2 sig2f])
disp([simwt wt matwt wtf])
disp(round([truene Eini matEi Ei]))
disp('estimated ener - true ener'); disp(Ei - truene)
disp('iter');disp(L)
disp('lam -- cooling schedule T=lam^(1:iter))');disp(lam)

figure(1)
simf0 = spectralfn(simwt,simsmoo,simt,simbcoef);
ff0=spectralfn(simwt,smoof,tf,bcoeff);
xx = linspace(0,2*simwt,100);
plot(xx,spectralfn(xx,simsmoo,simt,simbcoef)/simf0);hold on
plot(xx,spectralfn(xx,smoof,tf,bcoeff)/ff0,'r.');hold off
grid
xlabel('freq')
ylabel('f(w)')
title('freq')
legend('simulated','estimated')

figure(2)
covaf = calcCova(theta);
%simcova = hmatern(dista,simrango,simsmoo);
%plot(covaf(:),simcova(:),'.')
%hold on; plot([0 max(covaf(:))],[0 max(covaf(:))]);hold off
eigf = eig(covaf);simeig=eig(simcova);
plot(eig(covaf),eig(simcova),'.')
hold on; plot([0 max(eigf)],[0 max(eigf)]);hold off
xlabel('eig(cova) estimated')
ylabel('eig(cova) simulated')
title('eig of cova -- simulated vs. estiamted')
grid

figure(3);
%length( llvect( find(~isnan(llvect)) ) );
%evol = cumsum(~isnan(llvect));
%plot(evol); ylabel('cumulated number of changes in ener')
%xlabel('iteration number');title('evolution of energy')
[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);
plot(cnots,polmatern(cnots,u,v,simirango,simsmoo))
hold on; plot(cnots,covacnots,'r.')
grid
hold off
xlabel('distance');title('cov function')
legend('simulated','estimated')

figure(1);print -depsc plotspectral.eps
figure(2);print -depsc ploteig.eps
figure(3);print -depsc plotcovafn.eps

%[thetaff,fval,exitflag]=fminunc(@ener,thetaf);

%% TEST MATERN or polmatern
figure(1)
%simf0 = fmatern(0,simirango,simsmoo);
%simf0 = fpolmatern(0,u,v,simirango,simsmoo);
%ff0=spectralfn(0,smoof,tf,bcoeff);
xx = linspace(0,2*simwt,100);
%plot(xx,fmatern(xx,simirango,simsmoo)/simf0);hold on
plot(xx,fpolmatern(xx,u,v,simirango,simsmoo));hold on
plot(xx,spectralfn(xx,smoof,tf,bcoeff)/cova0,'r.');hold off
grid
xlabel('freq')
ylabel('f(w)')
title('freq')
legend('simulated','estimated')

figure(1);print -depsc plotspectral.eps
figure(2);print -depsc ploteig.eps
figure(3);print -depsc plotcovafn.eps


poltheta = log([smoof sig2f 500 .1 .1]);
[kkf fval exitflag] = fminsearch(@enerpolmatern,poltheta)
uf = exp(kkf(4));
vf=exp(kkf(end));
smf = exp(kkf(1));
af = 2*sqrt(smf)/exp(kkf(3));
figure(1);hold on;
plot(xx,fpolmatern(xx,uf,vf,af,smf),'k.');hold off