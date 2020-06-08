%% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt
global smoovect wtvect % tabulated values
global covario rsort rmax;


[lat,wnh4,lon,cwnh4]=textread('simul.dat','%f%f%f%f%*[^\n]');
%M = load('M.txt');
nobs = length(wnh4);
M = ones(nobs,1);
dista = cordist([lon,lat]);
indiag=find(eye(nobs)); dista(indiag)=0; %% clean up the diag
indupper = find(triu(dista,1)); 
rr = dista(indupper); rmin = min(rr); rmax = max(rr);
nnots=100; cnots = linspace(rmin,rmax,nnots);

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

% generate observations

simsmoo = 3;
indi = find(smoovect>=simsmoo); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simwt = 0.0100; %% simwt=0.0038
simwt = wtvect(find(wtvect <= simwt)); simwt = simwt(end);
simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
simsig2=1;
simll = 4;

% test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
u = 0.5*simirango; v=.1*simirango;
simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);

for ii=1:nt
  Za(:,ii) = simulachol(simcova);
end  

%% calc covariogram cloud and nadwat estimate of cov fn

dindupper = find(triu(ones(size(dista))));

ncov = length(dindupper);
covario = nan(ncov,nt);
[rsort,indsort] = sort(dista(dindupper));

zm = 0;

for ii=1:nt
  z = Za(:,ii); 
  zdif = repmat(z,1,nobs);
  zdif = (zdif -zm).*(zdif -zm)';
  zdiflin = zdif(dindupper);
  covario(:,ii) = zdiflin(indsort);
end
  
covario = sum(covario,2)/nt;

rohat = kcovario(cnots);

tic
%% hankel transform of rohat
wtf = 1/50;
xxx = linspace(0,2*wtf,100);
fcovario = xxx(2:end); 
for ii = 2:length(xxx)
  fcovario(ii-1) = myhankel('kcovario',xxx(ii));
end
toc

R = 500;

tic
[ww,gg] = sieghankel(@kcovario, R);
pp = spline(ww,gg);
fcovario2 = ppval(pp,xxx);
toc

plot(ww,gg);hold on;plot(xxx(2:end),fcovario,'r.','markersize',5);hold off

tic
  R = 300;
%% get pos def covariogram
kk = ktilde(cnots,R);
plot(cnots,kk)
hold on; plot(cnots,kcovario(cnots),'r.');hold off
toc


