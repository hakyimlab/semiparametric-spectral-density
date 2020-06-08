%% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values
global covario rsort rmax;

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
%save(strcat('tablahpol_',num2str(ll)),'Hpol')

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
%testsplinetail = false; testmatern = true; testpolmatern = false;
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
%%randn('state',0);
%------------------------------------
% generate observations
%------------------------------------
indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end);
simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
simsig2=1;
simll = 4;

if testsplinetail
  simt = nodos(simwt,simll);
  simknots = simt(4:end-3);
  %  simbcoef = [0.8 1.0 .8 .2  1.2 1 .8 .6 .4]; %% ll=6 should have ll+3 elements
  simbcoef = [0.8 1 .8  .2  .8 .6 .2]; %% ll=4
  simbcoef(end)=contderiv(simbcoef,simgm);
  simtheta = [simsmoo simsig2 simbcoef simwt];
  load(strcat('tablahpol_',num2str(simll)))
  simcova = calcCova(simtheta);
elseif testpolmatern 
  % test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
  u = 0.5*simirango; v=.1*simirango;
  simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);
elseif testmatern
  % test matern cova
  simcova = simsig2*hmatern(dista,simrango,simsmoo);
end

for ii=1:nt
  Za(:,ii) = simulachol(simcova);
end

if testsplinetail
  save(strcat('splinetail_',num2str(nt),'_',num2str(simll),'.mat'),...
       'Za','simbcoef','simll','simsig2','simsmoo','simgm','simwt',...
       'simt','simknots', 'simtheta','simcova')
elseif testpolmatern 
  save(strcat('polmatern_',num2str(nt),'_','.mat'),...
       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango','u','v',...
       'simcova')
elseif testmatern
  save(strcat('matern_',num2str(nt),'.mat'),...
       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango',...
       'simcova')
end