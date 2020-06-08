%% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values
global covario rsort rmax;

%------------------------------------
%% read sites locations
%% simul.dat: "lat"   "wnh4"  "lon"   "cwnh4" "site"
%------------------------------------
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

ll = 3 % number of polynomial pieces
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


%testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
testsplinetail = false; testmatern = true; testpolmatern = false;
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
%%randn('state',0);
%------------------------------------
% generate observations
%------------------------------------
%indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
%simgm=2*(simsmoo+1);
%simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end);
%simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
%simsig2=1;
%simll = 4;

%if testsplinetail
%  simt = nodos(simwt,simll);
%  simknots = simt(4:end-3);
%  %  simbcoef = [0.8 1.0 .8 .2  1.2 1 .8 .6 .4]; %% ll=6 should have ll+3 elements
%  simbcoef = [0.2 1 .2  2  .6 .4 .2]; %% ll=4
%  simbcoef(end)=contderiv(simbcoef,simgm);
%  simtheta = [simsmoo simsig2 simbcoef simwt];
%  simcova = calcCova(simtheta);
%elseif testpolmatern 
%  % test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
%  u = 0.5*simirango; v=.1*simirango;
%  simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);
%elseif testmatern
%  % test matern cova
%  simcova = simsig2*hmatern(dista,simrango,simsmoo);
%end

%for ii=1:nt
%  Za(:,ii) = simulachol(simcova);
%end

%if testsplinetail
%  save(strcat('splinetail_',num2str(nt),'_',num2str(ll),'.mat'),...
%       'Za','simbcoef','simll','simsig2','simsmoo','simgm','simwt',...
%       'simt','simknots', 'simtheta','simcova')
%elseif testpolmatern 
%  save(strcat('polmatern_',num2str(nt),'_','.mat'),...
%       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango','u','v',...
%       'simcova')
%elseif testmatern
%  save(strcat('matern_',num2str(nt),'.mat'),...
%       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango',...
%       'simcova')
%end

%% read observations
if testsplinetail
  simll = 4;
  load(strcat('splinetail_',num2str(nt),'_',num2str(simll),'.mat'));
elseif testpolmatern 
  load(strcat('polmatern_',num2str(nt),'_','.mat'))
elseif testmatern
  load(strcat('matern_',num2str(nt),'.mat'))
end

cc = chol(simcova);
invcc = inv(cc);
xx = invcc' * Za;
nobs = size(Za,1);
xx=reshape(xx,nt*nobs,1);
truene = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
clear('xx','cc');

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=0.2); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=3;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= .02); wt = wt(end);

t = nodos(wt,ll);
theta = [smoo,sig2,bcoef,wt];thetaini=theta;
Ei = ener(theta); Eini = Ei;


%------------------------------------
%% improve initial guess fitting matern
%------------------------------------
matheta = log([smoo sig2 1000]);
%[paramf fval exitflag] = fminunc(@enermatern,matheta)
[paramf fval exitflag] = fminsearch(@enermatern,matheta);
matheta = exp(paramf);
matsmoo = smoovect(smoovect<=matheta(1));matsmoo=matsmoo(end);
matsig2 = matheta(2);
matirango = 2*sqrt(matsmoo)/matheta(3);
matwt = wtvect(wtvect <= matirango); matwt = matwt(end);

theta = [matsmoo matsig2 bcoef matwt];
%% DEBUG

testwt = wtvect(wtvect <= 1.5*simwt); testwt = testwt(end);
theta = [matsmoo matsig2 bcoef testwt];


matEi = ener(theta); Ei = matEi;

%------------------------------------
%% loop for simulated annealing
%------------------------------------
%% minimize (-l(f)+penalty)
L = 3000 %%max number of iterations
Ti = 1000; %% scale temperature
lam = 0.99;
%tsca = Ti * lam.^(1:L);


tsca = ones(L,1);
a=300;
for kk=1:L;
  tsca(kk+1) = tsca(kk)/(1+a*tsca(kk));
end
tsca = tsca/tsca(1) * Ti;
exp(-.01/tsca(end))/(1+ exp(-.01/tsca(end)) )
plot(tsca(100:end))

indtheta = [1:2 4:length(theta)];indtheta(end-1)=[]; %% bcoef(1) and
                                                     %bcoef(end) not changed
np = length(indtheta);

fixwt = true;
conta=1;
llvect = nan((L+1)*np+1,1);
llvect(1) = Ei;
thetavect = nan((L+1)*np+1,length(theta));
thetavect(1,:)=theta;
gibbs = false;


for lind = 1:L
  T=tsca(lind);
  % generate new state %logf m-vector
  for i=1:np
    %% proposal theta
    if gibbs
      newtheta = bproposal(theta,fixwt,gibbs);  %% bproposal.m generates f and proposal.m log(f)
      thetaf = theta; thetaf(indtheta(i)) = newtheta(indtheta(i));
      thetaf(3)=thetaf(5);thetaf(end-1)=contderiv(thetaf(3:end-1),thetaf(1)*2+2);
    else
      thetaf = bproposal(theta,fixwt,gibbs);
    end
    %% evaluate energy 
    Ef = ener(thetaf);
    % accept or reject
    dE = Ef-Ei;
    if dE < 0
      theta = thetaf;Ei = Ef;%disp(theta)
      llvect(conta) = - Ef;
    else
      p = exp(-dE/T)/(1+exp(-dE/T));
      pu = rand(1);
      if pu<p
        theta = thetaf; Ei = Ef; %disp(theta)
        llvect(conta)=-Ef;
        %      disp('ascent')
      end
    end
    conta=conta+1;
    if(conta/2000==floor(conta/2000))
      disp(conta)
    end
    
    %% keep track of theta
    thetavect((lind-1)*np+i+1,:)=theta;
    
  end
end
thetaf = theta; %% revert back to last updated values
Ef = Ei; %% revert back to last updated values

smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll);

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


[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);

%% calc covariogram cloud and nadwat estimate of cov fn
dindupper = find(triu(ones(size(dista))));

ncov = length(dindupper);
covario = nan(ncov,nt);
[rsort,indsort] = sort(dista(dindupper));

zm = mean(Za(:));
for ii=1:nt
  z = Za(:,ii); 
  zdif = repmat(z,1,nobs);
  zdif = (zdif -zm).*(zdif -zm)';
  zdiflin = zdif(dindupper);
  covario(:,ii) = zdiflin(indsort);
end
covario = sum(covario,2)/nt;

rohat = kcovario(cnots);

%% hankel transform of rohat
xxx = linspace(0,2*simwt,100);
fcovario = xxx(2:end); 
for ii = 2:length(xxx)
  fcovario(ii-1) = hankel0('kcovariox',xxx(ii));
end

if testsplinetail
  load(strcat('tablahpol_',num2str(simll)))  
  [simcovacnots,simcova0] = bhankelx(cnots,simsmoo,simbcoef,simwt);
  simspectral = spectralfn(xxx,simsmoo,simt,simbcoef)/simcova0;
  prefijo = strcat(num2str(nt),'_',num2str(simll),'_',...
                    num2str(ll));
  epsfile = strcat('splinetail_',prefijo,'.eps');
  savefile = strcat('opt_st_',prefijo,'.mat');
elseif testpolmatern 
  simcovacnots = polmatern(cnots,u,v,simirango,simsmoo);
  simspectral = fpolmatern(xxx,u,v,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_', num2str(ll));
  epsfile = strcat('polmatern_',prefijo,'.eps');
  savefile = strcat('opt_polmat_',prefijo,'.mat');
elseif testmatern
  simcovacnots = hmatern(cnots,simrango,simsmoo);
  simspectral = fmatern(xxx,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_',num2str(ll));
  epsfile = strcat('matern_',prefijo,'.eps');
  savefile = strcat('opt_matern_',prefijo,'.mat');
end

mlespect = spectralfn(xxx,smoof,tf,bcoeff)/cova0;
rmsecovak = rmse(rohat,simcovacnots);
rmsecovam = rmse(covacnots,simcovacnots);
rmsespeck = rmse(fcovario/2/pi,simspectral(2:end));
rmsespecm = rmse(mlespect(2:end),simspectral(2:end));

%% calc likelihood using rohat
kcova = kcovario(dista);
cc = chol(kcova);
invcc = inv(cc);
xx = invcc' * Za;
xx=reshape(xx,nt*nobs,1);
klike = -( nt*sum(log(diag(cc))) + 0.5 * xx'*xx );

%% clear some variables and save
clear('xx','cc');
clear('Hpol','tailmat')
save(savefile);


%% plot results

fighandle = figure(1);
set(fighandle,'PaperPosition',[0.25 2.5 8 3.5]);
subplot(1,2,1)
plot(cnots,simcovacnots)
hold on; plot(cnots,covacnots*sig2f,'r.')
grid
%plot(rsort,covario,'.','markersize',2);
plot(cnots,nadwat(cnots,covario,rsort,50),'k+');hold off

rx =800;ry=0.82;dx=600;dy=-0.12;
text(rx,     ry,'sim')
text(rx+dx,  ry,'mle')
text(rx+2*dx,ry,'kernel')

text(rx-dx/2,ry+dy,'likeli')
text(rx-dx/2,ry+2*dy,'rmse')
text(rx-dx/2,ry+3*dy,'nu')
text(rx-dx/2,ry+4*dy,'sig2')
text(rx-dx/2,ry+5*dy,'wt')

text(rx,ry+dy,sprintf('%4.0f',-truene))
text(rx+dx,ry+dy,sprintf('%4.0f',-Ei))
text(rx+2*dx,ry+dy,sprintf('%4.0f',klike))

text(rx,ry+2*dy,sprintf('%4.3f',0))
text(rx+dx,ry+2*dy,sprintf('%4.3f',rmsecovam))
text(rx+2*dx,ry+2*dy,sprintf('%4.3f',rmsecovak))

text(rx,ry+3*dy,sprintf('%3.2f',simsmoo))
text(rx+dx,ry+3*dy,sprintf('%3.2f',smoof))
text(rx+2*dx,ry+3*dy,'NA')

text(rx,ry+4*dy,sprintf('%3.2f',simsig2))
text(rx+dx,ry+4*dy,sprintf('%3.2f',sig2f))
text(rx+2*dx,ry+4*dy,sprintf('%3.2f',kcovario(0)) )

text(rx,ry+5*dy,sprintf('%5.4f',simwt) )
text(rx+dx,ry+5*dy,sprintf('%5.4f',wtf) )
text(rx+2*dx,ry+5*dy,'NA')

legend('simulated','estimated','kernel estimator')
xlabel('distance');title('cov function')

subplot(1,2,2)
plot(xxx,simspectral);hold on
plot(xxx,mlespect,'r.')
plot(xxx(2:end),fcovario/2/pi,'k+');hold off
legend('simulated','estimated','kernel estimator')
grid
xlabel('freq')
ylabel('f(w)')
title('Spectral density')
ux = 1.4*simwt; uy= max([fcovario/2/pi simspectral mlespect])/2;
dux = .3*simwt;duy=-uy*.2;

text(ux-dux,uy+duy,'rmse')

text(ux,uy,'   mle')
text(ux+dux,uy,'kernel')

text(ux,uy+duy,sprintf('%5.0f',rmsespecm))
text(ux+dux,uy+duy,sprintf('%5.0f',rmsespeck))



disp('rohat')
disp(rmse(rohat,simcovacnots))
disp('mle')
disp(rmse(covacnots,simcovacnots))

print( '-depsc',epsfile)
