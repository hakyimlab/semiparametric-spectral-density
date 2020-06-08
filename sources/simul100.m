%% calls nogibbs to get Nsim simulations
%% of process and corresponding estimates
addpath('./mfunctions')

type simul100
type nogibbs

global Za dista M cnots indupper rr nt reml
global Hpol tailmat smoovect wtvect % tabulated values
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

%% interpolation points
%intpoints = [ quantile(lon,repmat([0.25 0.5 0.75],1,3)') ...
%              quantile(lat,reshape(repmat([0.25 0.5 0.75]',1,3)',9,1))];
nmk = 10-1; %% (nmk+1)^2 is number of points to interpolate
ninterp = (nmk+1)^2;
minlon = min(lon);maxlon=max(lon);
minlat = min(lat);maxlat=max(lat);
lonvect = minlon + (0:nmk)'/nmk * (maxlon - minlon);
latvect = minlat + (0:nmk)'/nmk  * (maxlat - minlat);
intpoints = [repmat(lonvect,nmk+1,1) reshape(repmat(latvect',nmk+1,1),(nmk+1)^2,1)];
dist2obs = cordist([lon,lat],[intpoints]);

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

ll = 3 % number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

nsimini=1;Nsim = 200;
%testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
testsplinetail = false; testmatern = true; testpolmatern = false;
nparam = 2+ll+3+1;
fixwt = false;
gibbs = false;
reml = false;

%% simulated annealing parameters
L = 20000; %%max number of iterations
Ti = 1000; %% scale temperature
%lam = 0.99;
%tsca = Ti * lam.^(1:L);
%% cooling schedule
tsca = ones(L,1);
a=30;
for kk=1:L;
  tsca(kk+1) = tsca(kk)/(1+a*tsca(kk));
end
tsca = tsca/tsca(1) * Ti;
%exp(-.01/tsca(end))/(1+ exp(-.01/tsca(end)) )
%plot(tsca(100:end))

if testsplinetail
  modelo = 'splinetail';
elseif testpolmatern 
  modelo = 'polmatern';
elseif testmatern
  modelo = 'matern';
end

if nsimini ==1
  inivect1   = nan(Nsim,nparam);
  inivect50  = nan(Nsim,nparam);
  inivect100 = nan(Nsim,nparam);
  
  estimate1  = nan(Nsim,nparam+8); %% thetaf
  estimate50 = nan(Nsim,nparam+8); %% thetaf
  estimate100= nan(Nsim,nparam+8); %% thetaf
  
  msemlevect = nan(Nsim,ninterp); %% mse true kriged and mle kriged
  msekervect = nan(Nsim,ninterp); %% mse true kriged and ker kriged
  msematvect = nan(Nsim,ninterp); %% mse true kriged and mat kriged
  
  estsimvect = nan(1,ninterp);   %% estimated variance with sim param
  estmlevect = nan(Nsim,ninterp);%% estimated variance with mle param
  estkervect = nan(Nsim,ninterp);%% estimated variance with ker param
  estmatvect = nan(Nsim,ninterp);%% estimated variance with mat param
  
  save(strcat('100sim_',modelo,'.mat'),'inivect1','inivect50',...
       'inivect100','estimate1','estimate50',...
       'estimate100','msemlevect','msekervect','msematvect',...
       'estsimvect','estmlevect','estkervect','estmatvect',...
       'lon','lat','intpoints')

else
  load(strcat('100sim_',modelo,'.mat'))
end

for nsim=nsimini:Nsim

%if length(find(nsim==[3 19 31]))==0
  testwt = wtvect(1);
  %  rand('state',sum(100*clock)); randn('state',sum(100*clock));
  rand('state',1000*nsim); randn('state',11000*nsim);
  newsim = true;
  pref100sim = strcat('100sim_wt001_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect1(nsim,:)=thetavect(1,:);
  estimate1(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  msemlevect(nsim,:) = msekrigmle;
  msekervect(nsim,:) = msekrigker;
  msematvect(nsim,:) = msekrigmat;
  estsimvect(1,:) = estvarsim;
  estmlevect(nsim,:) = estvarmle;
  estkervect(nsim,:) = estvarker;
  estmatvect(nsim,:) = estvarmat;
  dEmin = dE;
  disp(strcat('---------',num2str(nsim),'--','wtvect(1)---------'))
  
  testwt = wtvect(50);
  newsim = false;
  pref100sim = strcat('100sim_wt050_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect50(nsim,:)=thetavect(1,:);
  estimate50(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  if dE<dEmin
    dEmin = dE;    
    msemlevect(nsim,:) = msekrigmle;
    msekervect(nsim,:) = msekrigker;
    msematvect(nsim,:) = msekrigmat;
    estsimvect(1,:) = estvarsim;
    estmlevect(nsim,:) = estvarmle;
    estkervect(nsim,:) = estvarker;
    estmatvect(nsim,:) = estvarmat;
  end
  disp(strcat('---------',num2str(nsim),'--','wtvect(50)---------'))

  testwt = wtvect(100);
  newsim = false;
  pref100sim = strcat('100sim_wt100_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect100(nsim,:)=thetavect(1,:);
  estimate100(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  if dE<dEmin
    dEmin = dE;    
    msemlevect(nsim,:) = msekrigmle;
    msekervect(nsim,:) = msekrigker;
    msematvect(nsim,:) = msekrigmat;
    estsimvect(1,:) = estvarsim;
    estmlevect(nsim,:) = estvarmle;
    estkervect(nsim,:) = estvarker;
    estmatvect(nsim,:) = estvarmat;
  end
  disp(strcat('---------',num2str(nsim),'--','wtvect(100)---------'))

  %% save results for later analysis
  save(strcat('100sim_conv_',modelo,'.mat'),'inivect1','inivect50',...
       'inivect100','estimate1','estimate50',...
       'estimate100','msemlevect','msekervect','msematvect',...
       'estsimvect','estmlevect','estkervect','estmatvect',...
       'lon','lat','intpoints')

%end
end


