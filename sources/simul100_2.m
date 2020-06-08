% calls nogibbs to get Nsim simulations
%% of process and corresponding estimates
addpath('./mfunctions')

global Za dista M cnots indupper rr nt reml
global tailmat smoovect wtvect % tabulated values
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

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

nsimini=1;Nsim = 100;

testsplinetail = false; testmatern = false; testpolmatern =false; testexpon = false;
testsplinetail = true;
%testpolmatern = true;
%testmatern = true;
%testexpon = true;

%% simulated annealing parameters
L = 10000; %%max number of iterations
Ti = nt*100; %% scale temperature
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
  ll = 4-2;
elseif testpolmatern 
  modelo = 'polmatern';
  ll = 7-2;
elseif testmatern
  modelo = 'matern';
  ll = 3-1;
elseif testexpon
  modelo = 'expon';
  ll=3-1;
end

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat

load tablatailmat %% reads tailmat matrices
nparam = 2+ll+3+1;
fixwt = false;
gibbs = false; %% needed but not functional. it has to be false.
reml = false;

if nsimini ==1

  mlenuvect = nan(Nsim,1);
  mlesig2vect = nan(Nsim,1);
  mlebcoefvect = nan(Nsim,30);
  mlewtvect = nan(Nsim,1);

  likelivect = nan(Nsim,4);
  rmsecovavect = nan(Nsim,3);
  rmsespectvect = nan(Nsim,3);
  
  matnuvect = nan(Nsim,1);
  matsig2vect = nan(Nsim,1);
  matwtvect = nan(Nsim,1);
  
  kersig2vect = nan(Nsim,1);
  
  mlemspevect = nan(Nsim,ninterp); %% mse true kriged and mle kriged
  kermspevect = nan(Nsim,ninterp); %% mse true kriged and ker kriged
  matmspevect = nan(Nsim,ninterp); %% mse true kriged and mat kriged
  
  simkrigvarvect = nan(1,ninterp);   %% estimated variance with sim param
  mlekrigvarvect = nan(Nsim,ninterp);%% estimated variance with mle param
  kerkrigvarvect = nan(Nsim,ninterp);%% estimated variance with ker param
  matkrigvarvect = nan(Nsim,ninterp);%% estimated variance with mat param

  convcont = zeros(Nsim,1);
  save(strcat('100sim_conv_',modelo,'.mat'),'mlenuvect','mlesig2vect','mlebcoefvect',...
       'mlewtvect','kersig2vect','likelivect','rmsecovavect','rmsespectvect','matnuvect','matsig2vect','matwtvect',...
       'mlemspevect','kermspevect','matmspevect', 'simkrigvarvect',...
       'mlekrigvarvect','kerkrigvarvect','matkrigvarvect',  ...
       'lon','lat','intpoints','convcont','modelo')

else
  load(strcat('100sim_conv_',modelo,'.mat'))
end

%% make sure counter is set to zero for nsim>nsimini 
convcont(nsimini:end)=0;

ll0=ll; %% save initial ll in ll0

for nsim=nsimini:Nsim

  ll = ll0; %% reset ll
  bic = 1e100;
  newsim = true;
  
  for ii=1:7

    disp([nsim ll])
    testwt = wtvect(1);
    %  rand('state',sum(100*clock)); randn('state',sum(100*clock));
    rand('state',1000*nsim); randn('state',11000*nsim);
    pref100sim = strcat('100sim_wt001_',sprintf('%03d',nsim),'_'); 

    nogibbs2

    nza = numel(Za);
    penlike = - Ef - log(nza)*ll; %% BIC
    %penlike = - Ef - ll; %% AIC
    
    if - penlike < bic

      bic = - penlike;
      mlenuvect(nsim,1) = smoof;
      mlesig2vect(nsim,1) = sig2f;      
      mlebcoefvect(nsim,:) = nan;
      mlebcoefvect(nsim,1:(ll+3)) = bcoeff;
      mlewtvect(nsim,1) = wtf;

      likelivect(nsim,:) = [-truene -Ef klike -matEi];
      rmsecovavect(nsim,:)= [rmsecovam rmsecovak rmsecovamat];
      rmsespectvect(nsim,:) = [rmsespecm rmsespeck rmsespecmat];

      matnuvect(nsim,1) = matsmoo;
      matsig2vect(nsim,1) = matsig2;
      matwtvect(nsim,1) = matwt;
      
      kersig2vect(nsim,1) = sig2k;

      mlemspevect(nsim,:) = msekrigmle;
      kermspevect(nsim,:) = msekrigker;
      matmspevect(nsim,:) = msekrigmat;

      simkrigvarvect(nsim,:) = estvarsim;
      mlekrigvarvect(nsim,:) = estvarmle;
      kerkrigvarvect(nsim,:) = estvarker;
      matkrigvarvect(nsim,:) = estvarmat;
      
    end

    ll = ll+1;
    %% save results for later analysis
    save(strcat('100sim_conv_',modelo,'.mat'),'mlenuvect','mlesig2vect','mlebcoefvect',...
         'mlewtvect','likelivect','rmsecovavect','rmsespectvect','matnuvect','matsig2vect','matwtvect',...
         'mlemspevect','kersig2vect','kermspevect','matmspevect', 'simkrigvarvect',...
         'mlekrigvarvect','kerkrigvarvect','matkrigvarvect',  ...
         'lon','lat','intpoints','convcont','modelo')
  
  end
end




