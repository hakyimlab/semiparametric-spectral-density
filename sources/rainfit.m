clear all

addpath('./mfunctions')

global Za dista M cnots indupper rr nt reml
global tailmat smoovect wtvect % tabulated values
global covario rsort rmax;

  %% read data
  
  [ id1, lat1, lon1,   ele1, x1, y1    , mean1 , resi1] = ...
      textread('resi_train.txt', '%d%f%f%f%f%f%f%f%*[^\n]','headerlines',1);
  %%'%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*[^\n]');
  nobs = length(id1);
  dista = eudist([x1 y1]);
  indupper = find(triu(dista,1)); 
  rr = dista(indupper); rmin = min(rr); rmax = max(rr);
  nnots=100; cnots = linspace(rmin,rmax,nnots);
  
  %% interpolation points
  
  [ id2,lat2,   lon2,   ele2, x2, y2    , mean2 , resi2] = ...
      textread('resi_test.txt', '%d%f%f%f%f%f%f%f%*[^\n]','headerlines',1);
  ninterp = length(id2);
  dist2obs = eudist([x1 y1],[x2 y2]);
  
  Za = resi1;
  nt = 1; %number of replicates in time




%% fit matern y kernel

newdata = true;

if newdata

  maternkernel

else

  load rainmatker  
  
end


%% smoothness and wt values for tabulation

ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);



%% initial number of nodes

ll = 2;
ll0=ll; %% save initial ll in ll0

%% gen hpol tables
%for ii=2:13
%  gen_tables(ii,rmin,rmax,nnots,nwt)
%end

%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save raintablatailmat tailmat

load raintablatailmat %% reads tailmat matrices

nparam = 2+ll+3+1;
fixwt = false;
gibbs = false; %% needed but not functional. it has to be false.
reml = false;
nugget = true;

bic = 1e100;

for ii=1:5 %%11
  
  disp([ ll])
  nogibbsrain
  
  nza = numel(Za);
  %%    penlike = - Ef - log(nza)*ll; %% BIC
  penlike = - Ef - ll; %% AIC
  
  if - penlike < bic
    
    bic = - penlike;
    mlenu = smoof;
    mlesig2 = sig2f;
    mlebcoef = bcoeff;
    mlewt = wtf;
    
    mlelike = -Ef;
    kerlike = klike;
    matlike = -matEi;
    
    matnu = matsmoo;
    matsig2_aic = matsig2;
    matwt_aic = matwt;  
    
    kersig2 = sig2k;
    
    save rainfittemp
  end
  
  ll = ll+1;
  
  
end


clear tailmat Hpol
save rainfitrun



