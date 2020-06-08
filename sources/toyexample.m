%% calls nogibbs to get Nsim simulations
%% of process and corresponding estimates
clear all
addpath('./mfunctions')

%% observation points
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

rango = 200;
irango = 2*sqrt(0.5)/rango;
simcova = exp(-dista*irango);

ksim = exp(-dist2obs*irango);
kgauss = exp(-(dist2obs*irango).^2/2);
k2e = exp(-2*dist2obs*irango)/2;

invKsim = inv(simcova);
invKgauss = inv(exp(-(dista*irango).^2/2));
invK2e = inv(exp(-2*dista*irango)/2);

lambsim = (ksim' * invKsim)';
lambgauss = (kgauss' * invKgauss)';
lamb2e = (k2e' * invK2e)';

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

nsim=1;
rand('state',1000*nsim); randn('state',11000*nsim);
for ii=1:nt
  Za(:,ii) = simulachol(simcova);
end

%% interpolate to intpoints

krigZgauss = (lambgauss - lambsim)'* Za;
krigZ2e = (lamb2e - lambsim)'* Za;

estvarsim = 1 - diag(lambsim'*ksim);
estvargauss = 1 - diag(lambgauss'*kgauss);
estvar2e = .5 - diag(lamb2e'*k2e);

msekriggauss = mean(krigZgauss.^2,2);
msekrig2e = mean(krigZ2e.^2,2);

disp(['rango:' num2str(rango)])
disp(['e(-2t)/2: ' num2str(1+median(msekrig2e./estvarsim))])
disp(['gauss: ' num2str(1+median(msekriggauss./estvarsim))])


disp(['E1e1/E0e1:' num2str(rmse(estvar2e./(estvarsim+msekrig2e),1))] )
disp(['E2e2/E0e2 gauss:' num2str(rmse(estvargauss./(estvarsim+msekriggauss),1)) ] ) 
      

xx=(0:100)/100*5;
c0=exp(-xx);
c1=exp(-xx*2)/2;
c2=exp(-xx.^2/2);
rmse(c0,c1);
rmse(c0,c2);