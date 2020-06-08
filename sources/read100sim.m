clear all

%% of process and corresponding estimates
addpath('./mfunctions')

%load new100sim/100sim_splinetail.mat
%load mse_100sim/100sim_polmatern.mat
%load 100sim_conv_matern.mat
load mse_100sim/100sim_conv_expon.mat
indnan = find(~isnan(estimate1(:,1)));
nmax = length(indnan);
ebest = nan(nmax,size(estimate1,2));
inivect = nan(nmax,1);
%% ONLY 100!!
nmax = min(100,nmax);

%% choose best estimate from the three vectors (1,50,100) 
for ii=1:nmax
  estall=[estimate1(indnan(ii),:);estimate50(indnan(ii),:);estimate100(indnan(ii),:)];
  dE = estall(:,end-6) - estall(:,end-7);
  [temp,ix] = sort(dE);
  ebest(ii,:) = estall(ix(1),:); %% keep best estimate
  inivect(ii)=ix(1); %% which initial wt gave the min
end

%% sort by Emle - truene
dE = ebest(:,end-6)-ebest(:,end-7);
[dE,ix]=sort(dE);
ebest = ebest(ix,:);
%% good convergence if dE<2
conv = find(dE<20);

enu = ebest(conv,1);
esig2 = ebest(conv,2);
ebcoef = ebest(conv,3 : end-9);
ewt = ebest(conv,end-8);

etrue = ebest(conv,end-7);
emle = ebest(conv,end-6);
eker = -ebest(conv,end-5);

ermsecovm = ebest(conv,end-4);
ermsecovk = ebest(conv,end-3);
ermsespm = ebest(conv,end-2);
ermsespk = ebest(conv,end-1);

esig2k = ebest(conv,end);


%% corr nu a wt
%corrcoef(enu(conv),ewt(conv))

%% disp summary
disp(strcat('# of simu. conv = ',num2str(length(enu))))
disp('-----------')
disp(strcat('nu = ',num2str([mean(enu) std(enu)])))
disp(strcat('wt = ',num2str([mean(ewt) std(ewt)])))
disp(strcat('sig2 = ',num2str([mean(esig2) std(esig2)])))
disp(strcat('sig2k = ',num2str([mean(esig2k) std(esig2k)])))

disp(strcat('Etrue = ',num2str([mean(etrue) std(etrue)])))
disp(strcat('Emle = ',num2str([mean(emle) std(emle)])))
disp(strcat('Eker = ',num2str([mean(eker(~isnan(eker))) std(eker(~isnan(eker)))])))

disp(strcat('rmsecovm = ',num2str([mean(ermsecovm) std(ermsecovm)])))
disp(strcat('rmsecovk = ',num2str([mean(ermsecovk) std(ermsecovk)])))
disp(strcat('rmsespm = ',num2str([mean(ermsespm) std(ermsespm)])))
disp(strcat('rmsespk = ',num2str([mean(ermsespk) std(ermsespk)])))

%% plot summary
figure(1)
subplot(3,3,1);hist(enu);title('nu')
subplot(3,3,2);hist(ewt);title('wt')
subplot(3,3,3);hist(esig2);title('sig2')
subplot(3,3,4);hist(esig2k);title('sig2k')

subplot(3,3,5);plot(ermsecovk,ermsecovm,'.');title('rmse cov mle vs. ker')
subplot(3,3,6);plot(ermsespk,ermsespm,'.');title('rmse spec mle vs. ker')

subplot(3,3,7);hist(etrue);title('True Energy')
subplot(3,3,8);hist(emle);title('mle Energy')
subplot(3,3,9);hist(eker);title('kernel Energy')

print('-depsc','plot.eps')


msemle = msemlevect(conv,:); msemle = mean(msemle);
mseker = msekervect(conv,:); mseker = mean(mseker);
msemat = msematvect(conv,:); msemat = mean(msemat);
estvarsim = estsimvect(1,:);
estvarmle = estmlevect(conv,:); estvarmle = mean(estvarmle);
estvarker = estkervect(conv,:); estvarker = mean(estvarker);
estvarmat = estmatvect(conv,:); estvarmat = mean(estvarmat);

[lat,wnh4,lon,cwnh4]=textread('simul.dat','%f%f%f%f%*[^\n]');
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

ninterp = length(msemle);
figure(2)
plot(lon,lat,'.')
hold on;plot(intpoints(:,1),intpoints(:,2),'r+');hold off
for ii=1:ninterp
  text(intpoints(ii,1),intpoints(ii,2),num2str(round(mseker(ii)/msemle(ii))))
end

E0e1 = 1+msemle./estvarsim;
E0e2 = 1+mseker./estvarsim;
E0e3 = 1+msemat./estvarsim;
E1e1 = estvarmle./(estvarsim+msemle);
E2e2 = estvarker./(estvarsim+mseker);
E3e3 = estvarmat./(estvarsim+msemat);

%%  deviation from one
quantile([E0e1' E0e2' E0e3'],[.25 .5 .75])-1
quantile(abs([E1e1' E2e2' E3e3']-1),[.25 .5 .75])

