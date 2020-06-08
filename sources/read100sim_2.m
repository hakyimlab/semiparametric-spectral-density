clear all

%% of process and corresponding estimates
addpath('./mfunctions')

dire = 'remlnu3_ok_100sim/'; %% ok ordinary kriging = use est mean
%dire = '';

nmodelo = 1;
if nmodelo == 1 %'splinetail'
   load([dire '100sim_conv_splinetail.mat'])
elseif nmodelo == 2 %'polmatern'
   load([dire '100sim_conv_polmatern.mat'])
elseif nmodelo == 3 %'matern'
   load([dire '100sim_conv_matern.mat'])
elseif nmodelo == 4; %'expon'
   load([dire '100sim_conv_expon.mat'])
end
	    

%indnan = find(~isnan(estimate1(:,1)));
indnan = find(~isnan(estimate1(:,end-8-3))); % Kernel positive definite
nmax = length(indnan);
ebest = nan(nmax,size(estimate1,2));
inivect = nan(nmax,1);
%% ONLY 100!!%%nmax = min(100,nmax);

%% choose best estimate from the three vectors (1,50,100) 
for ii=1:nmax
  estall=[estimate1(indnan(ii),:);estimate50(indnan(ii),:);estimate100(indnan(ii),:)];
  dE = estall(:,end-6-3) - estall(:,end-7-3);
  [temp,ix] = sort(dE);
  ebest(ii,:) = estall(ix(1),:); %% keep best estimate
  inivect(ii)=ix(1); %% which initial wt gave the min
end

%% sort by Emle - truene
dE = ebest(:,end-9-3)-ebest(:,end-10-3);
[dE,ix]=sort(dE);
ebest = ebest(ix,:);
%% good convergence if dE<2
conv = find(dE<1e100);

enu = ebest(conv,1);
esig2 = ebest(conv,2);
ebcoef = ebest(conv,3 : end-12-3);
ewt = ebest(conv,end-11-3);
ematnu = ebest(conv,end-2);%mat
ematsig2 = ebest(conv,end-1);%mat
ematwt = ebest(conv,end);%mat

etrue = ebest(conv,end-10-3);
emle = ebest(conv,end-9-3);
eker = -ebest(conv,end-8-3);
emat = ebest(conv,end-7-3);

ermsecovm = ebest(conv,end-6-3);
ermsecovk = ebest(conv,end-5-3);
ermsecovmat = ebest(conv,end-4-3);
ermsespm = ebest(conv,end-3-3);
ermsespk = ebest(conv,end-2-3);
ermsespmat = ebest(conv,end-1-3);

esig2k = ebest(conv,end-3);


%% corr nu a wt
%corrcoef(enu(conv),ewt(conv))

%% disp summary
disp(modelo)
disp(strcat('# of simu. conv = ',num2str(length(enu))))
disp('-----------')

disp(strcat('nu = ',num2str([mean(enu) std(enu)])))
disp(strcat('matnu = ',num2str([mean(ematnu) std(ematnu)])))

disp(strcat('sig2 = ',num2str([mean(esig2) std(esig2)])))
disp(strcat('matsig2 = ',num2str([mean(ematsig2) std(ematsig2)])))
disp(strcat('sig2k = ',num2str([mean(esig2k) std(esig2k)])))

disp(strcat('wt = ',num2str([mean(ewt) std(ewt)])))
disp(strcat('matwt = ',num2str([mean(ematwt) std(ematwt)])))

disp(strcat('rmsecovm = ',num2str([mean(ermsecovm) std(ermsecovm)])))
disp(strcat('rmsecovmat = ',num2str([mean(ermsecovmat) std(ermsecovmat)])))
disp(strcat('rmsecovk = ',num2str([mean(ermsecovk) std(ermsecovk)])))

disp(strcat('rmsespm = ',num2str([mean(ermsespm) std(ermsespm)])))
disp(strcat('rmsespmat = ',num2str([mean(ermsespmat) std(ermsespmat)])))
disp(strcat('rmsespk = ',num2str([mean(ermsespk) std(ermsespk)])))

disp(strcat('liketrue = ',num2str([-mean(etrue) std(etrue)])))
disp(strcat('likemle-true = ',num2str([-mean(emle)+mean(etrue) std(emle)])))
disp(strcat('likemat-true = ',num2str([-mean(emat)+mean(etrue) std(emat)])))
disp(strcat('likeker-true = ',num2str([-mean(eker(~isnan(eker)))+mean(etrue) std(eker(~isnan(eker)))])))

%% plot summary
%figure(1)
%subplot(3,3,1);hist(enu);title('nu')
%subplot(3,3,2);hist(ewt);title('wt')
%subplot(3,3,3);hist(esig2);title('sig2')
%subplot(3,3,4);hist(esig2k);title('sig2k')

%subplot(3,3,5);plot(ermsecovk,ermsecovm,'.');title('rmse cov mle vs. ker')
%subplot(3,3,6);plot(ermsespk,ermsespm,'.');title('rmse spec mle vs. ker')

%subplot(3,3,7);hist(etrue);title('True Energy')
%subplot(3,3,8);hist(emle);title('mle Energy')
%subplot(3,3,9);hist(eker);title('kernel Energy')

%print('-depsc','plot.eps')


msemle = msemlevect(indnan,:); msemle = mean(msemle(conv,:));
mseker = msekervect(indnan,:); mseker = mean(mseker(conv,:));
msemat = msematvect(indnan,:); msemat = mean(msemat(conv,:));

estvarsim = estsimvect(1,:);
estvarmle = estmlevect(indnan,:); estvarmle = mean(estvarmle(conv,:));
estvarker = estkervect(indnan,:); estvarker = mean(estvarker(conv,:));
estvarmat = estmatvect(indnan,:); estvarmat = mean(estvarmat(conv,:));

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
%figure(2)
plot(lon,lat,'.')
%hold on;plot(intpoints(:,1),intpoints(:,2),'r+');hold off
for ii=1:ninterp
%  text(intpoints(ii,1),intpoints(ii,2),num2str(round(mseker(ii)/msemle(ii))))
  text(intpoints(ii,1),intpoints(ii,2),num2str(round(1000*estvarmle(ii)) ...
                                               ))
end

E0e0 = estvarsim;
E0e1 = 1+msemle./E0e0;
E0e2 = 1+mseker./E0e0;
E0e3 = 1+msemat./E0e0;
E1e1 = estvarmle./(E0e0+msemle);
E2e2 = estvarker./(E0e0+mseker);
E3e3 = estvarmat./(E0e0+msemat);

%% E0e12/E0e02 - deviation from one
disp([quantile(E0e1',.5)-1 quantile(E0e1',.75)-quantile(E0e1',.25)])
disp([quantile(E0e3',.5)-1 quantile(E0e3',.75)-quantile(E0e3',.25)])
disp([quantile(E0e2',.5)-1 quantile(E0e2',.75)-quantile(E0e2',.25)])


disp([ quantile(abs(E1e1'-1), .5) quantile(abs(E1e1' -1), .75)- quantile(abs(E1e1'-1),.25)] )
disp([ quantile(abs(E3e3'-1), .5) quantile(abs(E3e3' -1), .75)- quantile(abs(E3e3'-1),.25)] )
disp([ quantile(abs(E2e2'-1), .5) quantile(abs(E2e2' -1), .75)- quantile(abs(E2e2'-1),.25)] )

