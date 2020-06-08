%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------

global Hpol 

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=1); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=.1;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= 1); wt = wt(end);
t = nodos(wt,ll);

newtheta = matheta;
newsmoo = smoovect(smoovect<=newtheta(1));newsmoo=newsmoo(end);
newsig2 = newtheta(2);
newirango = 2*sqrt(newsmoo)/newtheta(3);
newwt = wtvect(wtvect <= newirango*2); newwt = newwt(end);

theta = [newsmoo newsig2 bcoef newwt];thetaini=theta;

load(strcat('raintablahpol_',num2str(ll)))
Ei = ener(theta); Eini = Ei;
    
%% Temperature
L = 1000; 

Ti = nt*3000; %% scale temperature
             %lam = 0.99;
             %tsca = Ti * lam.^(1:L);
             %% cooling schedule
tsca = ones(L,1);
a=30;
for kk=1:L;
  tsca(kk+1) = tsca(kk)/(1+a*tsca(kk));
end
tsca = tsca/tsca(1) * Ti;

%------------------------------------
%% loop for simulated annealing
%------------------------------------
tic
conta=1;
llvect = nan(L+1,1);
llvect(1) = Ei;
thetavect = nan(L+1,length(theta));
thetavect(1,:)=theta;

for lind = 1:L
  T=tsca(lind);
  % generate new state %logf m-vector
  thetaf = bproposal(theta,fixwt,gibbs);
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
  thetavect(lind+1,:)=theta;
  
end
toc

%% revert back to last updated values

thetaf = theta;
Ef = Ei; 
smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll);
covaf = calcCova(thetaf);

% display some values

disp('ini  matern  final')
disp([smoo matsmoo smoof ])
disp([sig2 matsig2 sig2f])
disp([wt matwt wtf])
disp(round([ Eini matEi Ei]))
disp(strcat('number of iterations =  ',num2str(L)));

%% calc cov and spect densities

[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);
mlespect = spectralfn(xxx,smoof,tf,bcoeff)/cova0;
matspect = fmatern(xxx,matirango,matsmoo);
matcovacnots = matsig2*hmatern(cnots,matrango,matsmoo);

%% plot results

fighandle = figure(1);
set(fighandle,'PaperPosition',[0.25 2.5 8 3.5]);
subplot(1,2,1)
plot(cnots,covacnots*sig2f,'r.'); hold on; 
plot(cnots,matcovacnots,'g+','markersize',3);
grid
plot(cnots,nadwat(cnots,covario,rsort,50),'kx','markersize',3);
plot(cnots,covacnots*sig2f,'r.');
hold off;

legend('ML S+T','ML Matern','Kernel')
xlabel('distance');title('cov function')

subplot(1,2,2)
plot(xxx,mlespect,'r.'); hold on
plot(xxx,matspect,'g.');
plot(xxx,fcovario/2/pi,'kx','markersize',3);
plot(xxx,mlespect,'r.'); hold off
legend('ML S+T','ML Matern','Kernel')
grid
xlabel('freq')
ylabel('sqrt(f(w))')
title('SQRT Spectral density')

%print( '-depsc',epsfile)

%% interpolate to intpoints

cc = spline(cnots,covacnots);
kmle = sig2f * ppval(cc,dist2obs);
kker = ppval(kpp,dist2obs);
kmat = matsig2 * hmatern(dist2obs,matrango,matsmoo);
invKmle = inv(covaf);
invKker = inv(kcova);
matcova = matsig2*hmatern(dista,matrango,matsmoo);
invKmat = inv(matcova);

M = ones(nobs,1);
if reml
  bmle = 1 - M'*invKmle*kmle; invWmle = inv(M'*invKmle*M);
  bker = 1 - M'*invKker*kker; invWker = inv(M'*invKker*M);
  bmat = 1 - M'*invKmat*kmat; invWmat = inv(M'*invKmat*M);
else
  bmle =zeros(1,ninterp);bker=zeros(1,ninterp);bmat=zeros(1,ninterp);
  invWmle=0;invWker=0;invWmat=0;
end

lambmle = kmle' * invKmle + bmle'*invWmle*M'*invKmle; lambmle=lambmle';
lambker = kker' * invKker + bker'*invWker*M'*invKker; lambker=lambker';
lambmat = kmat' * invKmat + bmat'*invWmat*M'*invKmat; lambmat=lambmat';

%% plugin kriging variance
estvarmle = sig2f - diag(lambmle'*kmle) + diag(bmle'*invWmle*bmle);
estvarker = sig2k - diag(lambker'*kker)+ diag(bker'*invWker*bker);
estvarmat = matsig2 - diag(lambmat'*kmat) + diag(bmat'*invWmat*bmat);






