%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------

%% diff with nogibbs.m = kcovario is replaced by ktilde which 
%% is positive definite

global Hpol 

%------------------------------------
% generate observations
%------------------------------------
if newsim == true
  simsmoo = 3;
  indi = find(smoovect>=simsmoo); simsmoo = smoovect(indi(1)); %% the first smoo>=3
  simgm=2*(simsmoo+1);
  simwt = 0.0100; %% simwt=0.0038
  simwt = wtvect(find(wtvect <= simwt)); simwt = simwt(end);
  simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
  simsig2=1;
  simll = 4;
  
  if testsplinetail
    load(strcat('tablahpol_',num2str(simll)))
    simt = nodos(simwt,simll);
    simknots = simt(4:end-3);
    %  simbcoef = [0.8 1.0 .8 .2  1.2 1 .8 .6 .4]; %% ll=6 should have ll+3 elements
    simbcoef = [0.2 1 .2  2  .6 .4 .2]; %% ll=4
    simbcoef(end)=contderiv(simbcoef,simgm);
    simtheta = [simsmoo simsig2 simbcoef simwt];
    simcova = calcCova(simtheta);
  elseif testpolmatern 
    % test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
    u = 0.5*simirango; v=.1*simirango;
    simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);
  elseif testmatern
    % test matern cova
    simcova = simsig2*hmatern(dista,simrango,simsmoo);
  elseif testexpon
    % test exponential spectral density
    simcova = simsig2*((dista*simirango).^2+1).^(-1.5);
  end
  
  for ii=1:nt
    Za(:,ii) = simulachol(simcova);
  end
  
  if testsplinetail
    save(strcat('splinetail_',sprintf('%03d',nt),'_',num2str(ll),'.mat'),...
         'Za','simbcoef','simll','simsig2','simsmoo','simgm','simwt',...
         'simt','simknots', 'simtheta','simcova')
  elseif testpolmatern 
    save(strcat('polmatern_',sprintf('%03d',nt),'_','.mat'),...
         'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango','u','v',...
         'simcova')
  elseif testmatern
    save(strcat('matern_',sprintf('%03d',nt),'.mat'),...
         'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango',...
         'simcova')
  elseif testexpon
    save(strcat('expon_',sprintf('%03d',nt),'.mat'),...
         'Za','simsig2','simwt','simrango','simirango',...
         'simcova')
  end
  
else
  
  %% read observations
  if testsplinetail
    simll = 4;
    load(strcat('splinetail_',sprintf('%03d',nt),'_',num2str(simll),'.mat'));
  elseif testpolmatern 
    load(strcat('polmatern_',sprintf('%03d',nt),'_','.mat'))
  elseif testmatern
    load(strcat('matern_',sprintf('%03d',nt),'.mat'))
  elseif testexpon
    load(strcat('expon_',sprintf('%03d',nt),'.mat'))
  end

end

cc = chol(simcova);
invcc = inv(cc);

if reml
  invcova = invcc * invcc';
  W = M' * invcova * M;
  ww = chol(W);
  invW = inv(W);
  y = nt * (sum(log(diag(cc))) + sum(log(diag(ww))));
  SS = invcova * (eye(nobs) - M * invW * M' * invcova);
  for ii = 1:nt
    z = Za(:,ii);
    y = y + 0.5 * z' * SS * z;
  end  
  truene = y;
else
  xx = invcc' * Za;
  nobs = size(Za,1);
  xx=reshape(xx,nt*nobs,1);
  truene = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
end

clear('xx','cc');

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=1); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=.1;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= 1); wt = wt(end);

t = nodos(wt,ll);
theta = [smoo,sig2,bcoef,wt];thetaini=theta;

load(strcat('tablahpol_',num2str(ll)))
Ei = ener(theta); Eini = Ei;


%------------------------------------
%% improve initial guess fitting matern
%------------------------------------
if newsim
  matparam = log([smoo sig2 1000]);
  %[paramf fval exitflag] = fminunc(@enermatern,matheta)
  maternini = true;
  dife = 1e10; epsi = 1;
  disp(num2str(dife));
  if maternini
    paramini = matparam;
    fvalini = Ei;
    while dife>epsi;
      [paramf fvalf exitflag] = fminsearch(@enermatern,paramini);
      dife = abs(fvalf - fvalini);
      disp(num2str(dife)); %DEBUG
      fvalini = fvalf;
      paramini = paramf;
    end
    newtheta = exp(paramf);
    newsmoo = smoovect(smoovect<=newtheta(1));newsmoo=newsmoo(end);
    newsig2 = newtheta(2);
    newirango = 2*sqrt(newsmoo)/newtheta(3);
    newwt = wtvect(wtvect <= newirango*2); newwt = newwt(end);
    theta = [newsmoo newsig2 bcoef newwt];
    matheta = exp(paramf);
    matsmoo = matheta(1);
    matsig2 = matheta(2);
    matrango = matheta(3);
    matirango = 2*sqrt(matsmoo)/matrango;
    matEi = enermatern(paramf);
    matwt = matirango;
    Ei = ener(theta); 
    
    %% DEBUG, wt is set not = matwt
    %testwt = wtvect(wtvect <= simwt); testwt = testwt(end);
    if ~testmatern
      theta(end) = testwt;
    else
      %theta(end) = testwt;
      %theta(end) = wtvect(50);
    end
    
  else
    matEi = nan;
    matwt = nan;
  end
end
%------------------------------------
%% loop for simulated annealing
%------------------------------------

conta=1;
llvect = nan(L+1,1);
llvect(1) = Ei;
thetavect = nan(L+1,length(theta));
thetavect(1,:)=theta;

for lind = 1:L
  T=tsca(lind);
  % generate new state %logf m-vector
  if ~testmatern
    thetaf = bproposal(theta,fixwt,gibbs);
  else
    thetaf = bproposalmatern(theta,fixwt,gibbs);
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
  thetavect(lind+1,:)=theta;
  
end
thetaf = theta; %% revert back to last updated values
Ef = Ei; %% revert back to last updated values

smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll);
covaf = calcCova(thetaf);

disp('number of indep replicates')
disp(nt)
disp('sim  ini  matern  final')
disp([simsmoo smoo matsmoo smoof ])
disp([simsig2 sig2 matsig2 sig2f])
disp([simwt wt matwt wtf])
disp(round([truene Eini matEi Ei]))
dE = Ei - truene;
disp('estimated ener - true ener'); disp(dE);
disp(strcat('number of iterations =  ',num2str(L)));
%disp('lam -- cooling schedule T=lam^(1:iter))');disp(lam)


[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);

%% calc covariogram cloud and nadwat estimate of cov fn
if newsim

  dindupper = find(triu(ones(size(dista))));
  
  ncov = length(dindupper);
  covario = nan(ncov,nt);
  [rsort,indsort] = sort(dista(dindupper));
  
  if reml
    zm = mean(Za(:));
  else
    zm = 0;
  end
  
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
  meanrmin = 50; 
  xxx = linspace(0,2*wtf,100);
  R = rmax/8; %% sieghankel will sample up to ~ 10*R
  [wexp,fcovario] = sieghankel(@kcovario,R);
  pp = spline(wexp,fcovario);
  fcovario = ppval(pp,xxx);
  
end

if testsplinetail
  load(strcat('tablahpol_',num2str(simll)))  
  [simcovacnots,simcova0] = bhankelx(cnots,simsmoo,simbcoef,simwt);
  simspectral = spectralfn(xxx,simsmoo,simt,simbcoef)/simcova0;
  prefijo = strcat(num2str(nt),'_',num2str(simll),'_',...
                    num2str(ll));
  cc = spline(cnots,simcovacnots);
  ksim = simsig2 * ppval(cc,dist2obs);
  savefile = strcat('opt_st_',prefijo,'.mat');
elseif testpolmatern 
  simcovacnots = polmatern(cnots,u,v,simirango,simsmoo);
  simspectral = fpolmatern(xxx,u,v,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_', num2str(ll));
  ksim = simsig2 * polmatern(dist2obs,u,v,simirango,simsmoo);
  savefile = strcat('opt_polmat_',prefijo,'.mat');
elseif testmatern
  simcovacnots = simsig2*hmatern(cnots,simrango,simsmoo);
  simspectral = fmatern(xxx,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_',num2str(ll));
  ksim = simsig2 * hmatern(dist2obs,simrango,simsmoo);
  savefile = strcat('opt_matern_',prefijo,'.mat');
elseif testexpon
  simcovacnots = simsig2*((cnots*simirango).^2 + 1).^(-1.5);
  constalfa = 1/simirango^2/2/pi;
  simspectral = constalfa * exp(-xxx /simirango);
  prefijo = strcat(num2str(nt),'_',num2str(ll));
  ksim = simsig2*( (dist2obs*simirango).^2 + 1 ).^(-1.5);
  savefile = strcat('opt_expon_',prefijo,'.mat');
end
epsfile = strcat(pref100sim,modelo,'_',prefijo,'.eps');

%% calc likelihood using rohat positivedefined
[kk,kpp] = ktilde(1,R); %kpp has spline info to interpolate pos def kcovario
if newsim
  kcova = ppval(kpp,dista);
  if min(eig(kcova))>0
    cc = chol(kcova);
    invcc = inv(cc);
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    klike = -( nt*sum(log(diag(cc))) + 0.5 * xx'*xx );
  else
    disp('ktilde interpolated is not positive definite')
    klike=nan;
    save nopositive.mat
  end
end

mlespect = spectralfn(xxx,smoof,tf,bcoeff)/cova0;
matspect = fmatern(xxx,matirango,matsmoo);
matcovacnots = matsig2*hmatern(cnots,matrango,matsmoo);
rmsecovak = rmse(ppval(kpp,cnots),simcovacnots);%%check this
rmsecovam = rmse(covacnots,simcovacnots);
rmsecovamat = rmse(matcovacnots,simcovacnots);
rmsespeck = rmse(fcovario/2/pi,simspectral);
rmsespecm = rmse(mlespect,simspectral);
rmsespecmat = rmse(matspect,simspectral);

sig2k = ppval(kpp,0);

%% plot results

fighandle = figure(1);
set(fighandle,'PaperPosition',[0.25 2.5 8 3.5]);
subplot(1,2,1)
plot(cnots,simcovacnots)
hold on; plot(cnots,covacnots*sig2f,'r.');
plot(cnots,matcovacnots,'g+','markersize',3);
grid
%plot(rsort,covario,'.','markersize',2);
plot(cnots,nadwat(cnots,covario,rsort,50),'kx','markersize',3);
hold on; plot(cnots,covacnots*sig2f,'r.');
plot(cnots,simcovacnots,'b');

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
text(rx+dx,ry+dy,sprintf('%4.0f',-truene+Ei))
text(rx+2*dx,ry+dy,sprintf('%4.0f',-truene-klike))

text(rx,ry+2*dy,sprintf('%4.3f',0))
text(rx+dx,ry+2*dy,sprintf('%4.3f',rmsecovam))
text(rx+2*dx,ry+2*dy,sprintf('%4.3f',rmsecovak))

text(rx,ry+3*dy,sprintf('%3.2f',simsmoo))
text(rx+dx,ry+3*dy,sprintf('%3.2f',smoof))
text(rx+2*dx,ry+3*dy,'NA')

text(rx,ry+4*dy,sprintf('%3.2f',simsig2))
text(rx+dx,ry+4*dy,sprintf('%3.2f',sig2f))
text(rx+2*dx,ry+4*dy,sprintf('%3.2f',sig2k) )

text(rx,ry+5*dy,sprintf('%5.4f',simwt) )
text(rx+dx,ry+5*dy,sprintf('%5.4f',wtf) )
text(rx+2*dx,ry+5*dy,'NA')

legend('True','ML S+T','ML Matern','Kernel')
xlabel('distance');title('cov function')

subplot(1,2,2)
plot(xxx,sqrt(simspectral));hold on
plot(xxx,sqrt(mlespect),'r.');
plot(xxx,sqrt(matspect),'g.');
plot(xxx,sign(fcovario).*sqrt(abs(fcovario/2/pi)),'kx','markersize',3);
plot(xxx,sqrt(mlespect),'r.');
plot(xxx,sqrt(simspectral));hold off
legend('True','ML S+T','ML Matern','Kernel')
grid
xlabel('freq')
ylabel('sqrt(f(w))')
title('SQRT Spectral density')

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
delete(fighandle)

%% interpolate to intpoints
% ksim defined before
cc = spline(cnots,covacnots);
kmle = sig2f * ppval(cc,dist2obs);
kker = ppval(kpp,dist2obs);
kmat = matsig2 * hmatern(dist2obs,matrango,matsmoo);
invKsim = inv(simcova);
invKmle = inv(covaf);
invKker = inv(kcova);
matcova = matsig2*hmatern(dista,matrango,matsmoo);
invKmat = inv(matcova);

M = ones(nobs,1);
if reml
  bsim = 1 - M'*invKsim*ksim; invWsim = inv(M'*invKsim*M);
  bmle = 1 - M'*invKmle*kmle; invWmle = inv(M'*invKmle*M);
  bker = 1 - M'*invKker*kker; invWker = inv(M'*invKker*M);
  bmat = 1 - M'*invKmat*kmat; invWmat = inv(M'*invKmat*M);
else
  bsim=zeros(1,ninterp);bmle =zeros(1,ninterp);bker=zeros(1,ninterp);bmat=zeros(1,ninterp);
  invWsim=0;invWmle=0;invWker=0;invWmat=0;
end

lambsim = ksim' * invKsim + bsim'*invWsim*M'*invKsim; lambsim=lambsim';
lambmle = kmle' * invKmle + bmle'*invWmle*M'*invKmle; lambmle=lambmle';
lambker = kker' * invKker + bker'*invWker*M'*invKker; lambker=lambker';
lambmat = kmat' * invKmat + bmat'*invWmat*M'*invKmat; lambmat=lambmat';

%% error for kriging with wrong cov function
krigZmle = (lambmle' - lambsim')* Za;
krigZker = (lambker' - lambsim')* Za;
krigZmat = (lambmat' - lambsim')* Za;

%% plugin kriging variance
estvarsim = simsig2 - diag(lambsim'*ksim) + diag(bsim'*invWsim*bsim);
estvarmle = sig2f - diag(lambmle'*kmle) + diag(bmle'*invWmle*bmle);
estvarker = sig2k - diag(lambker'*kker)+ diag(bker'*invWker*bker);
estvarmat = matsig2 - diag(lambmat'*kmat) + diag(bmat'*invWmat*bmat);

%% extra mse, mean over replications
msekrigmle = mean(krigZmle.^2,2);
msekrigker = mean(krigZker.^2,2);
msekrigmat = mean(krigZmat.^2,2);

%% clear some variables and save
if newsim
  clear('xx','cc');
  clear('Hpol','tailmat')
  save(savefile);
end





