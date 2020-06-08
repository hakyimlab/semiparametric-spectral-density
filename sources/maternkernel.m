  %------------------------------------
  %% fit matern
  %------------------------------------
  matparam = log([.2,.0018,800]); %%log([smoo sig2 1000]);
  dife = 1e10; epsi = 1;
  paramini = matparam;
  fvalini = enermatern(paramini)
  while dife>epsi;
    [paramf fvalf exitflag] = fminsearch(@enermatern,paramini);
    dife = abs(fvalf - fvalini);
    disp(num2str(dife)); %DEBUG
    fvalini = fvalf;
    paramini = paramf;
  end
  matparam = paramf;
  matheta = exp(paramf);
  matsmoo = matheta(1);
  matsig2 = matheta(2);
  matrango = matheta(3);
  matirango = 2*sqrt(matsmoo)/matrango;
  matEi = enermatern(paramf);
  matwt = matirango;
  

  %% calc covariogram cloud and nadwat estimate of cov fn
  
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
  xxx = linspace(0,2/100,100);
  R = rmax/8; %% sieghankel will sample up to ~ 10*R
  [wexp,fcovario] = sieghankel(@kcovario,R);
  pp = spline(wexp,fcovario);
  fcovario = ppval(pp,xxx);
  
  %% calc likelihood using rohat positivedefined
  [kk,kpp] = ktilde(1,R); %kpp has spline info to interpolate pos def kcovario
  kcova = ppval(kpp,dista);
  if min(eig(kcova))>0
    cc = chol(kcova);
    invcc = inv(cc);
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    klike = -( nt*sum(log(diag(cc))) + 0.5 * xx'*xx );  %% NOT REML
  else
    disp('ktilde interpolated is not positive definite')
    klike=nan;
    save nopositive.mat
  end
  sig2k = ppval(kpp,0);
  %% save so far
  save rainmatker
