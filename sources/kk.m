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
