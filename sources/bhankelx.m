function [covacnots,cova0] = bhankelx(cnots,smoo,bcoef,wt)
% res=bhankel(r,smoo,bcoef,wt)
% calc hankel tr of sum of bsplines fi*Bi
% and algebraic tail
% 2 pi is not included since I'm going to normalize by cova0
% knots are the breaks in deBoor's notation
  
global Hpol tailmat smoovect wtvect
  
%% some preliminary definitions and calculations
gm = 2*smoo+2;
gmvect = smoovect*2 + 2;
ll = length(bcoef)-3;
nnots = length(cnots);
t = nodos(wt,ll);
indwt = find(wtvect==wt);

sp = spmak(t,bcoef);
ft = spval(sp,wt);  % continuity ft = bcoef(ll+1)/6 + bcoef(ll+2)*2/3 + bcoef(ll+3)/6;
consta = ft * wt^gm;

%% coeff of polynomials from bspline coeffs
A = coefmat(t,bcoef); % coef of piecewise polynomials
Aext = repmat(A,[1,1,nnots+1]);

covacnots = sum(sum(Hpol(:,:,:,indwt).*Aext,2),1); %% calc hankel tr. of piecewise pol
cova0 = covacnots(end); 
covacnots = covacnots(1:(end-1)); 
covacnots = reshape(covacnots,1,nnots);

giir = find(gmvect>=gm);
giir = giir(1);
giil = find(gmvect<=gm);
giil = giil(end);

%% add tail integrals
cova0 = cova0 - consta*wt^(2-gm)/(2-gm); % add the tail

tailvect = ( tailmat(:,giil,indwt) + tailmat(:,giir,indwt) )' / 2 ;
covacnots = covacnots + consta * tailvect .* ( cnots.^(gm-2) ) ;

%%======================================
covacnots = covacnots/cova0; %% normalize so cova(0)=1

%%======================================
function res=coefmat(t,bcoef)
% calc ppform coef of from bspline form
% 
ll = length(t) - 7;
res = mexbsplpp(t,bcoef);
res = res./repmat(cumprod([1 1:3]'),1,ll);
res = res(4:-1:1,:)';
