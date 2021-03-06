bhankelx.m                                                                                          0000644 0002715 0000074 00000003014 10343675065 011513  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function [covacnots,cova0] = bhankelx(cnots,smoo,bcoef,wt)
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    bproposal.m                                                                                         0000644 0002715 0000074 00000005603 10304651703 011714  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function thetaf=bproposal(theta,fixwt,gibbs)
%% thetaf=bproposal(theta) 
%% fknots is in actual scale, not log scale

global smoovect wtvect

%randn('state',sum(100*clock))
%rand('state',sum(100*clock))


gmvect = 2*smoovect+2;  
ngm = length(smoovect);

%rangof = rrango(theta(1));
smoof = rsmoo(theta(1),ngm,smoovect,gibbs);
gmf = (smoof+1)*2;
sig2f = rsig2(theta(2),gibbs);

bcoef = theta(3:end-1);bcoeff = bcoef; 
bcoeff(2:end-1) = rbcoef(bcoef(2:end-1),gibbs);
bcoeff(1)=bcoeff(3);
bcoeff(end) = contderiv(bcoeff,gmf);

if fixwt
  wtf = theta(end);
else
  nwt = length(wtvect);
  wtf = rwt(theta(end),nwt,wtvect,gibbs);
end

thetaf = [smoof,sig2f,bcoeff, wtf];
    
%function y = rsmoo(smoo)
%    p1=1/3;p2=1/3;p3=1-p1-p2;s1=.1;s2=1;
%    y = exp(mixture(log(smoo),p1,p2,p3,s1,s2));
%    while (y > 5)
%      y = exp(mixture(log(smoo),p1,p2,p3,s1,s2));
%    end
 
function y = rsmoo(smoo,ngm,smoovect,gibbs)
  ns = find(smoovect==smoo);
  if gibbs
    p1=.9;p2=.1;p3=1-p1-p2;s1 = .1;  %% gibbs
  else
    %    p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;
    p1=.2;p2=.1;p3=1-p1-p2;s1 = .1; 
  end
  u = rand(1);
  if u<p1
    ini = max(1,ns-s1*ngm);
    fin = min(ngm, ns+s1*ngm);
    y =  sample(smoovect(ini:fin));
  elseif u < p2 + p1
    ini = 1;
    fin = ngm;
    y = sample(smoovect(ini:fin));
  else
    y = smoo;
  end

function y = rwt(wt,nwt,wtvect,gibbs)
  ns = find(wtvect==wt);
%  p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;
   if gibbs
     p1=.8;p2=.2;p3=1-p1-p2;s1 = .1;  
   else
     %p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;  
     p1=.2;p2=.1;p3=1-p1-p2;s1 = .1; 
   end
%  p1=.7;p2=.1;p3=1-p1-p2;s1 = .1; 
  u = rand(1);
  if u<p1
    ini = max(1,ns-s1*nwt);
    fin = min(nwt, ns+s1*nwt);
    y =  sample(wtvect(ini:fin));
  elseif u < p2 + p1
    ini = 1;
    fin = nwt;
    y = sample(wtvect(ini:fin));
  else
    y = wt;
  end

function y = rsig2(x,gibbs)
    %p1=1/3;p2=1/3;p3=1-p1-p2;s1=.1;s2=1;
    if gibbs
      p1=.7;p2=.3;p3=1-p1-p2;s1=.1;s2=1;  %gibbs
    else
      %p1=.3;p2=.2;p3=1-p1-p2;s1=.1;s2=1; 
      p1=.2;p2=.1;p3=1-p1-p2;s1 = .1;s2=1; 
    end
    y = exp(mixture(log(x),p1,p2,p3,s1,s2));

function y = rbcoef(x,gibbs)    
    %p1=1/3;p2=1/3+;p3=1-p1-p2;s1=.1;s2=1;
    %p1=.5;p2=0;p3=1-p1-p2;s1=1;s2=1;
    %p1=.3;p2=.5;p3=1-p1-p2;s1=.1;s2=1;
    %p1=.35;p2=.4;p3=1-p1-p2;s1=.1;s2=1; %well with polmatern
    %p1=.4;p2=.4;p3=1-p1-p2;s1=.1;s2=1;
    if gibbs
      p1=.6;p2=.4;p3=1-p1-p2;s1=.1;s2=1; %% gibbs
    else
      %p1=.3;p2=.2;p3=1-p1-p2;s1=.1;s2=1; 
      p1=.2;p2=.1;p3=1-p1-p2;s1 = .1;s2=1; 
    end
    y = x;
    for ii=1:length(y)
        y(ii) = exp(mixture(log(x(ii)),p1,p2,p3,s1,s2));    
    end

function y=mixture(x,p1,p2,p3,s1,s2,m2)
    if nargin==6
        m2=0;
    end
    u=rand(1);
    if u<p1
        y=randn(1)*s1 + x;
    elseif u<p1+p2
        y=randn(1)*s2 + m2;
    else
        y=x;
    end
    kkk=1;
    
function y = sample(vect)
  nvect=length(vect(:));
  y = vect(ceil( rand(1)*nvect ) );                                                                                                                             bproposalmatern.m                                                                                   0000664 0002715 0000074 00000006211 10306210545 013116  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function thetaf=bproposal(theta,fixwt,gibbs)
%% thetaf=bproposal(theta) 
%% fknots is in actual scale, not log scale

global smoovect wtvect

%randn('state',sum(100*clock))
%rand('state',sum(100*clock))


gmvect = 2*smoovect+2;  
ngm = length(smoovect);

%rangof = rrango(theta(1));
smoof = rsmoo(theta(1),ngm,smoovect,gibbs);
gmf = (smoof+1)*2;
sig2f = rsig2(theta(2),gibbs);

bcoef = theta(3:end-1);bcoeff = bcoef; 
bcoeff(2:end-1) = rbcoef(bcoef(2:end-1),gibbs);
bcoeff(1)=bcoeff(3);
bcoeff(end) = contderiv(bcoeff,gmf);

if fixwt
  wtf = theta(end);
else
  nwt = length(wtvect);
  wtf = rwt(theta(end),nwt,wtvect,gibbs);
end

thetaf = [smoof,sig2f,bcoeff, wtf];
    
%function y = rsmoo(smoo)
%    p1=1/3;p2=1/3;p3=1-p1-p2;s1=.1;s2=1;
%    y = exp(mixture(log(smoo),p1,p2,p3,s1,s2));
%    while (y > 5)
%      y = exp(mixture(log(smoo),p1,p2,p3,s1,s2));
%    end
 
function y = rsmoo(smoo,ngm,smoovect,gibbs)
  ns = find(smoovect==smoo);
  if gibbs
    p1=.9;p2=.1;p3=1-p1-p2;s1 = .1;  %% gibbs
  else
    %    p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;
    p1=.2;p2=.1;p3=1-p1-p2;s1 = .1; 
  end
  u = rand(1);
  if u<p1
    ini = max(1,ns-s1*ngm);
    fin = min(ngm, ns+s1*ngm);
    y =  sample(smoovect(ini:fin));
  elseif u < p2 + p1
    ini = 1;
    fin = ngm;
    y = sample(smoovect(ini:fin));
  else
    y = smoo;
  end

function y = rwt(wt,nwt,wtvect,gibbs)
  ns = find(wtvect==wt);
%  p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;
   if gibbs
     p1=.8;p2=.2;p3=1-p1-p2;s1 = .1;  
   else
     %p1=.3;p2=.2;p3=1-p1-p2;s1 = .1;  
%     p1=.2;p2=.1;p3=1-p1-p2;s1 = .1; % bproposa;
%     p1=.2;p2=.2;p3=1-p1-p2;s1 = .1; % bproposalmatern not good
%     p1=.1;p2=.5;p3=1-p1-p2;s1 = .1; % bproposalmatern better
%     p1=.05;p2=.3;p3=1-p1-p2;s1 = .1; % bproposalmatern similar
     p1=.5;p2=.05;p3=1-p1-p2;s1 = .1; % bproposalmatern
   end
%  p1=.7;p2=.1;p3=1-p1-p2;s1 = .1; 
  u = rand(1);
  if u<p1
    ini = max(1,ns-s1*nwt);
    fin = min(nwt, ns+s1*nwt);
    y =  sample(wtvect(ini:fin));
  elseif u < p2 + p1
    ini = 1;
    fin = nwt;
    y = sample(wtvect(ini:fin));
  else
    y = wt;
  end

function y = rsig2(x,gibbs)
    %p1=1/3;p2=1/3;p3=1-p1-p2;s1=.1;s2=1;
    if gibbs
      p1=.7;p2=.3;p3=1-p1-p2;s1=.1;s2=1;  %gibbs
    else
      %p1=.3;p2=.2;p3=1-p1-p2;s1=.1;s2=1; 
      p1=.2;p2=.1;p3=1-p1-p2;s1 = .1;s2=1; 
    end
    y = exp(mixture(log(x),p1,p2,p3,s1,s2));

function y = rbcoef(x,gibbs)    
    %p1=1/3;p2=1/3+;p3=1-p1-p2;s1=.1;s2=1;
    %p1=.5;p2=0;p3=1-p1-p2;s1=1;s2=1;
    %p1=.3;p2=.5;p3=1-p1-p2;s1=.1;s2=1;
    %p1=.35;p2=.4;p3=1-p1-p2;s1=.1;s2=1; %well with polmatern
    %p1=.4;p2=.4;p3=1-p1-p2;s1=.1;s2=1;
    if gibbs
      p1=.6;p2=.4;p3=1-p1-p2;s1=.1;s2=1; %% gibbs
    else
      %p1=.3;p2=.2;p3=1-p1-p2;s1=.1;s2=1; 
      p1=.2;p2=.1;p3=1-p1-p2;s1 = .1;s2=1; 
    end
    y = x;
    for ii=1:length(y)
        y(ii) = exp(mixture(log(x(ii)),p1,p2,p3,s1,s2));    
    end

function y=mixture(x,p1,p2,p3,s1,s2,m2)
    if nargin==6
        m2=0;
    end
    u=rand(1);
    if u<p1
        y=randn(1)*s1 + x;
    elseif u<p1+p2
        y=randn(1)*s2 + m2;
    else
        y=x;
    end
    kkk=1;
    
function y = sample(vect)
  nvect=length(vect(:));
  y = vect(ceil( rand(1)*nvect ) );
                                                                                                                                                                                                                                                                                                                                                                                       calcCova.m                                                                                          0000644 0002715 0000074 00000001246 10253636551 011434  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function cova=calcCova(theta)
%% y=calcCova(theta)

global Za dista M cnots indupper rr nt 

  %------------------------------------
  %% parse values of param
  %------------------------------------
  smoo = theta(1);
  sig2 = theta(2);
  bcoef = theta(3:end-1);
  wt = theta(end);

  covacnots = bhankelx(cnots,smoo,bcoef,wt); % bhankel forces covacnots(0)=1

  % define matrix cova
  cova = zeros(size(dista));

  % get spline coeff
  cc = spline(cnots,covacnots);
  % interpolate values of cova(rr) using spline coeff
  cova(indupper) = ppval(cc,rr);

  % complete diag and lower triangular part of cova matrix
  cova = cova + cova' + eye(size(cova));
  cova = cova * sig2;
                                                                                                                                                                                                                                                                                                                                                          calcHpol.m                                                                                          0000644 0002715 0000074 00000005022 10253635752 011444  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function Hpol = calcHpol(ll,cnots,wtvect)
%% calc hankel transform of polynomials
  nnots = length(cnots);
  nwt = length(wtvect);
  Hpol = zeros(ll,4,nnots+1,nwt);
  for jj=1:nwt
    wt = wtvect(jj);
    knots = knodos(wt,ll);
    for ii=1:nnots
      r = cnots(ii); 
      rknots = r*knots;
      Hpol(:,:,ii,jj) = mhat(rknots) .* repmat(1./[r^5  r^4 r^3 r^2],ll,1);
    end
    XXhat0 = zeros(ll,4);
    SS1 = knots(2:end)'; SS0 =knots(1:(end-1))';
    SSdif = SS1 - SS0;
    XXhat0(:,4) = SSdif    .* (SS0 +   SS1) / 2;
    XXhat0(:,3) = SSdif.^2 .* (SS0 + 2*SS1) / 6;
    XXhat0(:,2) = SSdif.^3 .* (SS0 + 3*SS1) / 12;
    XXhat0(:,1) = SSdif.^4 .* (SS0 + 4*SS1) / 20;
    Hpol(:,:,nnots+1,jj) = XXhat0;
  end
  

%%======================================
function res=mhat(rknots)
% calc matrix n*4 of hankel of monomials 1,x-a,(x-a)^2,(x-a)^3
% between knots
  n = length(rknots)-1;
  XXhat = zeros(n,4); % will hold hankel tr of all monomials
  TT = [besselj(0:1,rknots'),struveh0(rknots'),struveh1(rknots')]; %table of special functions at knots*rr
  for ik=1:n
    XXhat(ik,:)=xhat(ik,rknots,TT); %% knots(4) corresponde a 0
  end
  res = XXhat;

%%======================================
function res=xhat(k,rknots,TT)
%% calculates hankel transform of (1,x-k,(x-k)^2,(x-k)^3)I(i,i+1)
%% \int_a^b x (x-k)^j Jo(x) dx
%% k indicates the interval (a,b) where integration takes place
%% rknots = vector of knots*r; a=rknots(k), b=rknots(k+1)
  
  a = rknots(k); b = rknots(k+1);
  j0a = TT(k,1); j1a = TT(k,2); h0a = TT(k,3); h1a = TT(k,4);
  j0b = TT(k+1,1); j1b = TT(k+1,2); h0b = TT(k+1,3); h1b = TT(k+1,4);
  %%-a BesselJ[1, a] + b BesselJ[1, b]
  ki=a;
  
  xhat0 = - a * j1a + b * j1b;
  
  xhat1 = j1a*(-a^2 + a*(ki + (pi*h0a)/2.)) + ...
          j1b*(b^2 + b*(-ki - (pi*h0b)/2.)) - ...
          (a*pi*j0a*h1a)/2. + ...
          (b*pi*j0b*h1b)/2.;

  xhat2 = j1a*(-a^3 + 2*a^2*ki + a*(4 - ki^2 - ki*pi*h0a)) + ...
          j1b*(b^3 - 2*b^2*ki + b*(-4 + ki^2 + ki*pi*h0b)) + ...
          j0a*(-2*a^2 + a*ki*pi*h1a) + ...
          j0b*(2*b^2 - b*ki*pi*h1b);

  xhat3 = j1a*(-a^4 + 3*a^3*ki + a^2*(9 - 3*ki^2) + a*(-12*ki + ki^3 + ...
                                                    (-4.5 + (3*ki^2)/2.)*pi*h0a)) + ...
          j1b*(b^4 - 3*b^3*ki + b^2*(-9 + 3*ki^2) + b*(12*ki - ki^3 + (4.5 ...
                                                    - (3*ki^2)/2.)*pi*h0b)) + ...
          j0a*(-3*a^3 + 6*a^2*ki + a*(4.5 - (3*ki^2)/2.)*pi*h1a) + ... 
          j0b*(3*b^3 - 6*b^2*ki + b*(-4.5 + (3*ki^2)/2.)*pi*h1b);  
  
  res = [xhat3,xhat2,xhat1,xhat0];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                calcTailmat.m                                                                                       0000644 0002715 0000074 00000002760 10253647111 012133  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function tailmat = calcTailmat(cnots,smoovect,wtvect)
%% calc tailmat(nnots,ngm); ngm = 100
%% rwt = cnots*wt

T=50; %% if r*wt>T use asymptotic approx of tail integral
%if gm>14;%  T=40; %end   %% %%% gm is < 10
N = length(cnots);
gmvect = smoovect*2 +2;
ngm=length(gmvect);
nwt = length(wtvect);
tailmat=zeros(N,ngm,nwt);

for jj = 1:nwt

  rwt = cnots*wtvect(jj);
  
  indi = find(rwt>=T);
  nindi = length(indi);
  for gii = 1:length(gmvect) 
    if(nindi>0)
      for rwii=1:nindi
        rowii = indi(rwii);
        tailmat(rowii,gii,jj) = tailhatapp(rwt(rowii),gmvect(gii));
      end
    end
  end
  
  indi = find(rwt<T);
  nindi = length(indi);
  for gii = 1:ngm 
    if(nindi>0)
      for rwii=1:nindi
        rowii = indi(rwii);
        tailmat(rowii,gii,jj) = tailhat(rwt(rowii),gmvect(gii));
      end
    end
  end
  
end
%%======================================
function res=tailhatapp(st,gm)
%% res = tailhat(st,gm)
%% calc int(s^(1-gm) J0(s),{si,st,inf})
    pi4=pi/4;sq2pi=sqrt(2*pi);
    res= ( (-3 + 8*gm)*st^(-0.5 - gm)*cos(pi/4. - st)) / (4.*sq2pi) - ...
         ( (-15 + 16*gm + 128*gm^2)*st^(-1.5 - gm)*cos(pi/4. + st))/(64.*sq2pi) + ...
           sqrt(2/pi)*st^(0.5 - gm)*cos(pi/4. + st) ;
%%======================================
function res=tailhat(st,gm)
% calls Ken Wilder's implementation of 1F2 gh1f2.mexglx
% res = tailhat(st,gm)
% calc int(s^(1-gm) J0(s),{s,st,inf})
a = 1-gm/2;b=2-gm/2;z=-st^2/4;
res = - gamma(-gm/2)*gm/2^gm/gamma(gm/2) + ...
    st^(2-gm)* hg1f2(a,1,b,z) /(gm-2);                contderiv.m                                                                                         0000644 0002715 0000074 00000000416 10247110557 011710  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function flp1=contderiv(bcoef,gm)
%% fnp1 = contderiv
%% returns f_(n+1) such that deriv is cont at wt
%% called by bhankel.m and bproposal.m
    l=length(bcoef)-3;
    fl0 = bcoef(end-1);
    flm1 = bcoef(end-2);
    flp1 = ( flm1*(3*l - gm) - fl0*4*gm )/(3 * l + gm);
                                                                                                                                                                                                                                                  covariogram.m                                                                                       0000664 0002715 0000074 00000001675 10261065334 012234  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   ncov = length(indupper);
covario = nan(ncov,nt);
[kk,indsort] = sort(rr);
rsort = rr(indsort);

for ii=1:nt
  z = Za(:,ii); zm = mean(z);
  zdif = repmat(z,1,nobs);
  zdif = (zdif - zm).*(zdif - zm)';
  zdiflin = zdif(indupper);
  covario(:,ii) = zdiflin(indsort);
end
covario = sum(covario,2)/nt;

rohat = nadwat(cnots,covario,rsort,100);


figure(3);
%length( llvect( find(~isnan(llvect)) ) );
%evol = cumsum(~isnan(llvect));
%plot(evol); ylabel('cumulated number of changes in ener')
%xlabel('iteration number');title('evolution of energy')
[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);
%[simcovacnots,simcova0] = bhankelx(cnots,simsmoo,simbcoef,simwt);
plot(cnots,polmatern(cnots,u,v,simirango,simsmoo))
%plot(cnots,simcovacnots)
hold on; plot(cnots,covacnots,'r.')
grid
hold off
xlabel('distance');title('cov function')
legend('simulated','estimated')


figure(3);hold on
plot(rsort,covario,'.','markersize',2);
plot(cnots,rohat,'g.');hold off
                                                                   ener.m                                                                                              0000644 0002715 0000074 00000002303 10305667314 010644  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=ener(theta)
%% y=ener(theta)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,bcoef,wt] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt reml

cova = calcCova(theta);
mineig = min(eig(cova));
nobs = size(Za,1);


if( mineig > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  nobs = length(Za(:,1));
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
  else
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  end
  %%==========================================
else
  disp('WARNING: NEGATIVE EIGENVALUE')
  disp(mineig)
  disp('condition number = ')
  disp(cond(cova));
  disp('theta was ')
  disp(theta)
  y = 1e100;
end


                                                                                                                                                                                                                                                                                                                             enermatern.m                                                                                        0000664 0002715 0000074 00000002413 10314044660 012050  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=enermatern(theta)
%% y=ener(theta)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,bcoef,wt] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt reml

smoo = exp(theta(1));
sig2 = exp(theta(2));
rango = exp(theta(3));

%disp([smoo sig2 2*sqrt(smoo)/rango]) %% DEBUG

cova = hmatern(dista,rango,smoo)*sig2;
mineig = min(eig(cova));


if( mineig > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  nobs = length(Za(:,1));
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
  else
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  end

else
  disp('WARNING: NEGATIVE EIGENVALUE')
  disp(mineig)
  disp('condition number = ')
  disp(cond(cova));
  disp('theta was ')
  disp(theta)
  y = 1e100;
end


                                                                                                                                                                                                                                                     enerpolmatern.m                                                                                     0000664 0002715 0000074 00000001776 10260263430 012574  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=ener(theta)
%% y=ener(theta)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,bcoef,wt] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt 

smoo = exp(theta(1));
sig2 = exp(theta(2));
rango = exp(theta(3));
a = 2*sqrt(smoo)/rango;
u = exp(theta(4));
v = exp(theta(5));

cova = polmatern(dista,u,v,a,smoo)*sig2;
mineig = min(eig(cova));

if( mineig > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  xx = invcc' * Za;
  nobs = size(Za,1);
  xx=reshape(xx,nt*nobs,1);
  y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  %%==========================================
else
  disp('WARNING: NEGATIVE EIGENVALUE')
  disp(mineig)
  disp('condition number = ')
  disp(cond(cova));
  disp('theta was ')
  disp(theta)
  y = 1e100;
end


  eucldist.m                                                                                          0000664 0002715 0000074 00000000673 10313100344 011521  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=eucldist(loc1,loc2)
%% calculates cord distance bw two sets of points given x and y
%% loc1=[x1;y1] (n1x2), loc2=[x2,y2] (n2x2)
  if nargin==1
    loc2 = loc1;
  end
  x1 = loc1(:,1);x2 = loc2(:,1);
  y1 = loc1(:,2);y2 = loc2(:,2);
  n1 = length(x1); n2 = length(x2);
  y = (repmat(x1,1,n2)-repmat(x2',n1,1)).^2 + (repmat(y1,1,n2)- ...
                                               repmat(y2',n1,1)).^2 ;
  y = sqrt(y)*6400*pi/180;                                                                     fmatern.m                                                                                           0000644 0002715 0000074 00000000332 10106716134 011341  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=fmatern(x,irange,smoo)
%% calculates the matern spectral density 1/(a^2+w^2)^(nu+1/2)
    d=2;
    y = gamma(smoo+d/2)/gamma(smoo) * irange^(2*smoo) / pi^(d/2) ...
        ./ (irange^2 + x.^2 ).^(smoo+d/2);
                                                                                                                                                                                                                                                                                                      fpolmatern.m                                                                                        0000664 0002715 0000074 00000000401 10261313621 012047  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = fpolmatern(w,u,v,a,nu)
p0=(2*a^4 + 2*a^2*nu*(-u^2 + v^2) + nu*(1 + nu)*(u^2 + v^2)^2) / ...
    (2.*a^(2*(2 + nu))*nu*(1 + nu)*(2 + nu));

res = ( (w-u).^2 + v^2 ).*( (w+u).^2 + v^2) ./ ...
        ( w.^2 + a^2).^(nu+1+2);

res = res/p0/2/pi;                                                                                                                                                                                                                                                               gen_tables.m                                                                                        0000664 0002715 0000074 00000001056 10371744064 012025  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% generate tables for hankel transform of S+T
%% ------------------------------------
function res = gen_tables(ll,rmin,rmax,nnots,nwt)

wtvect = linspace(1/rmax,1/rmin,nwt);
cnots = linspace(rmin,rmax,nnots);

Hpol = calcHpol(ll,cnots,wtvect);

save(strcat('raintablahpol_',num2str(ll)),'Hpol')

% smoothness values for tail tabulation
%ngm = 100;
%smoovect = linspace(.05,5,ngm)+.0001;
%gmvect = smoovect*2 + 2;
%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  gibbsbsplines.m                                                                                     0000664 0002715 0000074 00000025300 10301451470 012533  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values
global covario rsort rmax;

%------------------------------------
%% read sites locations
%% simul.dat: "lat"   "wnh4"  "lon"   "cwnh4" "site"
%------------------------------------
[lat,wnh4,lon,cwnh4]=textread('simul.dat','%f%f%f%f%*[^\n]');
%M = load('M.txt');
nobs = length(wnh4);
M = ones(nobs,1);
dista = cordist([lon,lat]);
indiag=find(eye(nobs)); dista(indiag)=0; %% clean up the diag
indupper = find(triu(dista,1)); 
rr = dista(indupper); rmin = min(rr); rmax = max(rr);
nnots=100; cnots = linspace(rmin,rmax,nnots);

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

ll = 3 % number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

%% calculate or read tabulated values
%% hankel transform of polynomials
%tic
%Hpol = calcHpol(ll,cnots,wtvect);
%disp('hPol calc time')
%toc
%save(strcat('tablahpol_',num2str(ll)),'Hpol')

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


%testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
testsplinetail = false; testmatern = true; testpolmatern = false;
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
%%randn('state',0);
%------------------------------------
% generate observations
%------------------------------------
%indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
%simgm=2*(simsmoo+1);
%simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end);
%simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
%simsig2=1;
%simll = 4;

%if testsplinetail
%  simt = nodos(simwt,simll);
%  simknots = simt(4:end-3);
%  %  simbcoef = [0.8 1.0 .8 .2  1.2 1 .8 .6 .4]; %% ll=6 should have ll+3 elements
%  simbcoef = [0.2 1 .2  2  .6 .4 .2]; %% ll=4
%  simbcoef(end)=contderiv(simbcoef,simgm);
%  simtheta = [simsmoo simsig2 simbcoef simwt];
%  simcova = calcCova(simtheta);
%elseif testpolmatern 
%  % test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
%  u = 0.5*simirango; v=.1*simirango;
%  simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);
%elseif testmatern
%  % test matern cova
%  simcova = simsig2*hmatern(dista,simrango,simsmoo);
%end

%for ii=1:nt
%  Za(:,ii) = simulachol(simcova);
%end

%if testsplinetail
%  save(strcat('splinetail_',num2str(nt),'_',num2str(ll),'.mat'),...
%       'Za','simbcoef','simll','simsig2','simsmoo','simgm','simwt',...
%       'simt','simknots', 'simtheta','simcova')
%elseif testpolmatern 
%  save(strcat('polmatern_',num2str(nt),'_','.mat'),...
%       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango','u','v',...
%       'simcova')
%elseif testmatern
%  save(strcat('matern_',num2str(nt),'.mat'),...
%       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango',...
%       'simcova')
%end

%% read observations
if testsplinetail
  simll = 4;
  load(strcat('splinetail_',num2str(nt),'_',num2str(simll),'.mat'));
elseif testpolmatern 
  load(strcat('polmatern_',num2str(nt),'_','.mat'))
elseif testmatern
  load(strcat('matern_',num2str(nt),'.mat'))
end

cc = chol(simcova);
invcc = inv(cc);
xx = invcc' * Za;
nobs = size(Za,1);
xx=reshape(xx,nt*nobs,1);
truene = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
clear('xx','cc');

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=0.2); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=3;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= .02); wt = wt(end);

t = nodos(wt,ll);
theta = [smoo,sig2,bcoef,wt];thetaini=theta;
Ei = ener(theta); Eini = Ei;


%------------------------------------
%% improve initial guess fitting matern
%------------------------------------
matheta = log([smoo sig2 1000]);
%[paramf fval exitflag] = fminunc(@enermatern,matheta)
[paramf fval exitflag] = fminsearch(@enermatern,matheta);
matheta = exp(paramf);
matsmoo = smoovect(smoovect<=matheta(1));matsmoo=matsmoo(end);
matsig2 = matheta(2);
matirango = 2*sqrt(matsmoo)/matheta(3);
matwt = wtvect(wtvect <= matirango); matwt = matwt(end);

theta = [matsmoo matsig2 bcoef matwt];
%% DEBUG

testwt = wtvect(wtvect <= 1.5*simwt); testwt = testwt(end);
theta = [matsmoo matsig2 bcoef testwt];


matEi = ener(theta); Ei = matEi;

%------------------------------------
%% loop for simulated annealing
%------------------------------------
%% minimize (-l(f)+penalty)
L = 3000 %%max number of iterations
Ti = 1000; %% scale temperature
lam = 0.99;
%tsca = Ti * lam.^(1:L);


tsca = ones(L,1);
a=300;
for kk=1:L;
  tsca(kk+1) = tsca(kk)/(1+a*tsca(kk));
end
tsca = tsca/tsca(1) * Ti;
exp(-.01/tsca(end))/(1+ exp(-.01/tsca(end)) )
plot(tsca(100:end))

indtheta = [1:2 4:length(theta)];indtheta(end-1)=[]; %% bcoef(1) and
                                                     %bcoef(end) not changed
np = length(indtheta);

fixwt = true;
conta=1;
llvect = nan((L+1)*np+1,1);
llvect(1) = Ei;
thetavect = nan((L+1)*np+1,length(theta));
thetavect(1,:)=theta;
gibbs = false;


for lind = 1:L
  T=tsca(lind);
  % generate new state %logf m-vector
  for i=1:np
    %% proposal theta
    if gibbs
      newtheta = bproposal(theta,fixwt,gibbs);  %% bproposal.m generates f and proposal.m log(f)
      thetaf = theta; thetaf(indtheta(i)) = newtheta(indtheta(i));
      thetaf(3)=thetaf(5);thetaf(end-1)=contderiv(thetaf(3:end-1),thetaf(1)*2+2);
    else
      thetaf = bproposal(theta,fixwt,gibbs);
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
    thetavect((lind-1)*np+i+1,:)=theta;
    
  end
end
thetaf = theta; %% revert back to last updated values
Ef = Ei; %% revert back to last updated values

smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll);

disp('number of indep replicates')
disp(nt)
disp('sim  ini  after-matern  final')
disp([simsmoo smoo matsmoo smoof ])
disp([simsig2 sig2 matsig2 sig2f])
disp([simwt wt matwt wtf])
disp(round([truene Eini matEi Ei]))
disp('estimated ener - true ener'); disp(Ei - truene)
disp('iter');disp(L)
disp('lam -- cooling schedule T=lam^(1:iter))');disp(lam)


[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);

%% calc covariogram cloud and nadwat estimate of cov fn
dindupper = find(triu(ones(size(dista))));

ncov = length(dindupper);
covario = nan(ncov,nt);
[rsort,indsort] = sort(dista(dindupper));

zm = mean(Za(:));
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
xxx = linspace(0,2*simwt,100);
fcovario = xxx(2:end); 
for ii = 2:length(xxx)
  fcovario(ii-1) = hankel0('kcovariox',xxx(ii));
end

if testsplinetail
  load(strcat('tablahpol_',num2str(simll)))  
  [simcovacnots,simcova0] = bhankelx(cnots,simsmoo,simbcoef,simwt);
  simspectral = spectralfn(xxx,simsmoo,simt,simbcoef)/simcova0;
  prefijo = strcat(num2str(nt),'_',num2str(simll),'_',...
                    num2str(ll));
  epsfile = strcat('splinetail_',prefijo,'.eps');
  savefile = strcat('opt_st_',prefijo,'.mat');
elseif testpolmatern 
  simcovacnots = polmatern(cnots,u,v,simirango,simsmoo);
  simspectral = fpolmatern(xxx,u,v,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_', num2str(ll));
  epsfile = strcat('polmatern_',prefijo,'.eps');
  savefile = strcat('opt_polmat_',prefijo,'.mat');
elseif testmatern
  simcovacnots = hmatern(cnots,simrango,simsmoo);
  simspectral = fmatern(xxx,simirango,simsmoo);
  prefijo = strcat(num2str(nt),'_',num2str(ll));
  epsfile = strcat('matern_',prefijo,'.eps');
  savefile = strcat('opt_matern_',prefijo,'.mat');
end

mlespect = spectralfn(xxx,smoof,tf,bcoeff)/cova0;
rmsecovak = rmse(rohat,simcovacnots);
rmsecovam = rmse(covacnots,simcovacnots);
rmsespeck = rmse(fcovario/2/pi,simspectral(2:end));
rmsespecm = rmse(mlespect(2:end),simspectral(2:end));

%% calc likelihood using rohat
kcova = kcovario(dista);
cc = chol(kcova);
invcc = inv(cc);
xx = invcc' * Za;
xx=reshape(xx,nt*nobs,1);
klike = -( nt*sum(log(diag(cc))) + 0.5 * xx'*xx );

%% clear some variables and save
clear('xx','cc');
clear('Hpol','tailmat')
save(savefile);


%% plot results

fighandle = figure(1);
set(fighandle,'PaperPosition',[0.25 2.5 8 3.5]);
subplot(1,2,1)
plot(cnots,simcovacnots)
hold on; plot(cnots,covacnots*sig2f,'r.')
grid
%plot(rsort,covario,'.','markersize',2);
plot(cnots,nadwat(cnots,covario,rsort,50),'k+');hold off

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
text(rx+dx,ry+dy,sprintf('%4.0f',-Ei))
text(rx+2*dx,ry+dy,sprintf('%4.0f',klike))

text(rx,ry+2*dy,sprintf('%4.3f',0))
text(rx+dx,ry+2*dy,sprintf('%4.3f',rmsecovam))
text(rx+2*dx,ry+2*dy,sprintf('%4.3f',rmsecovak))

text(rx,ry+3*dy,sprintf('%3.2f',simsmoo))
text(rx+dx,ry+3*dy,sprintf('%3.2f',smoof))
text(rx+2*dx,ry+3*dy,'NA')

text(rx,ry+4*dy,sprintf('%3.2f',simsig2))
text(rx+dx,ry+4*dy,sprintf('%3.2f',sig2f))
text(rx+2*dx,ry+4*dy,sprintf('%3.2f',kcovario(0)) )

text(rx,ry+5*dy,sprintf('%5.4f',simwt) )
text(rx+dx,ry+5*dy,sprintf('%5.4f',wtf) )
text(rx+2*dx,ry+5*dy,'NA')

legend('simulated','estimated','kernel estimator')
xlabel('distance');title('cov function')

subplot(1,2,2)
plot(xxx,simspectral);hold on
plot(xxx,mlespect,'r.')
plot(xxx(2:end),fcovario/2/pi,'k+');hold off
legend('simulated','estimated','kernel estimator')
grid
xlabel('freq')
ylabel('f(w)')
title('Spectral density')
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
                                                                                                                                                                                                                                                                                                                                hankel0.m                                                                                           0000645 0002715 0000074 00000127646 10274764156 011267  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %
%  [z0,z1]=hankel01(funname,B)
%
%    funname       The name of the function to be transformed.
%    B             The transform argument.
%    z0            The value of the transform of order 0.
%    z1            The value of the transform of order 1.
%
%  Based on Walt Anderson's Fortran code, which was published in:
%
%  Anderson, W. L., 1979, Computer Program Numerical Integration of Related
%  Hankel Transforms of Orders 0 and 1 by Adaptive Digital Filtering.
%  Geophysic, 44(7):1287-1305.
%
function z0=hankel0(funname,B)
%
%  Points at which the function will be evaluated.  Note that these
%  are scaled by B before being used.
%
YBASE=[
   8.9170998013274418e-14
   9.8549193740052245e-14
   1.0891370292130841e-13
   1.2036825704856076e-13
   1.3302749714952345e-13
   1.4701812115404443e-13
   1.6248015192957209e-13
   1.7956833867707590e-13
   1.9845370571306282e-13
   2.1932526413842005e-13
   2.4239190352504162e-13
   2.6788448255287407e-13
   2.9605813952117967e-13
   3.2719484585839035e-13
   3.6160622818693727e-13
   3.9963668718722961e-13
   4.4166684447542096e-13
   4.8811735199247527e-13
   5.3945310203017795e-13
   5.9618788002944785e-13
   6.5888950671771899e-13
   7.2818552104963217e-13
   8.0476946082781583e-13
   8.8940780386232124e-13
   9.8294763913816720e-13
   1.0863251447666186e-12
   1.2005749575703847e-12
   1.3268405280766938e-12
   1.4663855645544968e-12
   1.6206066806315702e-12
   1.7910473730731204e-12
   1.9794134696161976e-12
   2.1875902014670364e-12
   2.4176610713286156e-12
   2.6719287057960000e-12
   2.9529379008172425e-12
   3.2635010908665675e-12
   3.6067264967338821e-12
   3.9860492336431492e-12
   4.4052656910401313e-12
   4.8685715281339743e-12
   5.3806036654647833e-12
   5.9464866927629101e-12
   6.5718841575654077e-12
   7.2630552479033658e-12
   8.0269174363595134e-12
   8.8711157124588670e-12
   9.8040990962934688e-12
   1.0835205199155280e-11
   1.1974753677488471e-11
   1.3234149515479672e-11
   1.4625997169973058e-11
   1.6164226720110948e-11
   1.7864233284247930e-11
   1.9743031099469828e-11
   2.1819423805797136e-11
   2.4114192639334462e-11
   2.6650304417866294e-11
   2.9453141400488783e-11
   3.2550755321790057e-11
   3.5974148143038491e-11
   3.9757582330231206e-11
   4.3938923764369768e-11
   4.8560020715924429e-11
   5.3667122676390673e-11
   5.9311343238745088e-11
   6.5549171659463761e-11
   7.2443038221987792e-11
   8.0061939059983487e-11
   8.8482126693838504e-11
   9.7787873191515278e-11
   1.0807231359173195e-10
   1.1943837803073370e-10
   1.3199982190169224e-10
   1.4588236435691521e-10
   1.6122494654737812e-10
   1.7818112219246310e-10
   1.9692059439719362e-10
   2.1763091409794871e-10
   2.4051935713527238e-10
   2.6581499874015355e-10
   2.9377100619593265e-10
   3.2466717262156566e-10
   3.5881271723520049e-10
   3.9654938012404428e-10
   4.3825484249401899e-10
   4.8434650663021337e-10
   5.3528567339924565e-10
   5.9158215910338558e-10
   6.5379939789346253e-10
   7.2256008080722362e-10
   7.9855238787053352e-10
   8.8253687563437822e-10
   9.7535408908045947e-10
   1.0779329740778886e-09
   1.1913001745856734e-09
   1.3165903076505279e-09
   1.4550573190356333e-09
   1.6080870331313015e-09
   1.7772110227512649e-09
   1.9641219376281762e-09
   2.1706904450210516e-09
   2.3989839519819520e-09
   2.6512872966606394e-09
   2.9301256157327411e-09
   3.2382896168163261e-09
   3.5788635088117363e-09
   3.9552558697009004e-09
   4.3712337607414382e-09
   4.8309604284818817e-09
   5.3390369719324459e-09
   5.9005483919104070e-09
   6.5211144834374116e-09
   7.2069460805369272e-09
   7.9649072163486862e-09
   8.8025838206794284e-09
   9.7283596425381269e-09
   1.0751500157513941e-08
   1.1882245299770154e-08
   1.3131911946747030e-08
   1.4513007182274982e-08
   1.6039353471673310e-08
   1.7726227001629019e-08
   1.9590510569407677e-08
   2.1650862551562963e-08
   2.3927903643240499e-08
   2.6444423237025737e-08
   2.9225607506844726e-08
   3.2299291479658129e-08
   3.5696237617766722e-08
   3.9450443699873716e-08
   4.3599483082281089e-08
   4.8184880745668262e-08
   5.3252528891055793e-08
   5.8853146244378082e-08
   6.5042785666539686e-08
   7.1883395149287240e-08
   7.9443437811532343e-08
   8.7798577101256825e-08
   9.7032434060731539e-08
   1.0723742423401342e-07
   1.1851568259277232e-07
   1.3098008573741623e-07
   1.4475538160404734e-07
   1.5997943798373573e-07
   1.7680462234971138e-07
   1.9539932680224871e-07
   2.1594965339340473e-07
   2.3866127669890702e-07
   2.6376150227843728e-07
   2.9150154162607256e-07
   3.2215902637935327e-07
   3.5604078695002664e-07
   3.9348592338593700e-07
   4.3486919919827997e-07
   4.8060479212078480e-07
   5.3115043933968360e-07
   5.8701201868132178e-07
   6.4874861160747570e-07
   7.1697809869053571e-07
   7.9238334356995170e-07
   8.7571902728105488e-07
   9.6781920135651654e-07
   1.0696056352944216e-06
   1.1820970419372223e-06
   1.3064192730922673e-06
   1.4438165874351013e-06
   1.5956641034684996e-06
   1.7634815621706371e-06
   1.9489485370736003e-06
   2.1539212439998211e-06
   2.3804511186939235e-06
   2.6308053482811660e-06
   2.9074895620382204e-06
   3.2132729085731429e-06
   3.5512157703953872e-06
   3.9247003932525885e-06
   4.3374647367828189e-06
   4.7936398852710157e-06
   5.2977913929290109e-06
   5.8549649774966190e-06
   6.4707370194807025e-06
   7.1512703724455683e-06
   7.9033760429228480e-06
   8.7345813572541237e-06
   9.6532052976029776e-06
   1.0668441761124589e-05
   1.1790451575578641e-05
   1.3030464192308714e-05
   1.4400890074365673e-05
   1.5915444904593194e-05
   1.7589286856791647e-05
   1.9439168303816350e-05
   2.1483603480955747e-05
   2.3743053782621043e-05
   2.6240132546858776e-05
   2.8999831377238596e-05
   3.2049770267221755e-05
   3.5420474030339066e-05
   3.9145677802784465e-05
   4.3262664675996812e-05
   4.7812638838370288e-05
   5.2841137960621061e-05
   5.8398488952101537e-05
   6.4540311649424626e-05
   7.1328075478483034e-05
   7.8829714661124184e-05
   8.7120308123675967e-05
   9.6282830912076270e-05
   1.0640898463402168e-04
   1.1760011523947924e-04
   1.2996822732501724e-04
   1.4363710511345378e-04
   1.5874355132796403e-04
   1.7543875635971472e-04
   1.9388981143211579e-04
   2.1428138090594562e-04
   2.3681755046234149e-04
   2.6172386966089197e-04
   2.8924960931543914e-04
   3.1967025628016629e-04
   3.5329027061462893e-04
   3.9044613272236350e-04
   4.3150971095986065e-04
   4.7689198342006658e-04
   5.2704715113927156e-04
   5.8247718389374341e-04
   6.4373684408196633e-04
   7.1143923897318674e-04
   7.8626195689103691e-04
   8.6895384874522254e-04
   9.6034252278312515e-04
   1.0613426275713101e-03
   1.1729650061058051e-03
   1.2963268126685603e-03
   1.4326626936829910e-03
   1.5833371444703617e-03
   1.7498581655775839e-03
   1.9338923553535471e-03
   2.1372815898255564e-03
   2.3620614568136901e-03
   2.6104816287778878e-03
   2.8850283782960702e-03
   3.1884494615157647e-03
   3.5237816186211822e-03
   3.8943809665496639e-03
   4.3039565881380203e-03
   4.7566076538702283e-03
   5.2568644477534125e-03
   5.8097337079228714e-03
   6.4207487357601564e-03
   7.0960247750331065e-03
   7.8423202153108801e-03
   8.6671042321983371e-03
   9.5786315413559676e-03
   1.0586025014468731e-02
   1.1699366984012178e-02
   1.2929800150624660e-02
   1.4289639103000504e-02
   1.5792493566432742e-02
   1.7453404613518231e-02
   1.9288995200267688e-02
   2.1317636534236600e-02
   2.3559631939745231e-02
   2.6037420060372587e-02
   2.8775799432443259e-02
   3.1802176677114019e-02
   3.5146840795050052e-02
   3.8843266308924096e-02
   4.2928448287690518e-02
   4.7443272605669898e-02
   5.2432925142121424e-02
   5.7947344016710048e-02
   6.4041719386992838e-02
   7.0777045810065886e-02
   7.8220732696592687e-02
   8.6447278966843177e-02
   9.5539018660927705e-02
   1.0558694496554391e-01
   1.1669162090437306e-01
   1.2896418580662142e-01
   1.4252746762678220e-01
   1.5751721224808804e-01
   1.7408344207293611e-01
   1.9239195749751564e-01
   2.1262599629790035e-01
   2.3498806753529980e-01
   2.5970197833480957e-01
   2.8701507382234348e-01
   3.1720071263778915e-01
   3.5056100280015512e-01
   3.8742982530616715e-01
   4.2817617572350458e-01
   4.7320785722246539e-01
   5.2297556200716211e-01
   5.7797738199458315e-01
   6.3876379388591276e-01
   7.0594316852237804e-01
   7.8018785966510817e-01
   8.6224093313756223e-01
   9.5292360367804285e-01
   1.0531434539328173e+00
   1.1639035178482904e+00
   1.2863123193718711e+00
   1.4215949669322265e+00
   1.5711054147362089e+00
   1.7363400135976372e+00
   1.9189524869191834e+00
   2.1207704817120212e+00
   2.3438138603014083e+00
   2.5903149157877352e+00
   2.8627407135861755e+00
   3.1638177826465683e+00
   3.4965594034715681e+00
   3.8642957660407120e+00
   4.2707072994710522e+00
   4.7198615069887930e+00
   5.2162536748687147e+00
   5.7648518627701284e+00
   6.3711466257477705e+00
   7.0412059655722290e+00
   7.7817360613311877e+00
   8.6001483871237632e+00
   9.5046338885843706e+00
   1.0504244960619703e+01
   1.1608986046819572e+01
   1.2829913767290970e+01
   1.4179247577028352e+01
   1.5670492062326328e+01
   1.7318572099218340e+01
   1.9139982226652432e+01
   2.1152951729381048e+01
   2.3377627082769912e+01
   2.5836273585494951e+01
   2.8553498198135060e+01
   3.1556495817904278e+01
   3.4875321454323611e+01
   3.8543191029858157e+01
   4.2596813816033411e+01
   4.7076759832163077e+01
   5.2027865883738443e+01
   5.7499684304247886e+01
   6.3546978891585546e+01
   7.0230273002547406e+01
   7.7616455290928698e+01
   8.5779449151653139e+01
   9.4800952570955843e+01
   1.0477125578728921e+02
   1.1579014494637693e+02
   1.2796790079449971e+02
   1.4142640240527066e+02
   1.5630034698636896e+02
   1.7273859797446769e+02
   1.9090567491054267e+02
   2.1098340000673559e+02
   2.3317271788416559e+02
   2.5769570669423729e+02
   2.8479780075142304e+02
   3.1475024692237560e+02
   3.4785281935573863e+02
   3.8443681972258412e+02
   4.2486839299489054e+02
   4.6955219194748827e+02
   5.1893542705903837e+02
   5.7351234234481581e+02
   6.3382916191693528e+02
   7.0048955677885772e+02
   7.7416068656769369e+02
   8.5557987671209173e+02
   9.4556199783295187e+02
   1.0450076212424869e+03
   1.1549120321646080e+03
   1.2763751908839718e+03
   1.4106127415182191e+03
   1.5589681785928965e+03
   1.7229262931862318e+03
   1.9041280332173003e+03
   2.1043869266043412e+03
   2.3257072316617105e+03
   2.5703039963907459e+03
   2.8406252274246667e+03
   3.1393763905017645e+03
   3.4695474876758481e+03
   3.8344429822617740e+03
   4.2377148710149695e+03
   4.6833992345424385e+03
   5.1759566317540530e+03
   5.7203167426353639e+03
   6.3219277061418234e+03
   6.9868106470046323e+03
   7.7216199371708199e+03
   8.5337097949943000e+03
   9.4312078887249972e+03
   1.0423096680944496e+04
   1.1519303328070666e+04
   1.2730799034675721e+04
   1.4069708856989137e+04
   1.5549433054535755e+04
   1.7184781204437102e+04
   1.8992120420636886e+04
   2.0989539161478522e+04
   2.3197028265075980e+04
   2.5636681024340771e+04
   2.8332914304083228e+04
   3.1312712913202311e+04
   3.4605899677722984e+04
   3.8245433917662871e+04
   4.2267741314984989e+04
   4.6713078474065944e+04
   5.1625935823323234e+04
   5.7055482890376610e+04
   6.3056060407206925e+04
   6.9687724170466376e+04
   7.7016846100076829e+04
   8.5116778511712779e+04
   9.4068588251431182e+04
   1.0396186803991429e+05
   1.1489563314653141e+05
   1.2697931236743495e+05
   1.4033384322573253e+05
   1.5509288235486683e+05
   1.7140414317912661e+05
   1.8943087427924512e+05
   2.0935349323906592e+05
   2.3137139232536239e+05
   2.5570493407266162e+05
   2.8259765674555639e+05
   3.1231871175151330e+05
   3.4516555739862355e+05
   3.8146693595832947e+05
   4.2158616382857127e+05
   4.6592476772641251e+05
   5.1492650330238225e+05
   5.6908179639617680e+05
   6.2893265138330159e+05
   6.9507807573703467e+05
   7.6818007509655319e+05
   8.4897027884187771e+05
   9.3825726248661662e+05
   1.0369346401734781e+06
   1.1459900082649642e+06
   1.2665148295397097e+06
   1.3997153569188234e+06
   1.5469247060505589e+06
   1.7096161975797976e+06
   1.8894181026362628e+06
   2.0881299391192668e+06
   2.3077404818776865e+06
   2.5504476670371005e+06
   2.8186805896832864e+06
   3.1151238150622859e+06
   3.4427442466117009e+06
   3.8048208197275074e+06
   4.2049773184515880e+06
   4.6472186435204167e+06
   5.1359708947577253e+06
   5.6761256689692009e+06
   6.2730890166874416e+06
   6.9328355477427216e+06
   7.6619682271663100e+06
   8.4677844598838333e+06
   9.3583491255965196e+06
   1.0342575294807941e+07
   1.1430313433829404e+07
   1.2632449991557652e+07
   1.3961016354714479e+07
   1.5429309262008933e+07
   1.7052023882367507e+07
   1.8845400889123969e+07
   2.0827389002136763e+07
   2.3017824624610133e+07
   2.5438630372484632e+07
   2.8114034483345896e+07
   3.1070813300769802e+07
   3.4338559260968812e+07
   3.7949977063839935e+07
   4.1941210992593758e+07
   4.6352206657889292e+07
   5.1227110786931656e+07
   5.6614713058756173e+07
   6.2568934407734923e+07
   6.9149366682411388e+07
   7.6421869060750201e+07
   8.4459227190926239e+07
   9.3341881654555663e+07
   1.0315873304307374e+08
   1.1400803170473446e+08
   1.2599836106711893e+08
   1.3924972437657478e+08
   1.5389474573104006e+08
   1.7007999742659190e+08
   1.8796746690225038e+08
   2.0773617796471399e+08
   2.2958398251878911e+08
   2.5372954073575363e+08
   2.8041450947784531e+08
   3.0990596088136274e+08
   3.4249905530437142e+08
   3.7851999539077419e+08
   4.1832929081601185e+08
   4.6232536638906646e+08
   5.1094854962186480e+08
   5.6468547767501700e+08
   6.2407396778608418e+08
   6.8970839992525887e+08
   7.6224566554988432e+08
   8.4241174199494874e+08
   9.3100895829826319e+08
   1.0289240251791439e+09
   1.1371369095373254e+09
   1.2567306422910707e+09
   1.3889021577146211e+09
   1.5349742727587159e+09
   1.6964089262472496e+09
   1.8748218104523966e+09
   2.0719985414859231e+09
   2.2899125303454008e+09
   2.5307447334747562e+09
   2.7969054805094066e+09
   3.0910585976653914e+09
   3.4161480682074847e+09
   3.7754274968232164e+09
   4.1724926727921586e+09
   4.6113175578536234e+09
   5.0962940589514427e+09
   5.6322759839148350e+09
   6.2246276199985800e+09
   6.8792774214728651e+09
   7.6027773435862408e+09
   8.4023684167359400e+09
   9.2860532171338844e+09
   1.0262675959279177e+10
   1.1342011011829447e+10
   1.2534860722767656e+10
   1.3853163532931507e+10
   1.5310113459941998e+10
   1.6920292148366428e+10
   1.8699814807718300e+10
   2.0666491498890625e+10
   2.2840005383231522e+10
   2.5242109718238716e+10
   2.7896845571472111e+10
   3.0830782431638401e+10
   3.4073284124964363e+10
   3.7656802698239258e+10
   4.1617203209806610e+10
   4.5994122679122765e+10
   5.0831366787370079e+10
   5.6177348299437775e+10
   6.2085571595145073e+10
   6.8615168159057838e+10
   7.5831488388260895e+10
   8.3806755641097107e+10
   9.2620789072812759e+10
   1.0236180249249139e+11
   1.1312728723650484e+11
   1.2502498789457556e+11
   1.3817398065384482e+11
   1.5270586505337646e+11
   1.6876608107657605e+11
   1.8651536476342874e+11
   2.0613135691081287e+11
   2.2781038096130206e+11
   2.5176940787416525e+11
   2.7824822764365344e+11
   3.0751184919785828e+11
   3.3985315269713715e+11
   3.7559582077719836e+11
   4.1509757807371277e+11
   4.5875377145070300e+11
   5.0700132676483929e+11
   5.6032312176626892e+11
   6.1925281890144031e+11
   6.8438020638623755e+11
   7.5635710100467944e+11
   8.3590387171037695e+11
   9.2381664932114575e+11
   1.0209752944638193e+12
   1.1283522035151340e+12
   1.2470220406715007e+12
   1.3781724935494902e+12
   1.5231161599626948e+12
   1.6833036848418264e+12
   1.8603382787767620e+12
   2.0559917634869844e+12
   2.2722223048088804e+12
   2.5111940106775947e+12
   2.7752985902466255e+12
   3.0671792909169141e+12
   3.3897573528452603e+12
   3.7462612456976733e+12
   4.1402589802589175e+12
   4.5756938182836924e+12
   5.0569237379856543e+12
   5.5887650501481416e+12
   6.1765406013813154e+12
   6.8261330469601016e+12
   7.5440437264154141e+12
   8.3374577311253535e+12
   9.2143158151247129e+12
   1.0183393868840340e+13
   1.1254390751152201e+13
   1.2438025358832957e+13
   1.3746143904869607e+13
   1.5191838479344713e+13
   1.6789578079474348e+13
   1.8555353420195434e+13
   2.0506836974615496e+13
   2.2663559846063445e+13
   2.5047107241936324e+13
   2.7681334505709973e+13
   3.0592605869234598e+13
   3.3810058314828449e+13
   3.7365893187990141e+13
   4.1295698479287656e+13
   4.5638805000929469e+13
   5.0438680022752680e+13
   5.5743362307269414e+13
   6.1605942897748391e+13
   6.8085096471220516e+13
   7.5245668574367812e+13
   8.3159324619549984e+13
   9.1905267136338875e+13
   1.0157102845705528e+14
   1.1225334676977153e+14
   1.2405913430661245e+14
   1.3710654735730897e+14
   1.5152616881705944e+14
   1.6746231510403516e+14
   1.8507448052659994e+14
   2.0453893355595603e+14
   2.2605048098024984e+14
   2.4982441759638447e+14
   2.7609868095271022e+14
   3.0513623270798212e+14
   3.3722769044002506e+14
   3.7269423624413281e+14
   4.1189083123143062e+14
   4.5520976809898188e+14
   5.0308459732695450e+14
   5.5599446629754788e+14
   6.1446891476304075e+14
   6.7909317465761662e+14
   7.5051402729526438e+14
   8.2944627657455900e+14
   9.1667990297633300e+14
   1.0130879699538496e+15
   1.1196353618452902e+15
   1.2373884407605195e+15
   1.3675257190914975e+15
   1.5113496544604105e+15
   1.6702996851533248e+15
   1.8459666365023652e+15
   2.0401086424003345e+15
   2.2546687412956410e+15
   2.4917943227741685e+15
   2.7538586193560145e+15
   3.0434844586042220e+15
   3.3635705132645935e+15
   3.7173203121568085e+15
   4.1082743021675935e+15
   4.5403452822331500e+15
   5.0178575639460460e+15
   5.5455902507190850e+15
   6.1288250686585730e+15
   6.7733992278544400e+15
   7.4857638431407750e+15
   8.2730484990213790e+15
   9.1431326049478160e+15
   1.0104724255097566e+16
   1.1167447381907442e+16
   1.2341938075624136e+16
   1.3639951033870320e+16
   1.5074477206609342e+16
   1.6659873813938872e+16
   1.8412008037975264e+16
   2.0348415826945328e+16
   2.2488477400850208e+16
   2.4853611215221080e+16
   2.7467488324221096e+16
   3.0356269288511564e+16
   3.3548865998935916e+16
   3.7077231036440888e+16
   4.0976677464246272e+16
   4.5286232252850760e+16
   5.0049026875070080e+16
   5.5312728980313968e+16
   6.1130019468443072e+16
   6.7559119737921448e+16
   7.4664374385141264e+16
   8.2516895186770448e+16
   9.1195272810315088e+16
   1.0078636337593507e+17
   1.1138615774168800e+17
   1.2310074221230024e+17
   1.3604736028656150e+17
   1.5035558606966758e+17
   1.6616862109441658e+17
   1.8364472753028080e+17
   2.0295881212439261e+17
   2.2430417672705789e+17
   2.4789445292164490e+17
   2.7396574012127472e+17
   3.0277896853110349e+17
   3.3462251062551731e+17
   3.6981506727678112e+17
   4.0870885742048762e+17
   4.5169314318104928e+17
   4.9919812573787520e+17
   5.5169925092337018e+17
   6.0972196764462810e+17
   6.7384698675270400e+17
   7.4471609299199475e+17
   8.2303856819767232e+17
   9.0959829002668813e+17
   1.0052615772688342e+18
   1.1109858602563712e+18
   1.2278292631485970e+18
   1.3569611939940810e+18
   1.4996740485594657e+18
   1.6573961450606881e+18
   1.8317060192517601e+18
   2.0243482229411579e+18
   2.2372507840526853e+18
   2.4725445029769687e+18
   2.7325842783379528e+18
   3.0199726756098365e+18
   3.3375859744670935e+18
   3.6886029555582029e+18
   4.0765367148108068e+18
   4.5052698236765440e+18
   4.9790931872111176e+18
   5.5027489888943135e+18
   6.0814781519961702e+18
   6.7210727924986010e+18
   7.4279341885389363e+18
   8.2091368465530675e+18
   9.0724993053136814e+18
   1.0026662386494198e+19
   1.1081175674916358e+19
   1.2246593094004847e+19
   1.3534578533000223e+19
   1.4958022583082809e+19
   1.6531171550741899e+19
   1.8269770039599454e+19
   2.0191218527695090e+19
   2.2314747517318812e+19
   2.4661610000341512e+19
   2.7255294165301002e+19
   3.0121758475087553e+19
   3.3289691467965432e+19
   3.6790798882106413e+19
   4.0660120977274061e+19
   4.4936383229520880e+19
   4.9662383908768727e+19
   5.4885422418279211e+19
   6.0657772682979369e+19
   6.7037206324472250e+19
   7.4087570858843619e+19
   8.1879428704062800e+19
   9.0490763392378634e+19
   1.0000776005572131e+20
   1.1052566799547061e+20
   1.2214975396947848e+20
   1.3499635573716302e+20
   1.4919404640690720e+20
   1.6488492123894242e+20
   1.8222601978247286e+20
   2.0139089758026668e+20
   2.2257136317086207e+20
   2.4597939777289005e+20
   2.7184927686435983e+20
   3.0043991489038549e+20
   3.3203745656597683e+20
   3.6695814070852361e+20
   4.0555146526217175e+20
   4.4820368519071852e+20
   4.9534167824711496e+20
   5.4743721730949612e+20
   6.0501169204271369e+20
   6.6864132714134700e+20
   7.3896294938012195e+20
   8.1668036119031775e+20
   9.0257138455105503e+20
   9.9749564569309793e+20
   1.1024031785271021e+21
   1.2183439329023096e+21
   1.3464782828575407e+21
   1.4880886400345899e+21
   1.6445922884849697e+21
   1.8175555693250643e+21
   2.0087095572044880e+21
   2.2199673854830119e+21
   2.4534433935122555e+21
   2.7114742876545719e+21
   2.9966425278257163e+21
   3.3118021736216768e+21
   3.6601074487063940e+21
   4.0450443093423623e+21
   4.4704653330125727e+21
   4.9406282763108614e+21
];
%
%  Next, setup the weights.
%
WT0=[
  0.21035620538389819885E-28
 -0.12644693616088940552E-13
  0.46157312567885668321E-13
 -0.27987033742576678494E-13
  0.54657649654108409156E-13
 -0.26529331099287291499E-13
  0.56749134340673213135E-13
 -0.21572768289772080733E-13
  0.58318460867739760925E-13
 -0.15465892848687829700E-13
  0.60573024556529743179E-13
 -0.85025312590830646706E-14
  0.63880180611476449908E-13
 -0.56596576350102877128E-15
  0.68485006047914070374E-13
  0.85728977321682762439E-14
  0.74650681546818133979E-13
  0.19208372932613381433E-13
  0.82693454289757706437E-13
  0.31701165629228998860E-13
  0.93000040396952081623E-13
  0.46490696394179916916E-13
  0.10604419444905640479E-12
  0.64112165895974571186E-13
  0.12240608340017008854E-12
  0.85217767515070225126E-13
  0.14279579404871719178E-12
  0.11060266069684630524E-12
  0.16808202030984049793E-12
  0.14123670281595459178E-12
  0.19932710117763080694E-12
  0.17830320429641190908E-12
  0.23782981967994220694E-12
  0.22324626027650970672E-12
  0.28517768171070337866E-12
  0.27782855757859949941E-12
  0.34331077352335570356E-12
  0.34420197641399489491E-12
  0.41459976123278392377E-12
  0.42499381982249688902E-12
  0.50194116319499513568E-12
  0.52341213107734220974E-12
  0.60887371932475596308E-12
  0.64337432538102752823E-12
  0.73972052807969328662E-12
  0.78966429788659195358E-12
  0.89976265596232089469E-12
  0.96812431295714248033E-12
  0.10954511874715580138E-11
  0.11858893754914550946E-11
  0.13346662261643230004E-11
  0.14516734901185850330E-11
  0.16270332417806469835E-11
  0.17761192965262209981E-11
  0.19843094598657728049E-11
  0.21722251127126358340E-11
  0.24208558013562638359E-11
  0.26558665246211218806E-11
  0.29542133130007581834E-11
  0.32464334551104291636E-11
  0.36058072230543221354E-11
  0.39676082798212998573E-11
  0.44018068787207206107E-11
  0.48483162182209038194E-11
  0.53741760778847699533E-11
  0.59238861421284167256E-11
  0.65619559488550109069E-11
  0.72374683888302430790E-11
  0.80128318647926483616E-11
  0.88417664804021047188E-11
  0.97850472788000974558E-11
  0.10801152249024990298E-10
  0.11949741288775620211E-10
  0.13194249255521529763E-10
  0.14593803746893311587E-10
  0.16117088182600741411E-10
  0.17823362499440759998E-10
  0.19686960839662181492E-10
  0.21768042712348459703E-10
  0.24047127453754399616E-10
  0.26586169224245549843E-10
  0.29372566166660639848E-10
  0.32471120715874496926E-10
  0.35876995485484108464E-10
  0.39659090711123993131E-10
  0.43821451522206353096E-10
  0.48438566886024388833E-10
  0.53524764256840081139E-10
  0.59161909123775483693E-10
  0.65376353273289546079E-10
  0.72259490983916676274E-10
  0.79851856505621554776E-10
  0.88256972132554113981E-10
  0.97532219231111761839E-10
  0.10779639493701439540E-09
  0.11912700941828918848E-09
  0.13166195190543400137E-09
  0.14550289515667290553E-09
  0.16081145810919829686E-09
  0.17771842706736448843E-09
  0.19641479168712890847E-09
  0.21706652163468219836E-09
  0.23990084518390418007E-09
  0.26512635046403231588E-09
  0.29301487204485379059E-09
  0.32382671796405940321E-09
  0.35788852978338560143E-09
  0.39552347102192978760E-09
  0.43712543089936042400E-09
  0.48309404739375797825E-09
  0.53390563500721587301E-09
  0.59005295736900699180E-09
  0.65211327580989620296E-09
  0.72069283339348376395E-09
  0.79649244503723412979E-09
  0.88025670846750671949E-09
  0.97283758951862838429E-09
  0.10751484374562240226E-08
  0.11882260626931209546E-08
  0.13131897062580609792E-08
  0.14513021636655630472E-08
  0.16039339435115999074E-08
  0.17726240632935599375E-08
  0.19590497332198700268E-08
  0.21650875406672578739E-08
  0.23927891159868709863E-08
  0.26444435360147638239E-08
  0.29225595734392348364E-08
  0.32299302912485969586E-08
  0.35696226515762041770E-08
  0.39450454481729238111E-08
  0.43599472612559631445E-08
  0.48184890913636877388E-08
  0.53252519017629532640E-08
  0.58853155833436851597E-08
  0.65042776355473761329E-08
  0.71883404192425699700E-08
  0.79443429030829021263E-08
  0.87798585629591097164E-08
  0.97032425780220928695E-08
  0.10723743227689399320E-07
  0.11851567478400949647E-07
  0.13098009332253010086E-07
  0.14475537424021579590E-07
  0.15997944513720238914E-07
  0.17680461540553181382E-07
  0.19539933354871179973E-07
  0.21594964684505000780E-07
  0.23866128306162058448E-07
  0.26376149610345229076E-07
  0.29150154762698309478E-07
  0.32215902055657983021E-07
  0.35604079260985092872E-07
  0.39348591789544360374E-07
  0.43486920453657868292E-07
  0.48060478694380069026E-07
  0.53115044437493142888E-07
  0.58701201380017479006E-07
  0.64874861635712681668E-07
  0.71697809408859342572E-07
  0.79238334805050121668E-07
  0.87571902294266750180E-07
  0.96781920558355529704E-07
  0.10696056312048640378E-06
  0.11820970459254829996E-06
  0.13064192692376749733E-06
  0.14438165911984879199E-06
  0.15956640998357920017E-06
  0.17634815657222329376E-06
  0.19489485336503760544E-06
  0.21539212473518470874E-06
  0.23804511154682969141E-06
  0.26308053514448537707E-06
  0.29074895589985508713E-06
  0.32132729115584617374E-06
  0.35512157675298091583E-06
  0.39247003960676882480E-06
  0.43374647340783931675E-06
  0.47936398879211592319E-06
  0.52977913903701991577E-06
  0.58549649799821895163E-06
  0.64707370170465826468E-06
  0.71512703747583724585E-06
  0.79033760405820345763E-06
  0.87345813593705873095E-06
  0.96532052953043588492E-06
  0.10668441762993300587E-05
  0.11790451573234470162E-05
  0.13030464193827640161E-05
  0.14400890071820729325E-05
  0.15915444905568059534E-05
  0.17589286853766969218E-05
  0.19439168303884759854E-05
  0.21483603476946079820E-05
  0.23743053781113989137E-05
  0.26240132540945578888E-05
  0.28999831372928098690E-05
  0.32049770257732388059E-05
  0.35420474020978068237E-05
  0.39145677786671023699E-05
  0.43262664657480666650E-05
  0.47812638810077579977E-05
  0.52841137925459449193E-05
  0.58398488901504216978E-05
  0.64540311583956847791E-05
  0.71328075387126483640E-05
  0.78829714540444842174E-05
  0.87120307957926337232E-05
  0.96282830690788715692E-05
  0.10640898433258409191E-04
  0.11760011483484689510E-04
  0.12996822677619209683E-04
  0.14363710437469910021E-04
  0.15874355032820411062E-04
  0.17543875501207770228E-04
  0.19388980961051210830E-04
  0.21428137844875080257E-04
  0.23681754714302638776E-04
  0.26172386518181301401E-04
  0.28924960326686671307E-04
  0.31967024811679323171E-04
  0.35329025959273332330E-04
  0.39044611784553647914E-04
  0.43150969087563300697E-04
  0.47689195631023510063E-04
  0.52704711454198507913E-04
  0.58247713449355698436E-04
  0.64373677739542765810E-04
  0.71143914895687183756E-04
  0.78626183537770216584E-04
  0.86895368472094227166E-04
  0.96034230136811580777E-04
  0.10613423286950590218E-03
  0.11729646026571620215E-03
  0.12963262680752149015E-03
  0.14326619585462989793E-03
  0.15833361521503999387E-03
  0.17498568260644121238E-03
  0.19338905472215898854E-03
  0.21372791490647819607E-03
  0.23620581621746139576E-03
  0.26104771814107751186E-03
  0.28850223750643329973E-03
  0.31884413578705622606E-03
  0.35237706800135917308E-03
  0.38943662007281590462E-03
  0.43039366566958368833E-03
  0.47565807487700027594E-03
  0.52568281303167982161E-03
  0.58096846836193664822E-03
  0.64206825611437718954E-03
  0.70959354469423893946E-03
  0.78421996374912559657E-03
  0.86669414659246905939E-03
  0.95784118347590099622E-03
  0.10585728435434579864E-02
  0.11698966653840010659E-02
  0.12929259749833240008E-02
  0.14288909657393959393E-02
  0.15791508896276298946E-02
  0.17452075486231769744E-02
  0.19287201026563420645E-02
  0.21315214730933669356E-02
  0.23556362776736309780E-02
  0.26033007312297349808E-02
  0.28769842725422829777E-02
  0.31794136300501061286E-02
  0.35135987233342810820E-02
  0.38828616256736492134E-02
  0.42908672540316900729E-02
  0.47416579734868603835E-02
  0.52396893384671194144E-02
  0.57898709851889693795E-02
  0.63976070720003401851E-02
  0.70688437831027441105E-02
  0.78101127966257229834E-02
  0.86285849757206656285E-02
  0.95321125348557418644E-02
  0.10529286971052299882E-01
  0.11629470426845669312E-01
  0.12842853010275019979E-01
  0.14180454004084899408E-01
  0.15654168426849128515E-01
  0.17276700258609969246E-01
  0.19061578784670331344E-01
  0.21022951738450229575E-01
  0.23175536220878610594E-01
  0.25534136813168170632E-01
  0.28113470567056909888E-01
  0.30927161344728969217E-01
  0.33987341082207328524E-01
  0.37302669043183162012E-01
  0.40877565846017260842E-01
  0.44708454582277110112E-01
  0.48782456555642082774E-01
  0.53070463888373879680E-01
  0.57525214796560822372E-01
  0.62068889249144526543E-01
  0.66590986905823254527E-01
  0.70926870977039285782E-01
  0.74857621553299974471E-01
  0.78074664211727651253E-01
  0.80188872113380424422E-01
  0.80676406709186576638E-01
  0.78917673064227769619E-01
  0.74124063016304961304E-01
  0.65458647531413310938E-01
  0.51957717257333460581E-01
  0.32847972741848592559E-01
  0.74970763258312700383E-02
 -0.23866128698945488634E-01
 -0.60174943784761181220E-01
 -0.98178997988506697125E-01
 -0.13281477972726110637E+00
 -0.15546285697725620301E+00
 -0.15639821579874499391E+00
 -0.12430498665290500016E+00
 -0.54868159863436967438E-01
  0.46862558991703072431E-01
  0.15112182958062059246E+00
  0.21193155344105990556E+00
  0.16951341358877961008E+00
  0.13861974203314799889E-01
 -0.18693504518381348634E+00
 -0.24558896069253360883E+00
 -0.53092694018998139172E-01
  0.25199984157985949595E+00
  0.19682248760574280744E+00
 -0.20143336189691599114E+00
 -0.24584508272586030886E+00
  0.34335593140766357267E+00
 -0.47700665106262918336E-01
 -0.20966284735074519618E+00
  0.25090851237942479734E+00
 -0.16764126613351920669E+00
  0.74951997868967626393E-01
 -0.16951351027511458308E-01
 -0.88285242013749955919E-02
  0.16167873561835501006E-01
 -0.15564117559925290737E-01
  0.12513111252163700016E-01
 -0.92993152469780984010E-02
  0.66505259967794931597E-02
 -0.46677074086718837662E-02
  0.32488855455568548675E-02
 -0.22555291085446750252E-02
  0.15668752111009150857E-02
 -0.10911324863680910667E-02
  0.76253830490409491537E-03
 -0.53525405778204085232E-03
  0.37771045672146697511E-03
 -0.26825282115352968703E-03
  0.19202652749967551236E-03
 -0.13882137232089440237E-03
  0.10159989178100049974E-03
 -0.75497250626312589804E-04
  0.57141729200529359656E-04
 -0.44191609676806338589E-04
  0.35017761995785832974E-04
 -0.28484969541399899885E-04
  0.23801048173241150913E-04
 -0.20412749240321660898E-04
  0.17933634897918600856E-04
 -0.16093677807779768610E-04
  0.14703946920177649808E-04
 -0.13632085400878389976E-04
  0.12785408514801159958E-04
 -0.12099103232765920005E-04
  0.11527807463720849231E-04
 -0.11039642499677090445E-04
  0.10612147244787119201E-04
 -0.10229546722558040478E-04
  0.98808019172690629615E-05
 -0.95581362683731644660E-05
  0.92559894628748428139E-05
 -0.89703680398349605195E-05
  0.86984392262176360162E-05
 -0.84381998993781744882E-05
  0.81881850397080960968E-05
 -0.79472759767526677708E-05
  0.77146204892770847038E-05
 -0.74895903346701209313E-05
  0.72717128455631313132E-05
 -0.70605952027200397989E-05
  0.68558909586342899260E-05
 -0.66573078263653524112E-05
  0.64646106359474119065E-05
 -0.62775967052488337714E-05
  0.60960690045641779074E-05
 -0.59198350149128057049E-05
  0.57487213891732378171E-05
 -0.55825764283683409239E-05
  0.54212557290731716876E-05
 -0.52646114832336199488E-05
  0.51124978277574427215E-05
 -0.49647806501954759280E-05
  0.48213365261592663842E-05
 -0.46820435304229967753E-05
  0.45467775878200409083E-05
 -0.44154180000760036035E-05
  0.42878525031129105819E-05
 -0.41639746681098469138E-05
  0.40436784078298272207E-05
 -0.39268575265170754164E-05
  0.38134098513603608910E-05
 -0.37032392398348510301E-05
  0.35962529271829669603E-05
 -0.34923586009061261710E-05
  0.33914651799141980567E-05
 -0.32934854312419398030E-05
  0.31983363564463991755E-05
 -0.31059371103484061739E-05
  0.30162076790696608685E-05
 -0.29290699115622908088E-05
  0.28444489540948809787E-05
 -0.27622729358204111026E-05
  0.26824715684158770123E-05
 -0.26049757091405561401E-05
  0.25297182336640620443E-05
 -0.24566346909850541750E-05
  0.23856627938558707926E-05
 -0.23167415846242317899E-05
  0.22498114158806991699E-05
 -0.21848145526023768863E-05
  0.21216953789870668849E-05
 -0.20603999269499001719E-05
  0.20008754358939798722E-05
 -0.19430704833474811314E-05
  0.18869353395796189952E-05
 -0.18324219538318030139E-05
  0.17794835990689859392E-05
 -0.17280746718615620423E-05
  0.16781508408187099176E-05
 -0.16296692208678659180E-05
  0.15825882749397249515E-05
 -0.15368675777174579453E-05
  0.14924677441823879966E-05
 -0.14493505440316709914E-05
  0.14074789632159439540E-05
 -0.13668170916979060623E-05
  0.13273299804037580486E-05
 -0.12889836283545319635E-05
  0.12517450527435769451E-05
 -0.12155822896230419416E-05
  0.11804642977313919266E-05
 -0.11463608782321180614E-05
  0.11132426817551759487E-05
 -0.10810812413606430626E-05
  0.10498489463079790445E-05
 -0.10195189699197359345E-05
  0.99006522696480728454E-06
 -0.96146238238294861780E-06
  0.93368585903962674737E-06
 -0.90671180391877650539E-06
  0.88051703794158972103E-06
 -0.85507903355065637431E-06
  0.83037591856657379114E-06
 -0.80638646950736030707E-06
  0.78309007959684342440E-06
 -0.76046672516600126341E-06
  0.73849695280747625374E-06
 -0.71716187755802708277E-06
  0.69644316956135183135E-06
 -0.67632302736954520768E-06
  0.65678415552512063648E-06
 -0.63780975556830963862E-06
  0.61938351997573715339E-06
 -0.60148961697207375561E-06
  0.58411266955661727328E-06
 -0.56723774005096348097E-06
  0.55085032232163946943E-06
 -0.53493633333997110944E-06
  0.51948209868244907623E-06
 -0.50447433644131324069E-06
  0.48990014588904722006E-06
 -0.47574699997612418437E-06
  0.46200273605516768571E-06
 -0.44865554313913170035E-06
  0.43569394953682061395E-06
 -0.42310681390672771474E-06
  0.41088331789010578938E-06
 -0.39901295702377770089E-06
  0.38748552999527151342E-06
 -0.37629112895780040410E-06
  0.36542013202237861199E-06
 -0.35486319617658597678E-06
  0.34461124894551899091E-06
 -0.33465547948734137742E-06
  0.32498733079349469753E-06
 -0.31559849314212400870E-06
  0.30648089749286697872E-06
 -0.29762670812339488137E-06
  0.28902831526969758359E-06
 -0.28067832866533619338E-06
  0.27256957173628210936E-06
 -0.26469507560502209121E-06
  0.25704807272840818962E-06
 -0.24962199077868072385E-06
  0.24241044716871418470E-06
 -0.23540724388966429376E-06
  0.22860636218196940360E-06
 -0.22200195709699679458E-06
  0.21558835236008589946E-06
 -0.20936003566018609291E-06
  0.20331165407816460064E-06
 -0.19743800942055008854E-06
  0.19173405358761880565E-06
 -0.18619488421725469657E-06
  0.18081574059851621267E-06
 -0.17559199964945570756E-06
  0.17051917187001359107E-06
 -0.16559289739464161294E-06
  0.16080894226704971095E-06
 -0.15616319488355259655E-06
  0.15165166247713490835E-06
 -0.14727046762624959541E-06
  0.14301584488191908927E-06
 -0.13888413756276189610E-06
  0.13487179465891340034E-06
 -0.13097536777476139733E-06
  0.12719150812468619494E-06
 -0.12351696364172060964E-06
  0.11994857620999400758E-06
 -0.11648327897320799497E-06
  0.11311809368621049480E-06
 -0.10985012813125349456E-06
  0.10667657363213940060E-06
 -0.10359470265995370341E-06
  0.10060186649751670646E-06
 -0.97695492950586510549E-07
  0.94873084124749734212E-07
 -0.92132214283314405060E-07
  0.89470527774724840543E-07
 -0.86885737009431156146E-07
  0.84375620484431720137E-07
 -0.81938020868743667602E-07
  0.79570843154679285033E-07
 -0.77272052863928245642E-07
  0.75039674297478696595E-07
 -0.72871788831488877790E-07
  0.70766533266929795021E-07
 -0.68722098232649864085E-07
  0.66736726633350816695E-07
 -0.64808712137167796518E-07
  0.62936397705674278282E-07
 -0.61118174170089166613E-07
  0.59352478851276259248E-07
 -0.57637794217708362227E-07
  0.55972646579199561646E-07
 -0.54355604818630417116E-07
  0.52785279162883091201E-07
 -0.51260319990189909104E-07
  0.49779416670227040800E-07
 -0.48341296436233499673E-07
  0.46944723290471018939E-07
 -0.45588496942848318593E-07
  0.44271451780258727462E-07
 -0.42992455864451837250E-07
  0.41750409958276353087E-07
 -0.40544246580819631899E-07
  0.39372929090650920317E-07
 -0.38235450795273419860E-07
  0.37130834085523768675E-07
 -0.36058129594866006882E-07
  0.35016415383563508384E-07
 -0.34004796146768310693E-07
  0.33022402445156628786E-07
 -0.32068389957362020693E-07
  0.31141938754091167181E-07
 -0.30242252593599009938E-07
  0.29368558237615870861E-07
 -0.28520104786762129308E-07
  0.27696163034962460481E-07
 -0.26896024842649531045E-07
  0.26119002528299780677E-07
 -0.25364428277531399452E-07
  0.24631653569082881086E-07
 -0.23920048617305701087E-07
  0.23229001830889051246E-07
 -0.22557919287333441164E-07
  0.21906224222551861542E-07
 -0.21273356535101449785E-07
  0.20658772304732728872E-07
 -0.20061943324939738367E-07
  0.19482356649058811450E-07
 -0.18919514149424371449E-07
  0.18372932089200169455E-07
 -0.17842140706600119178E-07
  0.17326683811178871537E-07
 -0.16826118391794460252E-07
  0.16340014235853439520E-07
 -0.15867953559529469194E-07
  0.15409530648689289756E-07
 -0.14964351510224189348E-07
  0.14532033533448680711E-07
 -0.14112205161252269515E-07
  0.13704505570744190592E-07
 -0.13308584363144380637E-07
  0.12924101262648650824E-07
 -0.12550725823983400048E-07
  0.12188137148392460070E-07
 -0.11836023607829919776E-07
  0.11494082577135149689E-07
 -0.11162020173950579471E-07
  0.10839551006143939185E-07
 -0.10526397926519000691E-07
  0.10222291794616780561E-07
 -0.99269712454067339727E-08
  0.96401824646609211080E-08
 -0.93616789708108062008E-08
  0.90912214031028235422E-08
 -0.88285773158782117222E-08
  0.85735209788006156203E-08
 -0.83258331828536520172E-08
  0.80853010519388850005E-08
 -0.78517178599159764754E-08
  0.76248828529317467947E-08
 -0.74046010768841071059E-08
  0.71906832098687622728E-08
 -0.69829453994641022983E-08
  0.67812091047173317459E-08
 -0.65853009426977313290E-08
  0.63950525394835678792E-08
 -0.62103003854524726628E-08
  0.60308856947513333758E-08
 -0.58566542688268000505E-08
  0.56874563638996400563E-08
 -0.55231465622676621705E-08
  0.53635836473256503502E-08
 -0.52086304821955339059E-08
  0.50581538918636123925E-08
 -0.49120245487234261852E-08
  0.47701168614250102390E-08
 -0.46323088669346301918E-08
  0.44984821257128528939E-08
 -0.43685216199214189850E-08
  0.42423156545711511770E-08
 -0.41197557615253826557E-08
  0.40007366062763786183E-08
 -0.38851558974150727748E-08
  0.37729142987165276171E-08
 -0.36639153437652961865E-08
  0.35580653530469988486E-08
 -0.34552733534349709264E-08
  0.33554510000030559535E-08
 -0.32585125000973940249E-08
  0.31643745396017359674E-08
 -0.30729562113327868796E-08
  0.29841789455041530200E-08
 -0.28979664421992758989E-08
  0.28142446057953132032E-08
 -0.27329414812814390995E-08
  0.26539871924168161407E-08
 -0.25773138816751668664E-08
  0.25028556519244291169E-08
 -0.24305485097912959499E-08
  0.23603303106619080563E-08
 -0.22921407052714058628E-08
  0.22259210878365399183E-08
 -0.21616145456867578372E-08
  0.20991658103504500059E-08
 -0.20385212100542671023E-08
  0.19796286235947108123E-08
 -0.19224374355423869606E-08
  0.18668984927404269626E-08
 -0.18129640620596530537E-08
  0.17605877893741759270E-08
 -0.17097246597221680675E-08
  0.16603309586176070024E-08
 -0.16123642344797300183E-08
  0.15657832621478869464E-08
 -0.15205480074504270167E-08
  0.14766195927971890460E-08
 -0.14339602637660350412E-08
  0.13925333566546899287E-08
 -0.13523032669700050132E-08
  0.13132354188275479691E-08
 -0.12752962352352419048E-08
  0.12384531092355140920E-08
 -0.12026743758811379070E-08
  0.11679292850206920657E-08
 -0.11341879748702300079E-08
  0.11014214463484779337E-08
 -0.10696015381534739059E-08
  0.10387009025592480307E-08
 -0.10086929819117270958E-08
  0.97955198580366660312E-09
 -0.95125286890900208679E-09
  0.92377130945756002788E-09
 -0.89708368833163302057E-09
  0.87116706876644941479E-09
 -0.84599917663709281674E-09
  0.82155838131493233922E-09
 -0.79782367707710818756E-09
  0.77477466505309300292E-09
 -0.75239153569281363935E-09
  0.73065505174126690031E-09
 -0.70954653170499806391E-09
  0.68904783379622916597E-09
 -0.66914134034083904798E-09
  0.64980994263679380996E-09
 -0.63103702625001653087E-09
  0.61280645673505631634E-09
 -0.59510256576828648337E-09
  0.57791013768171107199E-09
 -0.56121439638580773724E-09
  0.54500099267016873900E-09
 -0.52925599187102496775E-09
  0.51396586189505248004E-09
 -0.49911746158916965897E-09
  0.48469802944632955612E-09
 -0.47069517263760014005E-09
  0.45709685636110621144E-09
 -0.44389139349867930081E-09
  0.43106743457132480471E-09
 -0.41861395798487421376E-09
  0.40652026055743952578E-09
 -0.39477594832052842209E-09
  0.38337092758591272230E-09
 -0.37229539627057491552E-09
  0.36153983547227521766E-09
 -0.35109500128849910154E-09
  0.34095191687175519511E-09
 -0.33110186471439289704E-09
  0.32153637915631239045E-09
 -0.31224723910912398970E-09
  0.30322646099050710905E-09
 -0.29446629186269331473E-09
  0.28595920276917881651E-09
 -0.27769788226393698380E-09
  0.26967523012757238777E-09
 -0.26188435126501301881E-09
  0.25431854977949852574E-09
 -0.24697132321776888896E-09
  0.23983635698150977934E-09
 -0.23290751890024982077E-09
  0.22617885396104611209E-09
 -0.21964457919042880688E-09
  0.21329907868420528750E-09
 -0.20713689878085250208E-09
  0.20115274337434959303E-09
 -0.19534146936242159446E-09
  0.18969808222628370222E-09
 -0.18421773173808480806E-09
  0.17889570779236339873E-09
 -0.17372743635793169822E-09
  0.16870847554670939670E-09
 -0.16383451179612808855E-09
  0.15910135616182540748E-09
 -0.15450494071744330869E-09
  0.15004131505843449035E-09
 -0.14570664290687508993E-09
  0.14149719881436189859E-09
 -0.13740936496016368823E-09
  0.13343962804187350973E-09
 -0.12958457625588848896E-09
  0.12584089636512460714E-09
 -0.12220537085144369393E-09
  0.11867487515034671199E-09
 -0.11524637496555659739E-09
  0.11191692366118209830E-09
 -0.10868365972922199935E-09
  0.10554380433023220449E-09
 -0.10249465890504320142E-09
  0.99533602855475026004E-10
 -0.96658091292055187607E-10
  0.93865652846806834872E-10
 -0.91153887549225778384E-10
  0.88520464763624569567E-10
 -0.85963121186075071965E-10
  0.83479658899239781350E-10
 -0.81067943483446593655E-10
  0.78725902182445575131E-10
 -0.76451522122414838220E-10
  0.74242848583018786633E-10
 -0.72097983319807716229E-10
  0.70015082938317186447E-10
 -0.67992357322650702613E-10
  0.66028068126905105916E-10
 -0.64120527350693063962E-10
  0.62268096049913068230E-10
 -0.60469183303296924744E-10
  0.58722245716347976830E-10
 -0.57025788118346391824E-10
  0.55378366976800887211E-10
 -0.53778600071222512844E-10
  0.52225190653632309517E-10
 -0.50716985205207501023E-10
  0.49253109171572777739E-10
 -0.47833283755322262690E-10
  0.46458563164427748566E-10
 -0.45133048314954522821E-10
  0.43867868281155817508E-10
 -0.42690428488605550129E-10
  0.41665890740660430845E-10
 -0.40947061319725302450E-10
  0.40890256062850869212E-10
 -0.42324395200868817252E-10
  0.47175970286618397865E-10
 -0.59920514722972034865E-10
  0.90953607290146280299E-10
];
%
% Get the points at which to evaluate F.
%
Y=YBASE/B;
%
% Evaluate F at these points.
%
YF=feval(funname,Y);
%
% Compute the dot product.
%
z0=WT0.'*YF;
%
%  Finally correct for the factor of 1/B.
%
z0=z0./B;

                                                                                          kcovario.m                                                                                          0000664 0002715 0000074 00000000736 10372707673 011551  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = kcovario(x)
%% calculates nadwat(x,covario,rsort,50)

  global covario rsort rmax;
  nx = size(x,1);ny=size(x,2);
  T1 = 2000;
  T2 = rmax;
  res=x(:);
  indsmall = x<=T1;
  indlarge = x>T1 & x<=T2;
  indout = x>T2;
  xsmall = x(indsmall);
  xlarge = x(indlarge);
  res(indsmall) = nadwatloop(xsmall,covario,rsort,50);
  resT1 = nadwatloop(T1,covario,rsort,50);
  res(indlarge) = resT1 * (rmax - xlarge) / (rmax-T1);
  res(indout) = 0;
  res = reshape(res,nx,ny);                                  kcovariox.m                                                                                         0000664 0002715 0000074 00000000756 10300707471 011726  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = kcovariox(x)
%% calculates nadwat(x,covario,rsort,50).*x
  global covario rsort rmax;
  nx = size(x,1);ny=size(x,2);
  T1 = 1500;
  T2 = rmax;
  res=x(:);
  indsmall = x<=T1;
  indlarge = x>T1 & x<=T2;
  indout = x>T2;
  xsmall = x(indsmall);
  xlarge = x(indlarge);
  res(indsmall) = nadwat(xsmall,covario,rsort,50);
  resT1 = nadwat(T1,covario,rsort,50);
  res(indlarge) = resT1 * (rmax - xlarge) / (rmax-T1);
  res(indout) = 0;
  res = res .* x(:);
  res = reshape(res,nx,ny);                  keepathdef.m                                                                                        0000644 0002715 0000074 00000005715 10106474650 012023  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function p = pathdef
%PATHDEF Search path defaults.
%   PATHDEF returns a string that can be used as input to MATLABPATH
%   in order to set the path.

  
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2000/06/01 16:19:21 $

% PATH DEFINED HERE -- Don't change this line.

p = [...
'./mfunctions:',...
'.:',...
matlabroot,'/toolbox/matlab/general:',...
matlabroot,'/toolbox/matlab/ops:',...
matlabroot,'/toolbox/matlab/lang:',...
matlabroot,'/toolbox/matlab/elmat:',...
matlabroot,'/toolbox/matlab/elfun:',...
matlabroot,'/toolbox/matlab/specfun:',...
matlabroot,'/toolbox/matlab/matfun:',...
matlabroot,'/toolbox/matlab/datafun:',...
matlabroot,'/toolbox/matlab/audio:',...
matlabroot,'/toolbox/matlab/polyfun:',...
matlabroot,'/toolbox/matlab/funfun:',...
matlabroot,'/toolbox/matlab/sparfun:',...
matlabroot,'/toolbox/matlab/graph2d:',...
matlabroot,'/toolbox/matlab/graph3d:',...
matlabroot,'/toolbox/matlab/specgraph:',...
matlabroot,'/toolbox/matlab/graphics:',...
matlabroot,'/toolbox/matlab/uitools:',...
matlabroot,'/toolbox/matlab/strfun:',...
matlabroot,'/toolbox/matlab/iofun:',...
matlabroot,'/toolbox/matlab/timefun:',...
matlabroot,'/toolbox/matlab/datatypes:',...
matlabroot,'/toolbox/matlab/verctrl:',...
matlabroot,'/toolbox/matlab/demos:',...
matlabroot,'/toolbox/local:',...
matlabroot,'/toolbox/simulink/simulink:',...
matlabroot,'/toolbox/simulink/blocks:',...
matlabroot,'/toolbox/simulink/simdemos:',...
matlabroot,'/toolbox/simulink/simdemos/aerospace:',...
matlabroot,'/toolbox/simulink/simdemos/automotive:',...
matlabroot,'/toolbox/simulink/simdemos/simfeatures:',...
matlabroot,'/toolbox/simulink/simdemos/simgeneral:',...
matlabroot,'/toolbox/simulink/simdemos/simnew:',...
matlabroot,'/toolbox/simulink/dee:',...
matlabroot,'/toolbox/stateflow/stateflow:',...
matlabroot,'/toolbox/stateflow/sfdemos:',...
matlabroot,'/toolbox/stateflow/coder:',...
matlabroot,'/toolbox/compiler:',...
matlabroot,'/toolbox/finance/finance:',...
matlabroot,'/toolbox/finance/calendar:',...
matlabroot,'/toolbox/finance/findemos:',...
matlabroot,'/toolbox/finance/finsupport:',...
matlabroot,'/toolbox/ident/ident:',...
matlabroot,'/toolbox/ident/idobsolete:',...
matlabroot,'/toolbox/ident/idguis:',...
matlabroot,'/toolbox/ident/idutils:',...
matlabroot,'/toolbox/ident/iddemos:',...
matlabroot,'/toolbox/ident/idhelp:',...
matlabroot,'/toolbox/nnet/nnet:',...
matlabroot,'/toolbox/nnet/nnutils:',...
matlabroot,'/toolbox/nnet/nncontrol:',...
matlabroot,'/toolbox/nnet/nndemos:',...
matlabroot,'/toolbox/nnet/nnobsolete:',...
matlabroot,'/toolbox/optim:',...
matlabroot,'/toolbox/pde:',...
matlabroot,'/toolbox/sb2sl:',...
matlabroot,'/toolbox/signal/signal:',...
matlabroot,'/toolbox/signal/fdatoolgui:',...
matlabroot,'/toolbox/signal/sptoolgui:',...
matlabroot,'/toolbox/signal/sigdemos:',...
matlabroot,'/toolbox/stats:',...
matlabroot,'/toolbox/symbolic:',...
matlabroot,'/toolbox/wavelet/wavelet:',...
matlabroot,'/toolbox/wavelet/wavedemo:',...
     ...
];

p = [userpath,p];
                                                   kkcova.m                                                                                            0000664 0002715 0000074 00000003217 10350112327 011165  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   cova=[1 0.646687900976095 0.686815029195937 0.445716035532842 0.786425276263434 0.437232717857406 0.439780650605974 0.701311755457254 0.57635920393054 0.44296049137226
0.646687900976095 1 0.75047942607273 0.288435275272244 0.532729420752751 0.592551300967742 0.679944086962346 0.514017442619253 0.376113714512208 0.626191864394216
0.686815029195937 0.75047942607273 0.99999865585614 0.322393902222855 0.629857232229802 0.444791641593239 0.558205410470006 0.490198364032628 0.439737541688365 0.632309195554803
0.445716035532842 0.288435275272244 0.322393902222855 0.99999865585614 0.509678097657353 0.212625724751307 0.196123219823394 0.490316153564421 0.709675563736769 0.203855262783624
0.786425276263434 0.532729420752751 0.629857232229802 0.509678097657353 1 0.345050792049102 0.36642277089828 0.610681711438723 0.697346464882231 0.398778509439082
0.437232717857406 0.592551300967742 0.444791641593239 0.212625724751307 0.345050792049102 1 0.6118301657595 0.429412704536913 0.256136983635375 0.43804944687928
0.439780650605974 0.679944086962346 0.558205410470006 0.196123219823394 0.36642277089828 0.6118301657595 1 0.358954446562308 0.257035529239854 0.695009477970135
0.701311755457254 0.514017442619253 0.490198364032628 0.490316153564421 0.610681711438723 0.429412704536913 0.358954446562308 1 0.543841623480445 0.328702604955391
0.57635920393054 0.376113714512208 0.439737541688365 0.709675563736769 0.697346464882231 0.256136983635375 0.257035529239854 0.543841623480445 0.99999865585614 0.279005247101719
0.44296049137226 0.626191864394216 0.632309195554803 0.203855262783624 0.398778509439082 0.43804944687928 0.695009477970135 0.328702604955391 0.279005247101719 1];
                                                                                                                                                                                                                                                                                                                                                                                 kkk.m                                                                                               0000664 0002715 0000074 00000052433 10350112133 010466  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   Z=[-0.788493895276293 0.703673079805572 1.27457291385536 -0.544677048651387 -2.57843614037921 1.41796893895252 -0.929346478221447 0.084879893073441 -0.774729849597278 -0.78776654585701
-1.32478212503048 0.435001401641063 1.11981310689353 1.09553867346312 -1.27501211779994 1.30610048669975 -0.378260129500752 0.779594578739589 0.43952559976753 -1.53170402987404
-1.00136903763839 -0.159328458980459 1.83297709368828 -0.152729915899304 -1.50442039402892 -0.450236501280776 -0.558629843871903 -0.339876841792943 -0.915293815363893 -0.523326517450876
-0.919503156316905 -0.0094608729463225 2.82901885079767 -0.874580601074944 -2.33415914400688 -1.47957172666194 -1.50542244778464 0.690972975641348 -1.08991040840178 -0.107172086945141
-1.85924789367361 1.56919937552244 1.59705845093069 -0.0422235345211493 -2.57272525427863 0.96972820956773 2.16007629477509 0.581789419846204 -1.11025618585130 -1.85814988397718
-1.65089294947906 0.918841457136271 1.26432373241784 0.089571688851032 -1.93606806066669 0.00176988204495933 1.54115776125846 0.981162661127823 -0.229373172670579 -0.974137969945427
-0.696259017720218 0.91136003880062 0.396128007090506 -0.04662298361487 -2.22857830073966 -0.371226174309621 1.15793352916446 0.359942589862786 -1.08215828075584 -0.37419803801658
-0.604815792901481 -0.0493115055552741 0.699740247395695 0.801501680038555 -1.19917076663761 0.402883224192577 0.194338041265427 0.541207731973134 -0.727003034127508 -0.128861936751544
-0.955072334596684 1.18765952132572 -0.375820969251444 0.0362852649228608 -1.78607475840603 -1.14510535043337 1.69384106715120 1.09968620441966 -0.269725414777426 -1.32558844749239
-0.62039320094778 -0.193744620706567 0.229324541555850 -0.0534358669083121 -1.09887309618330 0.49455565604163 0.156227453892823 -0.830535226829382 -1.73852438315748 -1.99252534445159
-0.281002250941265 0.749657702227045 0.389131337618248 0.267738196040416 -2.43515351933685 0.554208535837819 1.38802541291162 1.51032010926201 0.358526877716164 -0.468878193020087
0.288847021325898 1.38471244364179 1.55215621230429 0.67861782133849 -1.74408645445487 0.083631000174013 2.46011697586854 0.05385849204358 -0.074384664731516 -1.01081391642636
-0.710413125784926 -0.483697494113187 0.257538260501924 0.51822688158313 -1.37217037911084 -1.01377368942824 1.40271160267503 1.17481615562206 -0.466437268898005 -1.48984913131736
-1.03344324403353 0.617847851664557 0.92732277359509 0.713711398810903 -1.46589032140105 -0.150181653578516 1.89177917230482 -1.68220038437630 1.02933716423106 -1.25905962975123
-1.01991874498473 -1.49494449587658 0.932214820918997 0.396721402556328 -1.19503961572755 -0.239542405702272 0.608362774926162 -2.11719254210186 1.00451291566068 0.327432586877264
0.37627478540549 -0.117089799135641 -0.6056332299589 0.308745714637596 -0.712476683807704 -1.39133239356956 0.867257735986296 -1.41298503564843 -0.454864963850137 -0.514626644566037
0.585459734389373 0.756359065609719 0.656574702023516 0.626860310711401 0.152137929865419 0.0134349712749856 0.72094272389864 -1.82956793106910 -0.965038638118149 -1.41318448364547
0.414429891318755 0.0814823555117243 1.24290830085455 0.649884535823153 0.620041256712626 -0.459021600055382 -0.150387977597983 -1.25485038737382 0.980196830071985 -0.886449493992028
-0.781820107438392 0.679511599920276 2.66634802278774 0.122170427311907 0.954339059348012 -0.360268603440747 -0.395854098452980 -2.50157728682052 1.00164831310170 0.218930463664936
0.463916636897944 -1.33388535027163 1.43578159198469 -0.0839939414414513 -0.237223296705353 -1.02354337953040 -0.611528808500623 -2.36212079103554 0.264535537504946 -0.147290684716003
1.82084601447108 0.707014699998518 2.31960228704755 0.716805429027164 -0.848890529227975 0.685603373265515 0.0538149871774566 -1.97788591086463 -0.717234546865008 -0.449869315560658
0.224274476460983 0.308755258464433 0.98411335988261 1.18077543429871 0.145998363766580 1.39827456290953 -0.104283612133143 -1.25887939177878 0.842610075394777 1.36258719328414
0.759367123244609 -1.46350211593832 1.78940866207486 -0.314006300484185 0.00226902161866542 -1.28979598929439 -0.191384585858511 -0.38687898145549 -0.67128363622201 1.52187611279426
-0.202994143620232 -1.85995221682814 -0.249892862870551 -0.695951103892229 0.460012297503655 -0.86654704041213 1.45446701815913 -1.56140218599337 0.380490438751839 1.46448802143683
1.20762184274553 0.211181208612445 -0.458241890469967 1.04472470176989 -0.834460210609009 1.34116348889835 0.0297747806898491 -0.720180960982304 0.0979403958145948 -0.097044731553732
-0.0649044526797247 1.30736701563849 2.12463113470602 0.364057646191854 0.536729572376107 0.278274250269733 -0.0767714336376494 -0.580638517564207 -1.30967364892428 1.01491963518752
0.994386842882293 -0.85693963801297 0.976988331854522 -0.0587334461276961 1.88238558843520 -0.165124510426159 -0.70370926058534 -0.0607472108626095 -0.576767444977683 -0.0271640941481722
0.146143734608190 -1.77230191985274 1.32176240767040 0.44677273837525 0.554543418988427 0.616727267503397 -0.305813184082976 -2.27001976985917 0.270844203367072 0.463446499652172
0.961164468492097 0.608411572346951 -0.655786402524531 -0.805888196413541 0.348096801484232 -0.097169321848472 0.203318902128522 -0.0531216830419947 0.0469040310587995 -0.160469571866869
1.53603869156626 0.74108138685757 0.238496588780542 0.60808629166119 1.49953712976665 0.568389285497634 0.885714113307944 -1.50976951405374 0.497408493247946 0.244116350022334
-0.90860198292105 -0.401468837159225 1.21058201018178 -1.50813994084041 2.49917000583559 -1.55412801148230 -0.704167260510139 0.0420728509957167 -0.201165091176519 -0.865937754496941
0.279724792254461 -1.05253293534557 -0.653070964560592 -0.189057222199969 2.06769985143937 -0.739913269966281 -1.27651331316043 0.438270315724929 -0.678995688038161 -1.04527157734682
-0.995528702375984 -0.216217369871992 1.36995621194852 -0.225467235936849 1.23934477053455 -0.298596166332391 -0.640617838363356 1.01539026108256 -1.00995136957492 -0.630588101639
0.0623905421118769 -0.324635071671329 -0.32505700986735 0.742138333878722 1.59725738000935 -1.20485537307836 -1.51196755641860 1.20865149761129 -1.64526811629296 -0.613148489159018
-0.0865170362299741 -0.3714830838792 -0.204594823263688 -2.85370085770399 1.68821249821923 -0.255077926605553 -1.03898100590696 0.92189253431435 -0.830468886234805 -1.70441537043112
-0.201757112035181 -0.266086019980574 -1.18005531163574 0.0355105551222931 0.88302731778461 -0.367819803137330 -1.05697187437981 0.731714706689743 -2.36836832409119 -0.881659933964649
0.752461298814913 -1.1149122806708 0.896572765080195 -2.02233371107790 2.09211829865169 -0.459826563065201 -0.560939599096234 1.30556438112356 -0.448803563486222 0.142798101587037
0.94058637342057 -1.30338354450569 0.528366098214045 -1.90424673473425 2.50152683807771 -0.602092159987931 0.0175119355085462 0.820609331382864 -1.41547040693117 1.02818297434545
0.253699072042711 0.365247800473644 0.692874891674611 -0.580883996017248 1.16249459085767 1.03391021242367 0.801770397553727 0.79980059779414 -1.56230075944043 -0.920400147962153
-0.0646376985358253 -0.240843750218413 0.582754310735184 -2.11586978591086 1.75085345032489 -2.66217965684922 1.51431944932574 0.178274044904885 0.0213117374184047 0.55769457952044
-0.230377158608337 -0.645556210123218 -0.212415070271665 -1.12026263327328 -1.71379236342553 -1.84852721072548 0.567528086495315 -1.20691764643695 0.524601119122924 0.95979057980284
1.06740288832981 -0.496260262542395 -0.438444391626934 1.17675749725425 -1.95325083077242 -1.86753528124117 0.212192813956959 -0.425025954269585 -0.291143610211395 -0.237688792550707
1.23493125593307 0.377654243344872 1.37849680084274 -0.417540795427515 -1.45680458121607 -2.52878983556371 0.460607145342866 2.46346960330801 -0.553362300078104 -1.02958219413597
0.0406645353091731 0.524802488530383 1.84319767102541 -0.585338467826228 -1.04949732751729 -2.08736181130532 0.902647269364156 -0.183875793989162 -0.332287971207814 -1.20079677199443
0.280527226893870 1.28212662782898 2.31201597008432 0.215788791733861 -2.63853737908734 -2.06789229135309 1.13009566289505 0.618442903263456 0.945275965074845 1.07571804298910
0.662989335164517 0.50854786573429 0.79166072586515 0.576789098851336 -0.823748502569833 -1.94788594717133 0.880471774200295 -0.62901696365591 -1.23792837131531 -1.33862437604865
0.525981362187104 -0.660879958374254 0.410108474303663 1.56674270771706 -1.49118722414311 -2.53076555938644 -0.7275602081407 -0.0245289183187077 0.452359420214995 0.57322876599616
-0.375915355649983 -0.379179721394808 0.646418051187915 -0.876335594968393 -0.890043101500624 -2.4181570127451 1.39484824563447 0.49615705090624 0.759978145080908 -0.138541057759067
-0.715349037019646 1.42108256657665 1.41577056186057 0.584678021159656 -0.689817912559597 -0.86558423616452 -0.574316294601922 -0.315784694596942 1.07050093368525 -0.266912033245575
1.16829795228854 1.44582257663954 0.782118059890904 1.00897667417184 -1.41027106735093 0.531870820918197 -0.486423552264882 -0.77285014558925 -0.59727082959116 -0.476740941644508
1.55892936505651 1.64951343247263 -0.272900648421469 0.336445897917337 0.410382807262394 -0.126502778639833 -0.963109757413506 0.73744221075116 0.684365293157238 0.00581536505251135
-0.69113433050422 -0.48954953133724 2.52161157416577 -1.81160146672612 0.916639345904796 -0.853679891460491 -0.876375653323389 0.272726500200932 -0.418868996203363 1.17957564796320
-0.084140353963087 0.460369126752912 -0.500758900563936 -2.52211741198144 0.410186054175025 0.130828278346727 -2.37680031270358 0.836249040118725 -1.52067934317415 -1.14968663051837
0.182367297293947 -0.0140856508936454 -1.00741344907048 -2.30920766892225 1.68273259690964 0.815980720100394 -1.42516774623716 -0.331216456281788 -0.0472810201325081 -0.436686418506921
0.571653948498554 0.800727561718109 -0.958358541603888 -1.82992497500018 0.658527263503354 -0.600235901317019 -1.52825041888485 -0.555556709377448 0.777169547521236 0.480479577228766
1.94768037366047 0.302838093218334 -0.636427181869509 -1.37528582405582 0.704583081588336 -0.276735534770669 -0.965354528648934 0.788599840414772 0.336507254284256 -1.48050157276964
2.34943420754662 -0.959473099580338 -0.278467518187692 -2.10878103033985 0.807448338531135 1.13479445640072 -0.933247342306544 -0.325232471501058 0.635479964430165 -3.10630457875390
-0.355006243930992 0.366185083897958 -0.640120794041414 -1.77063266560817 0.47719632251358 -0.033323338562499 -1.63390065799813 0.349677023641206 0.312369634710555 -2.14851777822185
1.06562948974201 1.10255253965448 0.816066252311692 -1.95965927887294 0.392980034037713 -0.364392389809701 -1.22125020631228 0.96625542821032 0.510734717993895 -0.0317415333979301
0.236394644921073 2.46662180322828 -0.0744415869794826 -1.60252249286453 -0.421837564516261 -0.155529909401838 -0.841786549606316 0.819392534097112 0.0183275642870916 -1.62980622875521
0.99399205070982 0.042061022912102 0.208689866822883 0.814353566935894 0.957759882980759 0.059036171248832 -0.0138008457355067 0.603223275925629 1.89310520011297 -1.39637570884192
0.508104364561355 1.34324739647198 0.562168659263243 0.959986632300074 -0.416284941974367 0.274786904777445 0.770427738511061 -0.755229690173805 0.852372838302998 -2.98785346496334
-0.0318810835877677 0.458336612443407 -0.447800269439252 0.873271966537103 -0.591516016866593 1.53303130420884 0.495919288469357 0.524258959047156 0.748814300078823 -1.70271855660717
0.542518800161961 -0.701901957680536 -0.0489747658609474 -1.09869727991547 0.840353166533945 0.487302698952012 -2.31614433077671 0.48120800179617 0.320230571623593 1.66010344114839
1.38173934906370 0.748714063148685 1.02327771397903 1.75586867846493 -0.334991644086096 0.979440009053994 -1.90370197125261 0.36311248414635 0.659886958730137 -0.73783624321504
1.23387992010073 0.796003711746704 1.22186657041013 1.37339497086423 0.968905910033959 2.08952401065463 -1.86361575039155 1.90029661886147 1.52999395683082 -1.8963832359471
0.39349465621243 0.545130758995853 -0.381718343601273 0.763810090359947 0.722432474084322 1.41976239262509 -0.0410001377831185 0.252175436874292 1.20999391121053 -2.52439417840336
-0.221221936936181 0.585454626168788 -0.333915358119306 -0.943806280542668 0.996020691071955 -0.406873753136822 -1.52130986610332 1.2335347386179 0.576979528181158 -0.0421061907342301
-0.33094377146155 -0.765890667451416 -1.12528176901212 -0.990647465066861 1.58875127137006 0.892230658745761 -0.653265987430149 1.31675889988689 0.681492300213563 1.13071507345378
0.149029932347682 0.820402014458334 -1.22871696733561 -0.0393788936625018 -0.106845805627417 1.68123994038917 -2.22645197605410 0.549903854033006 1.97579147819074 -1.74680971732985
0.88885972816151 -0.212994646086827 0.828654471396416 0.36417197805615 0.264989079483492 2.00536314356403 -0.625304775546514 2.24684467909164 0.170203486578434 -1.34400545594835
0.942008877727461 0.0955959253869856 -0.407291160341028 -1.74109109471500 0.911885420092523 -0.00632093055219914 -0.411354885645284 0.740790490539641 -0.670678096477916 -0.524244376651373
-0.282483005721801 0.174045503735345 -0.92900410187046 -1.00981752053125 -0.976080084928988 -0.435079716478918 -1.57551089792273 0.493103677265778 -0.654471101768679 -0.364293466412593
0.659608973820157 0.533731200421421 -1.53797590606662 0.230065686697437 0.296170826353167 1.41389241402524 0.113845206894840 0.611097024644335 0.743407836661837 -0.491998331554214
0.190447676128778 0.368448565615521 0.773773310393281 0.761873480011786 2.04638321626921 0.755326876503169 1.32473118091056 0.812456037400774 0.275738505786561 -0.835464770988384
-0.000880528579599139 0.284786506078641 -0.790829381550605 -0.418510937517212 -0.534066383786573 -0.56170304279466 -0.0619790406225042 1.41517596697113 -0.631698552023088 0.483572682470558
-0.745930605850684 -0.432569677399161 -2.26135817881265 0.851426038541 0.292965139051927 -0.190280980501961 0.34595173702483 -0.217204003000961 -0.344597321214711 0.765446966926517
0.539635946435506 -0.565720116972267 0.213122586163713 0.813378924817916 -0.307301061460563 -0.502423152370049 0.50846089158298 1.53613324560675 -1.07663141379610 -0.0883933139452832
-0.689522160529188 0.401464281580526 -0.242073505083123 1.38276820480909 -0.59889556880917 0.410478452747694 1.34297373920722 0.625518418442898 -1.56892698292472 -1.29856791122098
-1.16228762789092 -0.165287790292222 -1.01790202935334 -0.0677933853613085 -0.438766340413191 -0.520046449976325 1.49303100276316 1.18877035953818 -0.165745014838455 -0.2515166077663
-1.56418301602139 0.295434992981415 -0.0401633940352366 0.835265051031095 0.823372643544969 -0.603271647499787 -0.368875193828901 0.615528030483828 -0.306913601051777 0.956938160883928
-1.57830149459182 1.38965184610686 -1.16082646391755 1.03902149662082 0.0422511327785613 -0.351457779817942 1.84190347397384 0.568153119177459 0.738170385055213 0.0646518866145899
-0.705203061231291 1.47091633350400 0.463991764776935 1.47756874371350 -0.500733257975353 0.54829423349723 0.935398054725603 -0.89398619538804 0.184539589801195 0.902443915992114
-1.28414886739536 1.08690274200725 -0.233691190146298 0.226869949632643 -0.572438074597811 -1.16530985870940 -3.55599119012488 -1.45982057936757 -0.229512046159207 1.95616861828764
-0.269641149970506 -0.482874804819613 1.31109361792897 -0.686148113695446 0.0879075520838048 0.765392798122279 -1.26722040218017 -0.633204797639794 -0.0923234495588488 1.64865173989971
-1.36954254637995 0.766716062046485 1.45469539382638 2.06894864362377 -0.358719864989316 0.0249826850688207 -2.15230892892594 -0.127803847965862 0.915913414430536 1.81293413091805
-1.31089858705424 1.06871063928854 -1.12003765286769 0.496635005047796 -0.46454585164694 0.348022841217322 -1.10476371096586 -0.493873768935972 0.671995371397307 0.264229799637037
-0.743669608636617 0.803686110450666 -0.962817289986276 -1.28116093020535 -1.75416890477329 1.47901763750443 -3.66716072477172 -1.90740601905927 -0.259194249571728 1.27667954830418
-1.34317766958443 -0.965706560865458 0.0557956590957696 -0.756571885703914 1.21665852468195 1.44505534916850 -1.58095504875690 -0.828278504046493 0.4625264235579 1.66196558667910
-0.198658674737569 -1.07436480907111 -0.372714168774302 -0.114869841137163 -0.17845883973048 0.124819111825016 -0.587621255493976 -0.035634776160839 -0.109925681795669 -0.441785884730637
0.256748724589175 0.507703993560952 0.337095304630917 0.508903702093219 0.43325143167282 -0.40367500189962 -0.320843391980982 -0.728513230539846 -1.06735016166725 -0.393746320188533
-1.92497138122766 -0.380838431825787 1.06197947662071 -0.601587262073867 0.309430720477480 -0.498855486389624 -0.328910435052999 0.440712004993783 0.384352324621054 1.54984383926721
-0.629818660207998 1.10089876719004 0.121282959297534 -1.32252025816745 -0.512898226602611 -0.585698339317692 0.629097286072275 -1.57895716765022 -0.470516347326203 1.11835076224572
0.266854462116261 1.21396376906626 0.504224142296076 -1.75524452029695 0.32707731326942 -1.21622618359683 -0.744986401076596 -0.756113717404717 -1.32399608226916 0.90281252652788
0.247341313835710 1.05004592775391 0.6320728603509 -1.86383367945858 0.269462064297585 -1.13953408583799 -1.18716389119965 -0.661149649568322 0.444493578416841 -0.320513332792246
0.936080490321582 1.00813529320429 0.203044510546418 -2.18156761621396 0.955345926239357 0.756408584638223 -0.862399865218856 -0.407649947175573 -0.479162805748645 -0.749452204779827
0.0134441340172317 -0.46518676683434 0.199913901255210 -1.81207249694172 0.564308262381411 -0.572697675409941 -0.969434319591319 0.6175909533139 -1.52594679850368 0.00396190955886777
0.475817003175057 0.701800290724853 0.302544210235209 -2.87037073018937 0.743732774735474 -0.0667044727487769 -0.819680316602622 0.114357807207163 -0.808024206368864 -0.500371386542032
-0.674776287526495 0.655223198067882 0.092422244863333 -0.774375244309869 1.17410234675765 -0.974204483175368 1.17563734131020 1.07555814480442 -0.718506498269765 0.585749946054376
-0.395188024041528 -0.207132762082066 -0.120561610942403 -2.95004046654932 0.883289341842102 -1.08459965242902 0.746181055310485 -0.193341528330582 -0.75662129506846 -0.465847298623625
0.0196252728315817 0.886687887539179 -0.327459847858564 -0.816353848271954 0.153045850891014 -2.70981622912668 0.95911038559281 -1.9023480836967 -0.860227906026822 0.397482132991429
0.203409853324835 -0.256225072875711 0.91957915068254 0.180670031909999 1.32296743668724 -1.60505458629745 -0.388636633223224 -0.282973092404104 0.0282345455449358 0.0890798693346371
-0.126649786501909 0.314012801065664 0.433020739994493 -0.0181443310966560 1.63489943035167 -0.75249162317383 1.42944007012108 1.15227890294559 -0.882995039020258 0.374062714838295
-0.0269094253542966 -1.76248581027535 0.511388024941482 -0.412655188844004 -0.301796646081863 -0.526196366518792 1.49592000218619 1.51883975119453 -0.257514804500963 0.210306391048620
0.49812551693369 -1.55062389608717 0.128774855387068 0.0541082441767594 -0.322080694337649 -2.30097348671336 0.0735897445427817 1.27427964329434 -0.50023823698614 -0.307209664683821
-0.20683291398962 -1.98464013668522 0.812976234240209 -0.235573574042146 -0.367713237644013 -1.28411323452171 1.09917509499075 -0.375716934374667 -0.721897589306043 -0.0821591664183303
0.227679973800059 -0.285937422761897 0.660540938026578 1.02086841835877 0.93976520636416 -0.926517220364392 0.045060276151561 1.00589423307672 0.0946127590227907 0.832170574578696
0.237428731444062 -1.74834169628753 -0.191029396129022 0.380289605931682 0.793590921718829 -0.242454079389901 -1.04329402456887 1.30776927273839 0.199698284616362 0.0943070232050067
0.781884978491733 -1.59015754972931 0.763544243740902 -1.05321429279068 0.258189172764122 -1.49746590482746 -0.809863060913629 -0.137281234914450 1.18885053995576 -0.625082056144578
1.61183092470604 0.162266849521403 -0.111208841023641 0.0566221281968713 0.428336093775853 0.68561588445816 -0.748630910713866 0.404976312552999 2.25155167976550 0.179433142942350
1.01964126262773 -0.083402058775627 0.647313751689937 0.521636400969847 1.82811453228243 0.818140428439613 -1.14719870196744 0.63496537066156 1.63823150827067 0.153272917430012
1.51839706060119 -0.456813179266717 -0.41513525695058 -0.582473825750797 -0.46572844801414 -0.017888452694537 0.537502798286846 0.448811936036691 -0.87012587275299 -0.256349125619634
-0.289267137483731 -1.20436067607428 0.873970882752987 -2.01012984369007 1.16983817386097 -0.127782476526275 -0.386024240853474 1.37843512143946 0.89377651745078 -0.72044231047993
0.781544157559182 0.903590279543391 1.19330296951697 -0.993464382447492 0.624169536522855 -0.333254462956691 -1.04338530975646 0.872193387020996 -0.330411786361454 -0.928399602419626
0.453601434478319 0.824530094274635 -0.639615949045025 -0.581001534661772 1.52606367168691 0.671295717523005 -1.09627507763335 1.4263451813256 1.80102245912626 -0.955681598851364
-0.62559754295568 1.10756532735480 0.383037217982879 -0.415649844127594 -0.974398462177556 0.803710743861091 -0.987074437048068 0.483467695438718 1.04150085286466 -0.946784507189405
-0.340053104945521 1.35908928713120 0.558893997801316 -0.181880190164089 -1.32819714058955 0.169423223442671 -0.657469632498577 0.0240371201855713 -1.05441704453208 -0.614295560602156
0.376662011857016 1.09497148641317 1.16096898522823 -0.94144097248625 0.894149822846047 -0.09428000172514 -1.84774070672963 -0.32226571698513 -0.557058644665565 -0.0727297718581722
0.370146677679128 -0.262185232034980 1.87436380692619 0.236546751795295 -0.139137920730192 0.714940090107123 0.142188755326245 -0.464424965710215 0.642215081939469 -1.32488036218273
-0.674576357608925 -0.490190598702266 0.46528561801715 -0.515019070712407 -0.0416611628863148 -0.99788201296753 0.141767569700798 -1.02324526343877 -0.112290220739643 -0.334019599671732];
                                                                                                                                                                                                                                     kk.m                                                                                                0000664 0002715 0000074 00000003146 10362007434 010322  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% interpolate to intpoints
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
                                                                                                                                                                                                                                                                                                                                                                                                                          kksim.m                                                                                             0000664 0002715 0000074 00000002776 10350120454 011036  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   orig=[0.95   91.31    0.53 1988.85
1.11  111.84    0.49 2275.81
0.98   91.80    0.55 2097.17
0.94   93.13    0.50 2128.73
1.05  107.56    0.48 2205.36
1.01   96.93    0.52 2218.73
0.94   85.69    0.55 2181.16
0.98  102.05    0.47 2127.77
0.91   82.81    0.55 2120.48
1.01  101.11    0.50 2173.30
0.92   86.00    0.54 2035.08
1.07  109.61    0.50 2083.94
1.06  107.87    0.50 2074.91
1.03  106.23    0.49 2097.01
0.91   79.04    0.58 2203.51
1.03  111.58    0.47 2020.61
1.02  103.50    0.50 2085.31
0.97   94.73    0.53 1978.17
0.96   91.40    0.53 2093.31
0.95   92.43    0.53 1964.16];

resi=[1.05  126.09    0.44 1822.98
0.89   87.80    0.53 1798.51
0.91   91.49    0.51 1880.34
0.97   99.42    0.50 1960.72
0.94   95.23    0.50 1896.43
1.19  144.50    0.46 1763.14
1.10  123.35    0.48 1850.74
0.93   90.70    0.52 1978.43
0.88   86.74    0.52 1898.19
0.97   98.67    0.51 1915.25
1.05  106.27    0.52 1946.46
1.18  140.60    0.46 1848.80
0.99  100.74    0.50 1979.93
1.03  109.98    0.49 1948.63
1.13  129.48    0.47 1920.01
1.05  116.62    0.46 1971.18
0.93   97.69    0.50 1742.26
0.97  100.42    0.51 1854.26
0.96   93.37    0.51 2100.56
1.03  114.62    0.47 1928.21];



>> mean(resi(:,1))    1.0075
>> std(resi(:,1))    0.0906

>> mean(resi(:,2))  107.6890
>> std(resi(:,2))   17.3638

>> mean(resi(:,3))    0.4930
>> std(resi(:,3))    0.0254


>> mean(orig(:,1))    0.9900
>> std(orig(:,1))    0.0569

>> mean(orig(:,2))   97.3310
>> std(orig(:,2))    9.9836

>> mean(orig(:,3))    0.5155
>> std(orig(:,3))    0.0300

>> 

  knhankel.m                                                                                          0000664 0002715 0000074 00000000607 10340447217 011512  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   
n0=256; %n0=power(2,15);
m=2;
pp=2;
t=2;
del=t/(n0+0.0);
del1=ceil(log(pi/del)/log(2.0));
n= 4^del1; %(int) power(4, (int) del1);
if (n<n0) n=n0; end
mn=m*n;
rt = 5;
phi = rt/(mn-1)/del; 
phi=1;

xx = (0:(mn-1))*del*phi;  
%aa = kcovario(xx);

aa = exp(-xx);
aa(n0+1:end)=0;

[fw,del2]=knhankelc(n0, n,  mn, pp, del, aa);

ww = del2*(0:(mn-1))/phi;
fw = fw*phi^2;

plot(ww(1:n0),fw(1:n0))

                                                                                                                         knodos.m                                                                                            0000644 0002715 0000074 00000000066 10253635142 011210  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = knodos(wt,ll)
  res = (0:ll) / ll * wt;                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ktilde.m                                                                                            0000664 0002715 0000074 00000000663 10341165673 011201  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% calc positive definite kernel covariogram at r
%% by hankel tr, taking positive part and hankel 
%% transforming back

function [res,pp] = ktilde(r,R)
  [w,fw]=sieghankel(@kcovario,R);
  pp = spline(w,fw);
  [x,kx] = sieghankel(@fdens,1/R,pp);
  pp = spline(x,kx);
  res = ppval(pp,r);
  
  
function res = fdens(x,pp) 
  %% returns truncated spectral density
  %% corresp to kcovario
  res = ppval(pp,x);
  res(find(res<0))=0;
                                                                                 maternkernel.m                                                                                      0000664 0002715 0000074 00000003562 10373137662 012417  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                     %------------------------------------
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
                                                                                                                                              mplebsplines.m                                                                                      0000644 0002715 0000074 00000015271 10260317443 012413  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %%function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values

%------------------------------------
%% read sites locations
%% simul.dat: "lat"   "wnh4"  "lon"   "cwnh4" "site"
%------------------------------------
[lat,wnh4,lon,cwnh4]=textread('simul.dat','%f%f%f%f%*[^\n]');
M = load('M.txt');
nobs = length(wnh4);
dista = cordist([lon,lat]);
indupper = find(triu(dista,1)); 
rr = dista(indupper); rmin = min(rr); rmax = max(rr);
nnots=100; cnots = linspace(rmin,rmax,nnots);

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

ll = 6; %number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

%% calculate or read tabulated values
%% hankel transform of polynomials
%tic
%Hpol = calcHpol(ll,cnots,wtvect);
%disp('hPol calc time')
%toc
%save strcat('tablahpol_',num2str(ll)) Hpol

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


%------------------------------------
%% read observations or simulate
%------------------------------------
simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end)
simt = nodos(simwt,ll);
simknots = simt(4:end-3);
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simsig2=1;
simbcoef = [0.2 1.0 0.2 0.8 0.4 0.2];
simbcoef(end)=contderiv(simbcoef,simgm)
simtheta = [simsmoo simsig2 simbcoef simwt];
randn('state',0);
%simcova = calcCova(simtheta);

%% test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
simrango = 500; simirango = 2*sqrt(simsmoo)/simrango;
u = 0.5*simirango; v=.1*simirango;
simcova = polmatern(dista,u,v,simirango,simsmoo);
simwt = 1.5*simirango;

%% test matern cova
%simrango = 500; simirango = 2*sqrt(simsmoo)/simrango;
%simcova = hmatern(dista,simrango,simsmoo);

for ii=1:nt
  Za(:,ii) = simulachol(simcova*simsig2);
end

cc = chol(simcova);
invcc = inv(cc);
xx = invcc' * Za;
nobs = size(Za,1);
xx=reshape(xx,nt*nobs,1);
truene = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
disp(truene)

%------------------------------------
%% def initial f
%------------------------------------
indi = find(smoovect>=0.2); smoo=smoovect(indi(1));%smoo=.5;
gm=2*(smoo+1);sig2=3;
bcoef = ones(1,ll+3);
wt = wtvect(wtvect <= simwt*1); wt = wt(end)

t = nodos(wt,ll);
theta = [smoo,sig2,bcoef,wt];thetaini=theta;
Ei = ener(theta); Eini = Ei;


%------------------------------------
%% improve initial guess fitting matern
%------------------------------------
matheta = log([smoo sig2 500]);
%[paramf fval exitflag] = fminunc(@enermatern,matheta)
[paramf fval exitflag] = fminsearch(@enermatern,matheta);
matheta = exp(paramf);
matsmoo = smoovect(smoovect<=matheta(1));matsmoo=matsmoo(end);
matsig2 = matheta(2);
matirango = 2*sqrt(matsmoo)/matheta(3);
matwt = wtvect(wtvect <= matirango*.75); matwt = matwt(end);

theta = [matsmoo matsig2 bcoef matwt ];
matEi = ener(theta);

%------------------------------------
%% loop for simulated annealing
%------------------------------------
%% minimize (-l(f)+penalty)
L = 20000; %% max number of iterations
Ti = 1; %% scale temperature
lam = 0.995;
tsca = Ti * lam.^(1:L);

fixwt = false;
conta=1;
llvect = nan(L,1);
for T=tsca
  % generate new state %logf m-vector
  thetaf = bproposal(theta,fixwt);  %% bproposal.m generates f and proposal.m log(f)
  %% evaluate energy 
  Ef = ener(thetaf);
  %disp(Ef)
  % accept or reject
  dE = Ef-Ei;
  if dE < 0
    disp(conta)
    disp('descent')
    theta = thetaf;Ei = Ef;%disp(theta)
    llvect(conta) = - Ef;
  else
    p = exp(-dE/T)/(1+exp(-dE/T));
    pu = rand(1);
    if pu<p
      theta = thetaf; Ei = Ef; %disp(theta)
      llvect(conta)=-Ef;
      disp('ascent')
    end
  end
  conta=conta+1;
end
thetaf = theta;
Ef = Ei;

%save matlaboutput simtheta truene thetaf Ef;

smoof=thetaf(1); gmf=2*smoof+2;
sig2f=thetaf(2);
bcoeff=thetaf(3:end-1);
wtf = thetaf(end);
tf = nodos(wtf,ll)

disp('number of indep replicates')
disp(nt)
disp('sim  ini  after-matern  final')
disp([simsmoo smoo matsmoo smoof ])
disp([simsig2 sig2 matsig2 sig2f])
disp([simwt wt matwt wtf])
disp(round([truene Eini matEi Ei]))
disp('estimated ener - true ener'); disp(Ei - truene)
disp('iter');disp(L)
disp('lam -- cooling schedule T=lam^(1:iter))');disp(lam)

figure(1)
simf0 = spectralfn(simwt,simsmoo,simt,simbcoef);
ff0=spectralfn(simwt,smoof,tf,bcoeff);
xx = linspace(0,2*simwt,100);
plot(xx,spectralfn(xx,simsmoo,simt,simbcoef)/simf0);hold on
plot(xx,spectralfn(xx,smoof,tf,bcoeff)/ff0,'r.');hold off
grid
xlabel('freq')
ylabel('f(w)')
title('freq')
legend('simulated','estimated')

figure(2)
covaf = calcCova(theta);
%simcova = hmatern(dista,simrango,simsmoo);
%plot(covaf(:),simcova(:),'.')
%hold on; plot([0 max(covaf(:))],[0 max(covaf(:))]);hold off
eigf = eig(covaf);simeig=eig(simcova);
plot(eig(covaf),eig(simcova),'.')
hold on; plot([0 max(eigf)],[0 max(eigf)]);hold off
xlabel('eig(cova) estimated')
ylabel('eig(cova) simulated')
title('eig of cova -- simulated vs. estiamted')
grid

figure(3);
%length( llvect( find(~isnan(llvect)) ) );
%evol = cumsum(~isnan(llvect));
%plot(evol); ylabel('cumulated number of changes in ener')
%xlabel('iteration number');title('evolution of energy')
[covacnots,cova0] = bhankelx(cnots,smoof,bcoeff,wtf);
plot(cnots,polmatern(cnots,u,v,simirango,simsmoo))
hold on; plot(cnots,covacnots,'r.')
grid
hold off
xlabel('distance');title('cov function')
legend('simulated','estimated')

figure(1);print -depsc plotspectral.eps
figure(2);print -depsc ploteig.eps
figure(3);print -depsc plotcovafn.eps

%[thetaff,fval,exitflag]=fminunc(@ener,thetaf);

%% TEST MATERN or polmatern
figure(1)
%simf0 = fmatern(0,simirango,simsmoo);
%simf0 = fpolmatern(0,u,v,simirango,simsmoo);
%ff0=spectralfn(0,smoof,tf,bcoeff);
xx = linspace(0,2*simwt,100);
%plot(xx,fmatern(xx,simirango,simsmoo)/simf0);hold on
plot(xx,fpolmatern(xx,u,v,simirango,simsmoo));hold on
plot(xx,spectralfn(xx,smoof,tf,bcoeff)/cova0,'r.');hold off
grid
xlabel('freq')
ylabel('f(w)')
title('freq')
legend('simulated','estimated')

figure(1);print -depsc plotspectral.eps
figure(2);print -depsc ploteig.eps
figure(3);print -depsc plotcovafn.eps


poltheta = log([smoof sig2f 500 .1 .1]);
[kkf fval exitflag] = fminsearch(@enerpolmatern,poltheta)
uf = exp(kkf(4));
vf=exp(kkf(end));
smf = exp(kkf(1));
af = 2*sqrt(smf)/exp(kkf(3));
figure(1);hold on;
plot(xx,fpolmatern(xx,uf,vf,af,smf),'k.');hold off                                                                                                                                                                                                                                                                                                                                       myhankel.m                                                                                          0000645 0002715 0000074 00000130273 10320251537 011525  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %
%   based on Brian Borchers' code but modified so that
%   H(f(w))=Int[f[w] w besselj(0,w),{w,0,inf}] and vectorized
%   vectorized by eats up too much memory even with n=100
%
%  z0=hankel0(funname,B,varargin)
%
%    funname       The name of the function to be transformed.
%    B             The transform argument.
%    z0            The value of the transform of order 0.
%    z1            The value of the transform of order 1.
%
%  Based on Walt Anderson's Fortran code, which was published in:
%
%  Anderson, W. L., 1979, Computer Program Numerical Integration of Related
%  Hankel Transforms of Orders 0 and 1 by Adaptive Digital Filtering.
%  Geophysic, 44(7):1287-1305.
%
function z0=hankel0(funname,B,varargin)
%
%  Points at which the function will be evaluated.  Note that these
%  are scaled by B before being used.
%
YBASE=[
   8.9170998013274418e-14
   9.8549193740052245e-14
   1.0891370292130841e-13
   1.2036825704856076e-13
   1.3302749714952345e-13
   1.4701812115404443e-13
   1.6248015192957209e-13
   1.7956833867707590e-13
   1.9845370571306282e-13
   2.1932526413842005e-13
   2.4239190352504162e-13
   2.6788448255287407e-13
   2.9605813952117967e-13
   3.2719484585839035e-13
   3.6160622818693727e-13
   3.9963668718722961e-13
   4.4166684447542096e-13
   4.8811735199247527e-13
   5.3945310203017795e-13
   5.9618788002944785e-13
   6.5888950671771899e-13
   7.2818552104963217e-13
   8.0476946082781583e-13
   8.8940780386232124e-13
   9.8294763913816720e-13
   1.0863251447666186e-12
   1.2005749575703847e-12
   1.3268405280766938e-12
   1.4663855645544968e-12
   1.6206066806315702e-12
   1.7910473730731204e-12
   1.9794134696161976e-12
   2.1875902014670364e-12
   2.4176610713286156e-12
   2.6719287057960000e-12
   2.9529379008172425e-12
   3.2635010908665675e-12
   3.6067264967338821e-12
   3.9860492336431492e-12
   4.4052656910401313e-12
   4.8685715281339743e-12
   5.3806036654647833e-12
   5.9464866927629101e-12
   6.5718841575654077e-12
   7.2630552479033658e-12
   8.0269174363595134e-12
   8.8711157124588670e-12
   9.8040990962934688e-12
   1.0835205199155280e-11
   1.1974753677488471e-11
   1.3234149515479672e-11
   1.4625997169973058e-11
   1.6164226720110948e-11
   1.7864233284247930e-11
   1.9743031099469828e-11
   2.1819423805797136e-11
   2.4114192639334462e-11
   2.6650304417866294e-11
   2.9453141400488783e-11
   3.2550755321790057e-11
   3.5974148143038491e-11
   3.9757582330231206e-11
   4.3938923764369768e-11
   4.8560020715924429e-11
   5.3667122676390673e-11
   5.9311343238745088e-11
   6.5549171659463761e-11
   7.2443038221987792e-11
   8.0061939059983487e-11
   8.8482126693838504e-11
   9.7787873191515278e-11
   1.0807231359173195e-10
   1.1943837803073370e-10
   1.3199982190169224e-10
   1.4588236435691521e-10
   1.6122494654737812e-10
   1.7818112219246310e-10
   1.9692059439719362e-10
   2.1763091409794871e-10
   2.4051935713527238e-10
   2.6581499874015355e-10
   2.9377100619593265e-10
   3.2466717262156566e-10
   3.5881271723520049e-10
   3.9654938012404428e-10
   4.3825484249401899e-10
   4.8434650663021337e-10
   5.3528567339924565e-10
   5.9158215910338558e-10
   6.5379939789346253e-10
   7.2256008080722362e-10
   7.9855238787053352e-10
   8.8253687563437822e-10
   9.7535408908045947e-10
   1.0779329740778886e-09
   1.1913001745856734e-09
   1.3165903076505279e-09
   1.4550573190356333e-09
   1.6080870331313015e-09
   1.7772110227512649e-09
   1.9641219376281762e-09
   2.1706904450210516e-09
   2.3989839519819520e-09
   2.6512872966606394e-09
   2.9301256157327411e-09
   3.2382896168163261e-09
   3.5788635088117363e-09
   3.9552558697009004e-09
   4.3712337607414382e-09
   4.8309604284818817e-09
   5.3390369719324459e-09
   5.9005483919104070e-09
   6.5211144834374116e-09
   7.2069460805369272e-09
   7.9649072163486862e-09
   8.8025838206794284e-09
   9.7283596425381269e-09
   1.0751500157513941e-08
   1.1882245299770154e-08
   1.3131911946747030e-08
   1.4513007182274982e-08
   1.6039353471673310e-08
   1.7726227001629019e-08
   1.9590510569407677e-08
   2.1650862551562963e-08
   2.3927903643240499e-08
   2.6444423237025737e-08
   2.9225607506844726e-08
   3.2299291479658129e-08
   3.5696237617766722e-08
   3.9450443699873716e-08
   4.3599483082281089e-08
   4.8184880745668262e-08
   5.3252528891055793e-08
   5.8853146244378082e-08
   6.5042785666539686e-08
   7.1883395149287240e-08
   7.9443437811532343e-08
   8.7798577101256825e-08
   9.7032434060731539e-08
   1.0723742423401342e-07
   1.1851568259277232e-07
   1.3098008573741623e-07
   1.4475538160404734e-07
   1.5997943798373573e-07
   1.7680462234971138e-07
   1.9539932680224871e-07
   2.1594965339340473e-07
   2.3866127669890702e-07
   2.6376150227843728e-07
   2.9150154162607256e-07
   3.2215902637935327e-07
   3.5604078695002664e-07
   3.9348592338593700e-07
   4.3486919919827997e-07
   4.8060479212078480e-07
   5.3115043933968360e-07
   5.8701201868132178e-07
   6.4874861160747570e-07
   7.1697809869053571e-07
   7.9238334356995170e-07
   8.7571902728105488e-07
   9.6781920135651654e-07
   1.0696056352944216e-06
   1.1820970419372223e-06
   1.3064192730922673e-06
   1.4438165874351013e-06
   1.5956641034684996e-06
   1.7634815621706371e-06
   1.9489485370736003e-06
   2.1539212439998211e-06
   2.3804511186939235e-06
   2.6308053482811660e-06
   2.9074895620382204e-06
   3.2132729085731429e-06
   3.5512157703953872e-06
   3.9247003932525885e-06
   4.3374647367828189e-06
   4.7936398852710157e-06
   5.2977913929290109e-06
   5.8549649774966190e-06
   6.4707370194807025e-06
   7.1512703724455683e-06
   7.9033760429228480e-06
   8.7345813572541237e-06
   9.6532052976029776e-06
   1.0668441761124589e-05
   1.1790451575578641e-05
   1.3030464192308714e-05
   1.4400890074365673e-05
   1.5915444904593194e-05
   1.7589286856791647e-05
   1.9439168303816350e-05
   2.1483603480955747e-05
   2.3743053782621043e-05
   2.6240132546858776e-05
   2.8999831377238596e-05
   3.2049770267221755e-05
   3.5420474030339066e-05
   3.9145677802784465e-05
   4.3262664675996812e-05
   4.7812638838370288e-05
   5.2841137960621061e-05
   5.8398488952101537e-05
   6.4540311649424626e-05
   7.1328075478483034e-05
   7.8829714661124184e-05
   8.7120308123675967e-05
   9.6282830912076270e-05
   1.0640898463402168e-04
   1.1760011523947924e-04
   1.2996822732501724e-04
   1.4363710511345378e-04
   1.5874355132796403e-04
   1.7543875635971472e-04
   1.9388981143211579e-04
   2.1428138090594562e-04
   2.3681755046234149e-04
   2.6172386966089197e-04
   2.8924960931543914e-04
   3.1967025628016629e-04
   3.5329027061462893e-04
   3.9044613272236350e-04
   4.3150971095986065e-04
   4.7689198342006658e-04
   5.2704715113927156e-04
   5.8247718389374341e-04
   6.4373684408196633e-04
   7.1143923897318674e-04
   7.8626195689103691e-04
   8.6895384874522254e-04
   9.6034252278312515e-04
   1.0613426275713101e-03
   1.1729650061058051e-03
   1.2963268126685603e-03
   1.4326626936829910e-03
   1.5833371444703617e-03
   1.7498581655775839e-03
   1.9338923553535471e-03
   2.1372815898255564e-03
   2.3620614568136901e-03
   2.6104816287778878e-03
   2.8850283782960702e-03
   3.1884494615157647e-03
   3.5237816186211822e-03
   3.8943809665496639e-03
   4.3039565881380203e-03
   4.7566076538702283e-03
   5.2568644477534125e-03
   5.8097337079228714e-03
   6.4207487357601564e-03
   7.0960247750331065e-03
   7.8423202153108801e-03
   8.6671042321983371e-03
   9.5786315413559676e-03
   1.0586025014468731e-02
   1.1699366984012178e-02
   1.2929800150624660e-02
   1.4289639103000504e-02
   1.5792493566432742e-02
   1.7453404613518231e-02
   1.9288995200267688e-02
   2.1317636534236600e-02
   2.3559631939745231e-02
   2.6037420060372587e-02
   2.8775799432443259e-02
   3.1802176677114019e-02
   3.5146840795050052e-02
   3.8843266308924096e-02
   4.2928448287690518e-02
   4.7443272605669898e-02
   5.2432925142121424e-02
   5.7947344016710048e-02
   6.4041719386992838e-02
   7.0777045810065886e-02
   7.8220732696592687e-02
   8.6447278966843177e-02
   9.5539018660927705e-02
   1.0558694496554391e-01
   1.1669162090437306e-01
   1.2896418580662142e-01
   1.4252746762678220e-01
   1.5751721224808804e-01
   1.7408344207293611e-01
   1.9239195749751564e-01
   2.1262599629790035e-01
   2.3498806753529980e-01
   2.5970197833480957e-01
   2.8701507382234348e-01
   3.1720071263778915e-01
   3.5056100280015512e-01
   3.8742982530616715e-01
   4.2817617572350458e-01
   4.7320785722246539e-01
   5.2297556200716211e-01
   5.7797738199458315e-01
   6.3876379388591276e-01
   7.0594316852237804e-01
   7.8018785966510817e-01
   8.6224093313756223e-01
   9.5292360367804285e-01
   1.0531434539328173e+00
   1.1639035178482904e+00
   1.2863123193718711e+00
   1.4215949669322265e+00
   1.5711054147362089e+00
   1.7363400135976372e+00
   1.9189524869191834e+00
   2.1207704817120212e+00
   2.3438138603014083e+00
   2.5903149157877352e+00
   2.8627407135861755e+00
   3.1638177826465683e+00
   3.4965594034715681e+00
   3.8642957660407120e+00
   4.2707072994710522e+00
   4.7198615069887930e+00
   5.2162536748687147e+00
   5.7648518627701284e+00
   6.3711466257477705e+00
   7.0412059655722290e+00
   7.7817360613311877e+00
   8.6001483871237632e+00
   9.5046338885843706e+00
   1.0504244960619703e+01
   1.1608986046819572e+01
   1.2829913767290970e+01
   1.4179247577028352e+01
   1.5670492062326328e+01
   1.7318572099218340e+01
   1.9139982226652432e+01
   2.1152951729381048e+01
   2.3377627082769912e+01
   2.5836273585494951e+01
   2.8553498198135060e+01
   3.1556495817904278e+01
   3.4875321454323611e+01
   3.8543191029858157e+01
   4.2596813816033411e+01
   4.7076759832163077e+01
   5.2027865883738443e+01
   5.7499684304247886e+01
   6.3546978891585546e+01
   7.0230273002547406e+01
   7.7616455290928698e+01
   8.5779449151653139e+01
   9.4800952570955843e+01
   1.0477125578728921e+02
   1.1579014494637693e+02
   1.2796790079449971e+02
   1.4142640240527066e+02
   1.5630034698636896e+02
   1.7273859797446769e+02
   1.9090567491054267e+02
   2.1098340000673559e+02
   2.3317271788416559e+02
   2.5769570669423729e+02
   2.8479780075142304e+02
   3.1475024692237560e+02
   3.4785281935573863e+02
   3.8443681972258412e+02
   4.2486839299489054e+02
   4.6955219194748827e+02
   5.1893542705903837e+02
   5.7351234234481581e+02
   6.3382916191693528e+02
   7.0048955677885772e+02
   7.7416068656769369e+02
   8.5557987671209173e+02
   9.4556199783295187e+02
   1.0450076212424869e+03
   1.1549120321646080e+03
   1.2763751908839718e+03
   1.4106127415182191e+03
   1.5589681785928965e+03
   1.7229262931862318e+03
   1.9041280332173003e+03
   2.1043869266043412e+03
   2.3257072316617105e+03
   2.5703039963907459e+03
   2.8406252274246667e+03
   3.1393763905017645e+03
   3.4695474876758481e+03
   3.8344429822617740e+03
   4.2377148710149695e+03
   4.6833992345424385e+03
   5.1759566317540530e+03
   5.7203167426353639e+03
   6.3219277061418234e+03
   6.9868106470046323e+03
   7.7216199371708199e+03
   8.5337097949943000e+03
   9.4312078887249972e+03
   1.0423096680944496e+04
   1.1519303328070666e+04
   1.2730799034675721e+04
   1.4069708856989137e+04
   1.5549433054535755e+04
   1.7184781204437102e+04
   1.8992120420636886e+04
   2.0989539161478522e+04
   2.3197028265075980e+04
   2.5636681024340771e+04
   2.8332914304083228e+04
   3.1312712913202311e+04
   3.4605899677722984e+04
   3.8245433917662871e+04
   4.2267741314984989e+04
   4.6713078474065944e+04
   5.1625935823323234e+04
   5.7055482890376610e+04
   6.3056060407206925e+04
   6.9687724170466376e+04
   7.7016846100076829e+04
   8.5116778511712779e+04
   9.4068588251431182e+04
   1.0396186803991429e+05
   1.1489563314653141e+05
   1.2697931236743495e+05
   1.4033384322573253e+05
   1.5509288235486683e+05
   1.7140414317912661e+05
   1.8943087427924512e+05
   2.0935349323906592e+05
   2.3137139232536239e+05
   2.5570493407266162e+05
   2.8259765674555639e+05
   3.1231871175151330e+05
   3.4516555739862355e+05
   3.8146693595832947e+05
   4.2158616382857127e+05
   4.6592476772641251e+05
   5.1492650330238225e+05
   5.6908179639617680e+05
   6.2893265138330159e+05
   6.9507807573703467e+05
   7.6818007509655319e+05
   8.4897027884187771e+05
   9.3825726248661662e+05
   1.0369346401734781e+06
   1.1459900082649642e+06
   1.2665148295397097e+06
   1.3997153569188234e+06
   1.5469247060505589e+06
   1.7096161975797976e+06
   1.8894181026362628e+06
   2.0881299391192668e+06
   2.3077404818776865e+06
   2.5504476670371005e+06
   2.8186805896832864e+06
   3.1151238150622859e+06
   3.4427442466117009e+06
   3.8048208197275074e+06
   4.2049773184515880e+06
   4.6472186435204167e+06
   5.1359708947577253e+06
   5.6761256689692009e+06
   6.2730890166874416e+06
   6.9328355477427216e+06
   7.6619682271663100e+06
   8.4677844598838333e+06
   9.3583491255965196e+06
   1.0342575294807941e+07
   1.1430313433829404e+07
   1.2632449991557652e+07
   1.3961016354714479e+07
   1.5429309262008933e+07
   1.7052023882367507e+07
   1.8845400889123969e+07
   2.0827389002136763e+07
   2.3017824624610133e+07
   2.5438630372484632e+07
   2.8114034483345896e+07
   3.1070813300769802e+07
   3.4338559260968812e+07
   3.7949977063839935e+07
   4.1941210992593758e+07
   4.6352206657889292e+07
   5.1227110786931656e+07
   5.6614713058756173e+07
   6.2568934407734923e+07
   6.9149366682411388e+07
   7.6421869060750201e+07
   8.4459227190926239e+07
   9.3341881654555663e+07
   1.0315873304307374e+08
   1.1400803170473446e+08
   1.2599836106711893e+08
   1.3924972437657478e+08
   1.5389474573104006e+08
   1.7007999742659190e+08
   1.8796746690225038e+08
   2.0773617796471399e+08
   2.2958398251878911e+08
   2.5372954073575363e+08
   2.8041450947784531e+08
   3.0990596088136274e+08
   3.4249905530437142e+08
   3.7851999539077419e+08
   4.1832929081601185e+08
   4.6232536638906646e+08
   5.1094854962186480e+08
   5.6468547767501700e+08
   6.2407396778608418e+08
   6.8970839992525887e+08
   7.6224566554988432e+08
   8.4241174199494874e+08
   9.3100895829826319e+08
   1.0289240251791439e+09
   1.1371369095373254e+09
   1.2567306422910707e+09
   1.3889021577146211e+09
   1.5349742727587159e+09
   1.6964089262472496e+09
   1.8748218104523966e+09
   2.0719985414859231e+09
   2.2899125303454008e+09
   2.5307447334747562e+09
   2.7969054805094066e+09
   3.0910585976653914e+09
   3.4161480682074847e+09
   3.7754274968232164e+09
   4.1724926727921586e+09
   4.6113175578536234e+09
   5.0962940589514427e+09
   5.6322759839148350e+09
   6.2246276199985800e+09
   6.8792774214728651e+09
   7.6027773435862408e+09
   8.4023684167359400e+09
   9.2860532171338844e+09
   1.0262675959279177e+10
   1.1342011011829447e+10
   1.2534860722767656e+10
   1.3853163532931507e+10
   1.5310113459941998e+10
   1.6920292148366428e+10
   1.8699814807718300e+10
   2.0666491498890625e+10
   2.2840005383231522e+10
   2.5242109718238716e+10
   2.7896845571472111e+10
   3.0830782431638401e+10
   3.4073284124964363e+10
   3.7656802698239258e+10
   4.1617203209806610e+10
   4.5994122679122765e+10
   5.0831366787370079e+10
   5.6177348299437775e+10
   6.2085571595145073e+10
   6.8615168159057838e+10
   7.5831488388260895e+10
   8.3806755641097107e+10
   9.2620789072812759e+10
   1.0236180249249139e+11
   1.1312728723650484e+11
   1.2502498789457556e+11
   1.3817398065384482e+11
   1.5270586505337646e+11
   1.6876608107657605e+11
   1.8651536476342874e+11
   2.0613135691081287e+11
   2.2781038096130206e+11
   2.5176940787416525e+11
   2.7824822764365344e+11
   3.0751184919785828e+11
   3.3985315269713715e+11
   3.7559582077719836e+11
   4.1509757807371277e+11
   4.5875377145070300e+11
   5.0700132676483929e+11
   5.6032312176626892e+11
   6.1925281890144031e+11
   6.8438020638623755e+11
   7.5635710100467944e+11
   8.3590387171037695e+11
   9.2381664932114575e+11
   1.0209752944638193e+12
   1.1283522035151340e+12
   1.2470220406715007e+12
   1.3781724935494902e+12
   1.5231161599626948e+12
   1.6833036848418264e+12
   1.8603382787767620e+12
   2.0559917634869844e+12
   2.2722223048088804e+12
   2.5111940106775947e+12
   2.7752985902466255e+12
   3.0671792909169141e+12
   3.3897573528452603e+12
   3.7462612456976733e+12
   4.1402589802589175e+12
   4.5756938182836924e+12
   5.0569237379856543e+12
   5.5887650501481416e+12
   6.1765406013813154e+12
   6.8261330469601016e+12
   7.5440437264154141e+12
   8.3374577311253535e+12
   9.2143158151247129e+12
   1.0183393868840340e+13
   1.1254390751152201e+13
   1.2438025358832957e+13
   1.3746143904869607e+13
   1.5191838479344713e+13
   1.6789578079474348e+13
   1.8555353420195434e+13
   2.0506836974615496e+13
   2.2663559846063445e+13
   2.5047107241936324e+13
   2.7681334505709973e+13
   3.0592605869234598e+13
   3.3810058314828449e+13
   3.7365893187990141e+13
   4.1295698479287656e+13
   4.5638805000929469e+13
   5.0438680022752680e+13
   5.5743362307269414e+13
   6.1605942897748391e+13
   6.8085096471220516e+13
   7.5245668574367812e+13
   8.3159324619549984e+13
   9.1905267136338875e+13
   1.0157102845705528e+14
   1.1225334676977153e+14
   1.2405913430661245e+14
   1.3710654735730897e+14
   1.5152616881705944e+14
   1.6746231510403516e+14
   1.8507448052659994e+14
   2.0453893355595603e+14
   2.2605048098024984e+14
   2.4982441759638447e+14
   2.7609868095271022e+14
   3.0513623270798212e+14
   3.3722769044002506e+14
   3.7269423624413281e+14
   4.1189083123143062e+14
   4.5520976809898188e+14
   5.0308459732695450e+14
   5.5599446629754788e+14
   6.1446891476304075e+14
   6.7909317465761662e+14
   7.5051402729526438e+14
   8.2944627657455900e+14
   9.1667990297633300e+14
   1.0130879699538496e+15
   1.1196353618452902e+15
   1.2373884407605195e+15
   1.3675257190914975e+15
   1.5113496544604105e+15
   1.6702996851533248e+15
   1.8459666365023652e+15
   2.0401086424003345e+15
   2.2546687412956410e+15
   2.4917943227741685e+15
   2.7538586193560145e+15
   3.0434844586042220e+15
   3.3635705132645935e+15
   3.7173203121568085e+15
   4.1082743021675935e+15
   4.5403452822331500e+15
   5.0178575639460460e+15
   5.5455902507190850e+15
   6.1288250686585730e+15
   6.7733992278544400e+15
   7.4857638431407750e+15
   8.2730484990213790e+15
   9.1431326049478160e+15
   1.0104724255097566e+16
   1.1167447381907442e+16
   1.2341938075624136e+16
   1.3639951033870320e+16
   1.5074477206609342e+16
   1.6659873813938872e+16
   1.8412008037975264e+16
   2.0348415826945328e+16
   2.2488477400850208e+16
   2.4853611215221080e+16
   2.7467488324221096e+16
   3.0356269288511564e+16
   3.3548865998935916e+16
   3.7077231036440888e+16
   4.0976677464246272e+16
   4.5286232252850760e+16
   5.0049026875070080e+16
   5.5312728980313968e+16
   6.1130019468443072e+16
   6.7559119737921448e+16
   7.4664374385141264e+16
   8.2516895186770448e+16
   9.1195272810315088e+16
   1.0078636337593507e+17
   1.1138615774168800e+17
   1.2310074221230024e+17
   1.3604736028656150e+17
   1.5035558606966758e+17
   1.6616862109441658e+17
   1.8364472753028080e+17
   2.0295881212439261e+17
   2.2430417672705789e+17
   2.4789445292164490e+17
   2.7396574012127472e+17
   3.0277896853110349e+17
   3.3462251062551731e+17
   3.6981506727678112e+17
   4.0870885742048762e+17
   4.5169314318104928e+17
   4.9919812573787520e+17
   5.5169925092337018e+17
   6.0972196764462810e+17
   6.7384698675270400e+17
   7.4471609299199475e+17
   8.2303856819767232e+17
   9.0959829002668813e+17
   1.0052615772688342e+18
   1.1109858602563712e+18
   1.2278292631485970e+18
   1.3569611939940810e+18
   1.4996740485594657e+18
   1.6573961450606881e+18
   1.8317060192517601e+18
   2.0243482229411579e+18
   2.2372507840526853e+18
   2.4725445029769687e+18
   2.7325842783379528e+18
   3.0199726756098365e+18
   3.3375859744670935e+18
   3.6886029555582029e+18
   4.0765367148108068e+18
   4.5052698236765440e+18
   4.9790931872111176e+18
   5.5027489888943135e+18
   6.0814781519961702e+18
   6.7210727924986010e+18
   7.4279341885389363e+18
   8.2091368465530675e+18
   9.0724993053136814e+18
   1.0026662386494198e+19
   1.1081175674916358e+19
   1.2246593094004847e+19
   1.3534578533000223e+19
   1.4958022583082809e+19
   1.6531171550741899e+19
   1.8269770039599454e+19
   2.0191218527695090e+19
   2.2314747517318812e+19
   2.4661610000341512e+19
   2.7255294165301002e+19
   3.0121758475087553e+19
   3.3289691467965432e+19
   3.6790798882106413e+19
   4.0660120977274061e+19
   4.4936383229520880e+19
   4.9662383908768727e+19
   5.4885422418279211e+19
   6.0657772682979369e+19
   6.7037206324472250e+19
   7.4087570858843619e+19
   8.1879428704062800e+19
   9.0490763392378634e+19
   1.0000776005572131e+20
   1.1052566799547061e+20
   1.2214975396947848e+20
   1.3499635573716302e+20
   1.4919404640690720e+20
   1.6488492123894242e+20
   1.8222601978247286e+20
   2.0139089758026668e+20
   2.2257136317086207e+20
   2.4597939777289005e+20
   2.7184927686435983e+20
   3.0043991489038549e+20
   3.3203745656597683e+20
   3.6695814070852361e+20
   4.0555146526217175e+20
   4.4820368519071852e+20
   4.9534167824711496e+20
   5.4743721730949612e+20
   6.0501169204271369e+20
   6.6864132714134700e+20
   7.3896294938012195e+20
   8.1668036119031775e+20
   9.0257138455105503e+20
   9.9749564569309793e+20
   1.1024031785271021e+21
   1.2183439329023096e+21
   1.3464782828575407e+21
   1.4880886400345899e+21
   1.6445922884849697e+21
   1.8175555693250643e+21
   2.0087095572044880e+21
   2.2199673854830119e+21
   2.4534433935122555e+21
   2.7114742876545719e+21
   2.9966425278257163e+21
   3.3118021736216768e+21
   3.6601074487063940e+21
   4.0450443093423623e+21
   4.4704653330125727e+21
   4.9406282763108614e+21
];
%
%  Next, setup the weights.
%
WT0=[
  0.21035620538389819885E-28
 -0.12644693616088940552E-13
  0.46157312567885668321E-13
 -0.27987033742576678494E-13
  0.54657649654108409156E-13
 -0.26529331099287291499E-13
  0.56749134340673213135E-13
 -0.21572768289772080733E-13
  0.58318460867739760925E-13
 -0.15465892848687829700E-13
  0.60573024556529743179E-13
 -0.85025312590830646706E-14
  0.63880180611476449908E-13
 -0.56596576350102877128E-15
  0.68485006047914070374E-13
  0.85728977321682762439E-14
  0.74650681546818133979E-13
  0.19208372932613381433E-13
  0.82693454289757706437E-13
  0.31701165629228998860E-13
  0.93000040396952081623E-13
  0.46490696394179916916E-13
  0.10604419444905640479E-12
  0.64112165895974571186E-13
  0.12240608340017008854E-12
  0.85217767515070225126E-13
  0.14279579404871719178E-12
  0.11060266069684630524E-12
  0.16808202030984049793E-12
  0.14123670281595459178E-12
  0.19932710117763080694E-12
  0.17830320429641190908E-12
  0.23782981967994220694E-12
  0.22324626027650970672E-12
  0.28517768171070337866E-12
  0.27782855757859949941E-12
  0.34331077352335570356E-12
  0.34420197641399489491E-12
  0.41459976123278392377E-12
  0.42499381982249688902E-12
  0.50194116319499513568E-12
  0.52341213107734220974E-12
  0.60887371932475596308E-12
  0.64337432538102752823E-12
  0.73972052807969328662E-12
  0.78966429788659195358E-12
  0.89976265596232089469E-12
  0.96812431295714248033E-12
  0.10954511874715580138E-11
  0.11858893754914550946E-11
  0.13346662261643230004E-11
  0.14516734901185850330E-11
  0.16270332417806469835E-11
  0.17761192965262209981E-11
  0.19843094598657728049E-11
  0.21722251127126358340E-11
  0.24208558013562638359E-11
  0.26558665246211218806E-11
  0.29542133130007581834E-11
  0.32464334551104291636E-11
  0.36058072230543221354E-11
  0.39676082798212998573E-11
  0.44018068787207206107E-11
  0.48483162182209038194E-11
  0.53741760778847699533E-11
  0.59238861421284167256E-11
  0.65619559488550109069E-11
  0.72374683888302430790E-11
  0.80128318647926483616E-11
  0.88417664804021047188E-11
  0.97850472788000974558E-11
  0.10801152249024990298E-10
  0.11949741288775620211E-10
  0.13194249255521529763E-10
  0.14593803746893311587E-10
  0.16117088182600741411E-10
  0.17823362499440759998E-10
  0.19686960839662181492E-10
  0.21768042712348459703E-10
  0.24047127453754399616E-10
  0.26586169224245549843E-10
  0.29372566166660639848E-10
  0.32471120715874496926E-10
  0.35876995485484108464E-10
  0.39659090711123993131E-10
  0.43821451522206353096E-10
  0.48438566886024388833E-10
  0.53524764256840081139E-10
  0.59161909123775483693E-10
  0.65376353273289546079E-10
  0.72259490983916676274E-10
  0.79851856505621554776E-10
  0.88256972132554113981E-10
  0.97532219231111761839E-10
  0.10779639493701439540E-09
  0.11912700941828918848E-09
  0.13166195190543400137E-09
  0.14550289515667290553E-09
  0.16081145810919829686E-09
  0.17771842706736448843E-09
  0.19641479168712890847E-09
  0.21706652163468219836E-09
  0.23990084518390418007E-09
  0.26512635046403231588E-09
  0.29301487204485379059E-09
  0.32382671796405940321E-09
  0.35788852978338560143E-09
  0.39552347102192978760E-09
  0.43712543089936042400E-09
  0.48309404739375797825E-09
  0.53390563500721587301E-09
  0.59005295736900699180E-09
  0.65211327580989620296E-09
  0.72069283339348376395E-09
  0.79649244503723412979E-09
  0.88025670846750671949E-09
  0.97283758951862838429E-09
  0.10751484374562240226E-08
  0.11882260626931209546E-08
  0.13131897062580609792E-08
  0.14513021636655630472E-08
  0.16039339435115999074E-08
  0.17726240632935599375E-08
  0.19590497332198700268E-08
  0.21650875406672578739E-08
  0.23927891159868709863E-08
  0.26444435360147638239E-08
  0.29225595734392348364E-08
  0.32299302912485969586E-08
  0.35696226515762041770E-08
  0.39450454481729238111E-08
  0.43599472612559631445E-08
  0.48184890913636877388E-08
  0.53252519017629532640E-08
  0.58853155833436851597E-08
  0.65042776355473761329E-08
  0.71883404192425699700E-08
  0.79443429030829021263E-08
  0.87798585629591097164E-08
  0.97032425780220928695E-08
  0.10723743227689399320E-07
  0.11851567478400949647E-07
  0.13098009332253010086E-07
  0.14475537424021579590E-07
  0.15997944513720238914E-07
  0.17680461540553181382E-07
  0.19539933354871179973E-07
  0.21594964684505000780E-07
  0.23866128306162058448E-07
  0.26376149610345229076E-07
  0.29150154762698309478E-07
  0.32215902055657983021E-07
  0.35604079260985092872E-07
  0.39348591789544360374E-07
  0.43486920453657868292E-07
  0.48060478694380069026E-07
  0.53115044437493142888E-07
  0.58701201380017479006E-07
  0.64874861635712681668E-07
  0.71697809408859342572E-07
  0.79238334805050121668E-07
  0.87571902294266750180E-07
  0.96781920558355529704E-07
  0.10696056312048640378E-06
  0.11820970459254829996E-06
  0.13064192692376749733E-06
  0.14438165911984879199E-06
  0.15956640998357920017E-06
  0.17634815657222329376E-06
  0.19489485336503760544E-06
  0.21539212473518470874E-06
  0.23804511154682969141E-06
  0.26308053514448537707E-06
  0.29074895589985508713E-06
  0.32132729115584617374E-06
  0.35512157675298091583E-06
  0.39247003960676882480E-06
  0.43374647340783931675E-06
  0.47936398879211592319E-06
  0.52977913903701991577E-06
  0.58549649799821895163E-06
  0.64707370170465826468E-06
  0.71512703747583724585E-06
  0.79033760405820345763E-06
  0.87345813593705873095E-06
  0.96532052953043588492E-06
  0.10668441762993300587E-05
  0.11790451573234470162E-05
  0.13030464193827640161E-05
  0.14400890071820729325E-05
  0.15915444905568059534E-05
  0.17589286853766969218E-05
  0.19439168303884759854E-05
  0.21483603476946079820E-05
  0.23743053781113989137E-05
  0.26240132540945578888E-05
  0.28999831372928098690E-05
  0.32049770257732388059E-05
  0.35420474020978068237E-05
  0.39145677786671023699E-05
  0.43262664657480666650E-05
  0.47812638810077579977E-05
  0.52841137925459449193E-05
  0.58398488901504216978E-05
  0.64540311583956847791E-05
  0.71328075387126483640E-05
  0.78829714540444842174E-05
  0.87120307957926337232E-05
  0.96282830690788715692E-05
  0.10640898433258409191E-04
  0.11760011483484689510E-04
  0.12996822677619209683E-04
  0.14363710437469910021E-04
  0.15874355032820411062E-04
  0.17543875501207770228E-04
  0.19388980961051210830E-04
  0.21428137844875080257E-04
  0.23681754714302638776E-04
  0.26172386518181301401E-04
  0.28924960326686671307E-04
  0.31967024811679323171E-04
  0.35329025959273332330E-04
  0.39044611784553647914E-04
  0.43150969087563300697E-04
  0.47689195631023510063E-04
  0.52704711454198507913E-04
  0.58247713449355698436E-04
  0.64373677739542765810E-04
  0.71143914895687183756E-04
  0.78626183537770216584E-04
  0.86895368472094227166E-04
  0.96034230136811580777E-04
  0.10613423286950590218E-03
  0.11729646026571620215E-03
  0.12963262680752149015E-03
  0.14326619585462989793E-03
  0.15833361521503999387E-03
  0.17498568260644121238E-03
  0.19338905472215898854E-03
  0.21372791490647819607E-03
  0.23620581621746139576E-03
  0.26104771814107751186E-03
  0.28850223750643329973E-03
  0.31884413578705622606E-03
  0.35237706800135917308E-03
  0.38943662007281590462E-03
  0.43039366566958368833E-03
  0.47565807487700027594E-03
  0.52568281303167982161E-03
  0.58096846836193664822E-03
  0.64206825611437718954E-03
  0.70959354469423893946E-03
  0.78421996374912559657E-03
  0.86669414659246905939E-03
  0.95784118347590099622E-03
  0.10585728435434579864E-02
  0.11698966653840010659E-02
  0.12929259749833240008E-02
  0.14288909657393959393E-02
  0.15791508896276298946E-02
  0.17452075486231769744E-02
  0.19287201026563420645E-02
  0.21315214730933669356E-02
  0.23556362776736309780E-02
  0.26033007312297349808E-02
  0.28769842725422829777E-02
  0.31794136300501061286E-02
  0.35135987233342810820E-02
  0.38828616256736492134E-02
  0.42908672540316900729E-02
  0.47416579734868603835E-02
  0.52396893384671194144E-02
  0.57898709851889693795E-02
  0.63976070720003401851E-02
  0.70688437831027441105E-02
  0.78101127966257229834E-02
  0.86285849757206656285E-02
  0.95321125348557418644E-02
  0.10529286971052299882E-01
  0.11629470426845669312E-01
  0.12842853010275019979E-01
  0.14180454004084899408E-01
  0.15654168426849128515E-01
  0.17276700258609969246E-01
  0.19061578784670331344E-01
  0.21022951738450229575E-01
  0.23175536220878610594E-01
  0.25534136813168170632E-01
  0.28113470567056909888E-01
  0.30927161344728969217E-01
  0.33987341082207328524E-01
  0.37302669043183162012E-01
  0.40877565846017260842E-01
  0.44708454582277110112E-01
  0.48782456555642082774E-01
  0.53070463888373879680E-01
  0.57525214796560822372E-01
  0.62068889249144526543E-01
  0.66590986905823254527E-01
  0.70926870977039285782E-01
  0.74857621553299974471E-01
  0.78074664211727651253E-01
  0.80188872113380424422E-01
  0.80676406709186576638E-01
  0.78917673064227769619E-01
  0.74124063016304961304E-01
  0.65458647531413310938E-01
  0.51957717257333460581E-01
  0.32847972741848592559E-01
  0.74970763258312700383E-02
 -0.23866128698945488634E-01
 -0.60174943784761181220E-01
 -0.98178997988506697125E-01
 -0.13281477972726110637E+00
 -0.15546285697725620301E+00
 -0.15639821579874499391E+00
 -0.12430498665290500016E+00
 -0.54868159863436967438E-01
  0.46862558991703072431E-01
  0.15112182958062059246E+00
  0.21193155344105990556E+00
  0.16951341358877961008E+00
  0.13861974203314799889E-01
 -0.18693504518381348634E+00
 -0.24558896069253360883E+00
 -0.53092694018998139172E-01
  0.25199984157985949595E+00
  0.19682248760574280744E+00
 -0.20143336189691599114E+00
 -0.24584508272586030886E+00
  0.34335593140766357267E+00
 -0.47700665106262918336E-01
 -0.20966284735074519618E+00
  0.25090851237942479734E+00
 -0.16764126613351920669E+00
  0.74951997868967626393E-01
 -0.16951351027511458308E-01
 -0.88285242013749955919E-02
  0.16167873561835501006E-01
 -0.15564117559925290737E-01
  0.12513111252163700016E-01
 -0.92993152469780984010E-02
  0.66505259967794931597E-02
 -0.46677074086718837662E-02
  0.32488855455568548675E-02
 -0.22555291085446750252E-02
  0.15668752111009150857E-02
 -0.10911324863680910667E-02
  0.76253830490409491537E-03
 -0.53525405778204085232E-03
  0.37771045672146697511E-03
 -0.26825282115352968703E-03
  0.19202652749967551236E-03
 -0.13882137232089440237E-03
  0.10159989178100049974E-03
 -0.75497250626312589804E-04
  0.57141729200529359656E-04
 -0.44191609676806338589E-04
  0.35017761995785832974E-04
 -0.28484969541399899885E-04
  0.23801048173241150913E-04
 -0.20412749240321660898E-04
  0.17933634897918600856E-04
 -0.16093677807779768610E-04
  0.14703946920177649808E-04
 -0.13632085400878389976E-04
  0.12785408514801159958E-04
 -0.12099103232765920005E-04
  0.11527807463720849231E-04
 -0.11039642499677090445E-04
  0.10612147244787119201E-04
 -0.10229546722558040478E-04
  0.98808019172690629615E-05
 -0.95581362683731644660E-05
  0.92559894628748428139E-05
 -0.89703680398349605195E-05
  0.86984392262176360162E-05
 -0.84381998993781744882E-05
  0.81881850397080960968E-05
 -0.79472759767526677708E-05
  0.77146204892770847038E-05
 -0.74895903346701209313E-05
  0.72717128455631313132E-05
 -0.70605952027200397989E-05
  0.68558909586342899260E-05
 -0.66573078263653524112E-05
  0.64646106359474119065E-05
 -0.62775967052488337714E-05
  0.60960690045641779074E-05
 -0.59198350149128057049E-05
  0.57487213891732378171E-05
 -0.55825764283683409239E-05
  0.54212557290731716876E-05
 -0.52646114832336199488E-05
  0.51124978277574427215E-05
 -0.49647806501954759280E-05
  0.48213365261592663842E-05
 -0.46820435304229967753E-05
  0.45467775878200409083E-05
 -0.44154180000760036035E-05
  0.42878525031129105819E-05
 -0.41639746681098469138E-05
  0.40436784078298272207E-05
 -0.39268575265170754164E-05
  0.38134098513603608910E-05
 -0.37032392398348510301E-05
  0.35962529271829669603E-05
 -0.34923586009061261710E-05
  0.33914651799141980567E-05
 -0.32934854312419398030E-05
  0.31983363564463991755E-05
 -0.31059371103484061739E-05
  0.30162076790696608685E-05
 -0.29290699115622908088E-05
  0.28444489540948809787E-05
 -0.27622729358204111026E-05
  0.26824715684158770123E-05
 -0.26049757091405561401E-05
  0.25297182336640620443E-05
 -0.24566346909850541750E-05
  0.23856627938558707926E-05
 -0.23167415846242317899E-05
  0.22498114158806991699E-05
 -0.21848145526023768863E-05
  0.21216953789870668849E-05
 -0.20603999269499001719E-05
  0.20008754358939798722E-05
 -0.19430704833474811314E-05
  0.18869353395796189952E-05
 -0.18324219538318030139E-05
  0.17794835990689859392E-05
 -0.17280746718615620423E-05
  0.16781508408187099176E-05
 -0.16296692208678659180E-05
  0.15825882749397249515E-05
 -0.15368675777174579453E-05
  0.14924677441823879966E-05
 -0.14493505440316709914E-05
  0.14074789632159439540E-05
 -0.13668170916979060623E-05
  0.13273299804037580486E-05
 -0.12889836283545319635E-05
  0.12517450527435769451E-05
 -0.12155822896230419416E-05
  0.11804642977313919266E-05
 -0.11463608782321180614E-05
  0.11132426817551759487E-05
 -0.10810812413606430626E-05
  0.10498489463079790445E-05
 -0.10195189699197359345E-05
  0.99006522696480728454E-06
 -0.96146238238294861780E-06
  0.93368585903962674737E-06
 -0.90671180391877650539E-06
  0.88051703794158972103E-06
 -0.85507903355065637431E-06
  0.83037591856657379114E-06
 -0.80638646950736030707E-06
  0.78309007959684342440E-06
 -0.76046672516600126341E-06
  0.73849695280747625374E-06
 -0.71716187755802708277E-06
  0.69644316956135183135E-06
 -0.67632302736954520768E-06
  0.65678415552512063648E-06
 -0.63780975556830963862E-06
  0.61938351997573715339E-06
 -0.60148961697207375561E-06
  0.58411266955661727328E-06
 -0.56723774005096348097E-06
  0.55085032232163946943E-06
 -0.53493633333997110944E-06
  0.51948209868244907623E-06
 -0.50447433644131324069E-06
  0.48990014588904722006E-06
 -0.47574699997612418437E-06
  0.46200273605516768571E-06
 -0.44865554313913170035E-06
  0.43569394953682061395E-06
 -0.42310681390672771474E-06
  0.41088331789010578938E-06
 -0.39901295702377770089E-06
  0.38748552999527151342E-06
 -0.37629112895780040410E-06
  0.36542013202237861199E-06
 -0.35486319617658597678E-06
  0.34461124894551899091E-06
 -0.33465547948734137742E-06
  0.32498733079349469753E-06
 -0.31559849314212400870E-06
  0.30648089749286697872E-06
 -0.29762670812339488137E-06
  0.28902831526969758359E-06
 -0.28067832866533619338E-06
  0.27256957173628210936E-06
 -0.26469507560502209121E-06
  0.25704807272840818962E-06
 -0.24962199077868072385E-06
  0.24241044716871418470E-06
 -0.23540724388966429376E-06
  0.22860636218196940360E-06
 -0.22200195709699679458E-06
  0.21558835236008589946E-06
 -0.20936003566018609291E-06
  0.20331165407816460064E-06
 -0.19743800942055008854E-06
  0.19173405358761880565E-06
 -0.18619488421725469657E-06
  0.18081574059851621267E-06
 -0.17559199964945570756E-06
  0.17051917187001359107E-06
 -0.16559289739464161294E-06
  0.16080894226704971095E-06
 -0.15616319488355259655E-06
  0.15165166247713490835E-06
 -0.14727046762624959541E-06
  0.14301584488191908927E-06
 -0.13888413756276189610E-06
  0.13487179465891340034E-06
 -0.13097536777476139733E-06
  0.12719150812468619494E-06
 -0.12351696364172060964E-06
  0.11994857620999400758E-06
 -0.11648327897320799497E-06
  0.11311809368621049480E-06
 -0.10985012813125349456E-06
  0.10667657363213940060E-06
 -0.10359470265995370341E-06
  0.10060186649751670646E-06
 -0.97695492950586510549E-07
  0.94873084124749734212E-07
 -0.92132214283314405060E-07
  0.89470527774724840543E-07
 -0.86885737009431156146E-07
  0.84375620484431720137E-07
 -0.81938020868743667602E-07
  0.79570843154679285033E-07
 -0.77272052863928245642E-07
  0.75039674297478696595E-07
 -0.72871788831488877790E-07
  0.70766533266929795021E-07
 -0.68722098232649864085E-07
  0.66736726633350816695E-07
 -0.64808712137167796518E-07
  0.62936397705674278282E-07
 -0.61118174170089166613E-07
  0.59352478851276259248E-07
 -0.57637794217708362227E-07
  0.55972646579199561646E-07
 -0.54355604818630417116E-07
  0.52785279162883091201E-07
 -0.51260319990189909104E-07
  0.49779416670227040800E-07
 -0.48341296436233499673E-07
  0.46944723290471018939E-07
 -0.45588496942848318593E-07
  0.44271451780258727462E-07
 -0.42992455864451837250E-07
  0.41750409958276353087E-07
 -0.40544246580819631899E-07
  0.39372929090650920317E-07
 -0.38235450795273419860E-07
  0.37130834085523768675E-07
 -0.36058129594866006882E-07
  0.35016415383563508384E-07
 -0.34004796146768310693E-07
  0.33022402445156628786E-07
 -0.32068389957362020693E-07
  0.31141938754091167181E-07
 -0.30242252593599009938E-07
  0.29368558237615870861E-07
 -0.28520104786762129308E-07
  0.27696163034962460481E-07
 -0.26896024842649531045E-07
  0.26119002528299780677E-07
 -0.25364428277531399452E-07
  0.24631653569082881086E-07
 -0.23920048617305701087E-07
  0.23229001830889051246E-07
 -0.22557919287333441164E-07
  0.21906224222551861542E-07
 -0.21273356535101449785E-07
  0.20658772304732728872E-07
 -0.20061943324939738367E-07
  0.19482356649058811450E-07
 -0.18919514149424371449E-07
  0.18372932089200169455E-07
 -0.17842140706600119178E-07
  0.17326683811178871537E-07
 -0.16826118391794460252E-07
  0.16340014235853439520E-07
 -0.15867953559529469194E-07
  0.15409530648689289756E-07
 -0.14964351510224189348E-07
  0.14532033533448680711E-07
 -0.14112205161252269515E-07
  0.13704505570744190592E-07
 -0.13308584363144380637E-07
  0.12924101262648650824E-07
 -0.12550725823983400048E-07
  0.12188137148392460070E-07
 -0.11836023607829919776E-07
  0.11494082577135149689E-07
 -0.11162020173950579471E-07
  0.10839551006143939185E-07
 -0.10526397926519000691E-07
  0.10222291794616780561E-07
 -0.99269712454067339727E-08
  0.96401824646609211080E-08
 -0.93616789708108062008E-08
  0.90912214031028235422E-08
 -0.88285773158782117222E-08
  0.85735209788006156203E-08
 -0.83258331828536520172E-08
  0.80853010519388850005E-08
 -0.78517178599159764754E-08
  0.76248828529317467947E-08
 -0.74046010768841071059E-08
  0.71906832098687622728E-08
 -0.69829453994641022983E-08
  0.67812091047173317459E-08
 -0.65853009426977313290E-08
  0.63950525394835678792E-08
 -0.62103003854524726628E-08
  0.60308856947513333758E-08
 -0.58566542688268000505E-08
  0.56874563638996400563E-08
 -0.55231465622676621705E-08
  0.53635836473256503502E-08
 -0.52086304821955339059E-08
  0.50581538918636123925E-08
 -0.49120245487234261852E-08
  0.47701168614250102390E-08
 -0.46323088669346301918E-08
  0.44984821257128528939E-08
 -0.43685216199214189850E-08
  0.42423156545711511770E-08
 -0.41197557615253826557E-08
  0.40007366062763786183E-08
 -0.38851558974150727748E-08
  0.37729142987165276171E-08
 -0.36639153437652961865E-08
  0.35580653530469988486E-08
 -0.34552733534349709264E-08
  0.33554510000030559535E-08
 -0.32585125000973940249E-08
  0.31643745396017359674E-08
 -0.30729562113327868796E-08
  0.29841789455041530200E-08
 -0.28979664421992758989E-08
  0.28142446057953132032E-08
 -0.27329414812814390995E-08
  0.26539871924168161407E-08
 -0.25773138816751668664E-08
  0.25028556519244291169E-08
 -0.24305485097912959499E-08
  0.23603303106619080563E-08
 -0.22921407052714058628E-08
  0.22259210878365399183E-08
 -0.21616145456867578372E-08
  0.20991658103504500059E-08
 -0.20385212100542671023E-08
  0.19796286235947108123E-08
 -0.19224374355423869606E-08
  0.18668984927404269626E-08
 -0.18129640620596530537E-08
  0.17605877893741759270E-08
 -0.17097246597221680675E-08
  0.16603309586176070024E-08
 -0.16123642344797300183E-08
  0.15657832621478869464E-08
 -0.15205480074504270167E-08
  0.14766195927971890460E-08
 -0.14339602637660350412E-08
  0.13925333566546899287E-08
 -0.13523032669700050132E-08
  0.13132354188275479691E-08
 -0.12752962352352419048E-08
  0.12384531092355140920E-08
 -0.12026743758811379070E-08
  0.11679292850206920657E-08
 -0.11341879748702300079E-08
  0.11014214463484779337E-08
 -0.10696015381534739059E-08
  0.10387009025592480307E-08
 -0.10086929819117270958E-08
  0.97955198580366660312E-09
 -0.95125286890900208679E-09
  0.92377130945756002788E-09
 -0.89708368833163302057E-09
  0.87116706876644941479E-09
 -0.84599917663709281674E-09
  0.82155838131493233922E-09
 -0.79782367707710818756E-09
  0.77477466505309300292E-09
 -0.75239153569281363935E-09
  0.73065505174126690031E-09
 -0.70954653170499806391E-09
  0.68904783379622916597E-09
 -0.66914134034083904798E-09
  0.64980994263679380996E-09
 -0.63103702625001653087E-09
  0.61280645673505631634E-09
 -0.59510256576828648337E-09
  0.57791013768171107199E-09
 -0.56121439638580773724E-09
  0.54500099267016873900E-09
 -0.52925599187102496775E-09
  0.51396586189505248004E-09
 -0.49911746158916965897E-09
  0.48469802944632955612E-09
 -0.47069517263760014005E-09
  0.45709685636110621144E-09
 -0.44389139349867930081E-09
  0.43106743457132480471E-09
 -0.41861395798487421376E-09
  0.40652026055743952578E-09
 -0.39477594832052842209E-09
  0.38337092758591272230E-09
 -0.37229539627057491552E-09
  0.36153983547227521766E-09
 -0.35109500128849910154E-09
  0.34095191687175519511E-09
 -0.33110186471439289704E-09
  0.32153637915631239045E-09
 -0.31224723910912398970E-09
  0.30322646099050710905E-09
 -0.29446629186269331473E-09
  0.28595920276917881651E-09
 -0.27769788226393698380E-09
  0.26967523012757238777E-09
 -0.26188435126501301881E-09
  0.25431854977949852574E-09
 -0.24697132321776888896E-09
  0.23983635698150977934E-09
 -0.23290751890024982077E-09
  0.22617885396104611209E-09
 -0.21964457919042880688E-09
  0.21329907868420528750E-09
 -0.20713689878085250208E-09
  0.20115274337434959303E-09
 -0.19534146936242159446E-09
  0.18969808222628370222E-09
 -0.18421773173808480806E-09
  0.17889570779236339873E-09
 -0.17372743635793169822E-09
  0.16870847554670939670E-09
 -0.16383451179612808855E-09
  0.15910135616182540748E-09
 -0.15450494071744330869E-09
  0.15004131505843449035E-09
 -0.14570664290687508993E-09
  0.14149719881436189859E-09
 -0.13740936496016368823E-09
  0.13343962804187350973E-09
 -0.12958457625588848896E-09
  0.12584089636512460714E-09
 -0.12220537085144369393E-09
  0.11867487515034671199E-09
 -0.11524637496555659739E-09
  0.11191692366118209830E-09
 -0.10868365972922199935E-09
  0.10554380433023220449E-09
 -0.10249465890504320142E-09
  0.99533602855475026004E-10
 -0.96658091292055187607E-10
  0.93865652846806834872E-10
 -0.91153887549225778384E-10
  0.88520464763624569567E-10
 -0.85963121186075071965E-10
  0.83479658899239781350E-10
 -0.81067943483446593655E-10
  0.78725902182445575131E-10
 -0.76451522122414838220E-10
  0.74242848583018786633E-10
 -0.72097983319807716229E-10
  0.70015082938317186447E-10
 -0.67992357322650702613E-10
  0.66028068126905105916E-10
 -0.64120527350693063962E-10
  0.62268096049913068230E-10
 -0.60469183303296924744E-10
  0.58722245716347976830E-10
 -0.57025788118346391824E-10
  0.55378366976800887211E-10
 -0.53778600071222512844E-10
  0.52225190653632309517E-10
 -0.50716985205207501023E-10
  0.49253109171572777739E-10
 -0.47833283755322262690E-10
  0.46458563164427748566E-10
 -0.45133048314954522821E-10
  0.43867868281155817508E-10
 -0.42690428488605550129E-10
  0.41665890740660430845E-10
 -0.40947061319725302450E-10
  0.40890256062850869212E-10
 -0.42324395200868817252E-10
  0.47175970286618397865E-10
 -0.59920514722972034865E-10
  0.90953607290146280299E-10
];
%
% Get the points at which to evaluate F.
%
NN = size(B,1)*size(B,2);
BB = zeros(NN,1);
BB(1:NN) = B;
Y=YBASE * (1 ./BB');
%
% Evaluate F at these points.
%
YF=feval(funname,Y,varargin{:}).*Y;
%
% Compute the dot product.
%
z0=WT0.'*YF;
%
%  Finally correct for the factor of 1/B.
%
z0=z0' ./ BB;
 
                                                                                                                                                                                                                                                                                                                                     nadwatloop.m                                                                                        0000664 0002715 0000074 00000000433 10372710633 012064  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res=natwat(r,covaij,rij,h)
%  res=natwat(r,covaij,rij,h)
  nr = numel(r);
  npairs = numel(covaij);
  res = r; res(:)=0;
  
  for ii = 1:nr
    res(ii) = sum( covaij .* exp(-(rij-r(ii)).^2 / 2 / h^2 ) ) / ...
              sum( exp(-(rij-r(ii)).^2 / 2 / h^2 ) );
  end
  
                                                                                                                                                                                                                                       nadwat.m                                                                                            0000664 0002715 0000074 00000000442 10274761766 011210  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res=natwat(r,covaij,rij,h)
%  res=natwat(r,covaij,rij,h)
  nr = numel(r);
  npairs = numel(covaij);
  covaij = repmat(covaij,1,nr);
  rij = repmat(rij,1,nr);
  r=reshape(r,1,nr);
  r = repmat(r,npairs,1);
  res=sum(covaij.*exp(-(rij-r).^2/2/h^2)) ./ sum(exp(-(rij-r).^2/2/h^2));
                                                                                                                                                                                                                                newbhankelx.m                                                                                       0000644 0002715 0000074 00000005460 10252061372 012221  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function covacnots = bhankelx(cnots,gm,bcoef)
% res=bhankel(r,gm,knots,bcoef)
% calc hankel tr of sum of bsplines fi*Bi
% and algebraic tail
% 2 pi is not included since I'm going to normalize by cova0
% knots are the breaks in deBoor's notation
  
%global Za dista M cnots indupper rr nt knots Hpol t
global knots Hpol t
  
%% some preliminary definitions and calculations
l = length(knots)-1;
N = length(cnots);
wt = knots(l+1); %% or l*h, h spacing between knots

sp = spmak(t,bcoef);
%continuity of bspline and 1/w^gm
ft = spval(sp,wt);
%ft = bcoef(l+1)/6 + bcoef(l+2)*2/3 + bcoef(l+3)/6;
consta = ft * wt^gm;

%% coeff of polynomials from bspline coeffs
A = coefmat(t,bcoef); % coeff of piecewise polynomials
Aext = repmat(A,[1,1,N+1]);

covacnots = sum(sum(Hpol.*Aext,2),1); %% calc hankel tr. of piecewise pol
cova0 = covacnots(end); 
covacnots = covacnots(1:(end-1));
covacnots = reshape(covacnots,1,N);

%% add tail integrals
cova0 = cova0 - consta*wt^(2-gm)/(2-gm); % add the tail

%% calc tailmat(nnots,ngm); ngm = 100
T=50; %% if r*wt>T use asymptotic approx of tail integral
%if gm>14;%  T=40; %end   %% %%% gm is < 10
ngm=10;
tailmat=zeros(N,ngm);
gmvect = 2.0001 + linspace(.1,10,ngm);
rwt = cnots*wt;

indi = find(rwt>=T);
nindi = length(indi);
for gii = 1:length(gmvect) 
  if(nindi>0)
    for rwii=1:nindi
      rowii = indi(rwii);
      tailmat(indi(rowii),gii) = tailhatapp(rwt(rowii),gmvect(gii));
    end
  end
end

indi = find(rwt<T);
nindi = length(indi);
for gii = 1:length(gmvect) 
  if(nindi>0)
    for rwii=1:nindi
      rowii = indi(rwii);
      tailmat(rowii,gii) = tailhat(rwt(rowii),gmvect(gii));
    end
  end
end



for ii=1:N
  %calc tailint(cnots(ii)*wt,gm), row ii of tailmat
  pp=spline(gmvect,tailmat(ii,:));
  covacnots(ii) = covacnots(ii) + ppval(pp,gm)*cnots(ii)^(gm-2)*consta;
end

%%======================================
covacnots = covacnots/cova0; %% normalize so cova(0)=1

%%======================================
function res=tailhatapp(st,gm)
%% res = tailhat(st,gm)
%% calc int(s^(1-gm) J0(s),{si,st,inf})
    pi4=pi/4;sq2pi=sqrt(2*pi);
    res= ( (-3 + 8*gm)*st^(-0.5 - gm)*cos(pi/4. - st)) / (4.*sq2pi) - ...
         ( (-15 + 16*gm + 128*gm^2)*st^(-1.5 - gm)*cos(pi/4. + st))/(64.*sq2pi) + ...
           sqrt(2/pi)*st^(0.5 - gm)*cos(pi/4. + st) ;
%%======================================
function res=tailhat(st,gm)
% calls Ken Wilder's implementation of 1F2 gh1f2.mexglx
% res = tailhat(st,gm)
% calc int(s^(1-gm) J0(s),{s,st,inf})
a = 1-gm/2;b=2-gm/2;z=-st^2/4;
res = - gamma(-gm/2)*gm/2^gm/gamma(gm/2) + ...
    st^(2-gm)* hg1f2(a,1,b,z) /(gm-2);
%%======================================
function res=coefmat(t,bcoef)
% calc ppform coef of from bspline form
% t = [-3 -2 -1 0:l l+1 l+2 l+3]
% 
ll = length(t) - 7;
res = mexbsplpp(t,bcoef);
res = res./repmat(cumprod([1 1:3]'),1,ll);
res = res(4:-1:1,:)';
                                                                                                                                                                                                                newener.m                                                                                           0000644 0002715 0000074 00000003117 10252343213 011350  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=ener(theta)
%% y=ener(theta,data,dista,M)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,dfstart,bcoef] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt knots Hpol t

  %------------------------------------
  %% values of param
  %------------------------------------
  smoo = theta(1);
%% FIX this with approx for integer smoothness
%  if(smoo-floor(smoo)<0.01)
%    smoo = floor(smoo)+.01;
%  end
%  if(ceil(smoo)-smoo<0.01)
%    smoo = ceil(smoo)-.01;
%  end
   
  gm = 2*(smoo+1);
%  irango = 2*sqrt(smoo)/rango;
  sig2 = theta(2);
  bcoef = theta(3:end);
  ll = length(bcoef) - 3;

  covacnots = bhankelx(cnots,smoo,bcoef); % bhankel forces covacnots(0)=1

  % define matrix cova
  cova = zeros(size(dista));

  % get spline coeff
  cc = spline(cnots,covacnots);
  % interpolate values of cova(rr) using spline coeff
  cova(indupper) = ppval(cc,rr);

  % complete diag and lower triangular part of cova matrix
  cova = cova + cova' + eye(size(cova));
  cova = cova * sig2;
  
  %%DEBUG
%  disp('cond number');disp(cond(cova))
%  plot(cnots,covacnots)

  
if( min(eig(cova)) > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  cc = chol(cova);
  invcc = inv(cc);
  xx = invcc' * Za;
  nobs = size(Za,1);
  xx=reshape(xx,nt*nobs,1);
  y = nt*sum(log(diag(cc))) + 0.5 * xx'*xx;
  %%==========================================
else
  y = 1e300;
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                 nodos.m                                                                                             0000644 0002715 0000074 00000000113 10253637242 011031  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = nodos(wt,ll)
  res = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll * wt;                                                                                                                                                                                                                                                                                                                                                                                                                                                     nogibbs2.m                                                                                          0000664 0002715 0000074 00000032031 10345330221 011407  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% finds mple spectral density values at knots 
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





                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       nogibbs.m                                                                                           0000664 0002715 0000074 00000031677 10321270553 011351  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------

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
  xxx = linspace(0,2*wtf,100);
  %fcovario = myhankel('kcovario',xxx(2:end));
  fcovario = xxx(2:end); 
  for ii = 2:length(xxx)
    fcovario(ii-1) = myhankel('kcovario',xxx(ii));
  end
  
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

mlespect = spectralfn(xxx,smoof,tf,bcoeff)/cova0;
matspect = fmatern(xxx,matirango,matsmoo);
matcovacnots = matsig2*hmatern(cnots,matrango,matsmoo);
rmsecovak = rmse(rohat,simcovacnots);
rmsecovam = rmse(covacnots,simcovacnots);
rmsecovamat = rmse(matcovacnots,simcovacnots);
rmsespeck = rmse(fcovario/2/pi,simspectral(2:end));
rmsespecm = rmse(mlespect(2:end),simspectral(2:end));
rmsespecmat = rmse(matspect(2:end),simspectral(2:end));

%% calc likelihood using rohat
if newsim
  kcova = kcovario(dista);
  if min(eig(kcova))>0
    cc = chol(kcova);
    invcc = inv(cc);
    xx = invcc' * Za;
    xx=reshape(xx,nt*nobs,1);
    klike = -( nt*sum(log(diag(cc))) + 0.5 * xx'*xx );
  else
    disp('kcovario is not positive definite')
    klike=nan;
    save nopositive.mat
  end
end


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
sig2k = kcovario(0);
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
plot(xxx(2:end),sign(fcovario).*sqrt(abs(fcovario/2/pi)),'kx','markersize',3);
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
kker = kcovario(dist2obs);
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
  bsim=zeros(1,100);bmle =zeros(1,100);bker=zeros(1,100);bmat=zeros(1,100);
  invWsim=0;invWmle=0;invWker=0;invWmat=0;
end

lambsim = ksim' * invKsim + bsim'*invWsim*M'*invKsim; lambsim=lambsim';
lambmle = kmle' * invKmle + bmle'*invWmle*M'*invKmle; lambmle=lambmle';
lambker = kker' * invKker + bker'*invWker*M'*invKker; lambker=lambker';
lambmat = kmat' * invKmat + bmat'*invWmat*M'*invKmat; lambmat=lambmat';

%% detect negative fcovario
if newsim
  if (length(find(fcovario<0))>0 | isnan(klike))
    disp(['fker<0 or klike=' num2str(klike)])
    lambker=nan;bker=nan;invWker=nan;kker=nan;
    %  klike = nan;
  end
end

krigZmle = (lambmle' - lambsim')* Za;
krigZker = (lambker' - lambsim')* Za;
krigZmat = (lambmat' - lambsim')* Za;

estvarsim = simsig2 - diag(lambsim'*ksim) + diag(bsim'*invWsim*bsim);
estvarmle = sig2f - diag(lambmle'*kmle) + diag(bmle'*invWmle*bmle);
estvarker = sig2k - diag(lambker'*kker)+ diag(bker'*invWker*bker);
estvarmat = matsig2 - diag(lambmat'*kmat) + diag(bmat'*invWmat*bmat);

msekrigmle = mean(krigZmle.^2,2);
msekrigker = mean(krigZker.^2,2);
msekrigmat = mean(krigZmat.^2,2);

%% clear some variables and save
if newsim
  clear('xx','cc');
  clear('Hpol','tailmat')
  save(savefile);
end





                                                                 nogibbsrain.m                                                                                       0000664 0002715 0000074 00000010202 10373140412 012175  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% finds mple spectral density values at knots 
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






                                                                                                                                                                                                                                                                                                                                                                                              oldbhankelx.m                                                                                       0000644 0002715 0000074 00000005114 10252061551 012201  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function covacnots = bhankelx(cnots,gm,bcoef)
% res=bhankel(r,gm,knots,bcoef)
% calc hankel tr of sum of bsplines fi*Bi
% and algebraic tail
% 2 pi is not included since I'm going to normalize by cova0
% knots are the breaks in deBoor's notation
  
%global Za dista M cnots indupper rr nt knots Hpol t
global knots Hpol t
  
%% some preliminary definitions and calculations
l = length(knots)-1;
N = length(cnots);

%% CHANGE to ft = mexbsplvb()
ft = bcoef(l+1)/6 + bcoef(l+2)*2/3 + bcoef(l+3)/6;
%continuity of bspline and 1/w^gm

wt = knots(l+1); %% or l*h, h spacing between knots

%% coeff of polynomials from bspline coeffs
A = coefmat(t,bcoef); % coeff of piecewise polynomials
Aext = repmat(A,[1,1,N+1]);

covacnots = sum(sum(Hpol.*Aext,2),1);

cova0 = covacnots(end); 
covacnots = covacnots(1:(end-1));
covacnots = reshape(covacnots,1,N);

%% TAIL 
T=50; %% if r*wt>T use asymptotic approx of tail integral
if gm>14
  T=40;
end
consta = ft * wt^gm;
cova0 = cova0 - consta*wt^(2-gm)/(2-gm); % add the tail
%%======================================
%%% calculate tail with asymp approx for r > 100/wt
indi = find(cnots>=T/wt);
if size(indi,1)>0
    indT=min(indi);
    for ik=indT:N
        inte = tailhatapp(cnots(ik)*wt,gm);
        covacnots(ik)=covacnots(ik) + inte * cnots(ik)^(gm-2)*consta;
    end
end
%%======================================
%% calculate tail exactly
%% use 
indi = find(cnots<T/wt);
if size(indi,1)>0
    indT=max(indi);
    for ik=indT:-1:1
        inte = tailhat(cnots(ik)*wt,gm);
        covacnots(ik) = covacnots(ik) + inte * cnots(ik)^(gm-2)*consta;
    end
end
%%======================================
covacnots = covacnots/cova0; %% normalize so cova(0)=1

%%======================================
function res=tailhatapp(st,gm)
%% res = tailhat(st,gm)
%% calc int(s^(1-gm) J0(s),{si,st,inf})
    pi4=pi/4;sq2pi=sqrt(2*pi);
    res= ( (-3 + 8*gm)*st^(-0.5 - gm)*cos(pi/4. - st)) / (4.*sq2pi) - ...
         ( (-15 + 16*gm + 128*gm^2)*st^(-1.5 - gm)*cos(pi/4. + st))/(64.*sq2pi) + ...
           sqrt(2/pi)*st^(0.5 - gm)*cos(pi/4. + st) ;
%%======================================
function res=tailhat(st,gm)
% calls Ken Wilder's implementation of 1F2 gh1f2.mexglx
% res = tailhat(st,gm)
% calc int(s^(1-gm) J0(s),{s,st,inf})
a = 1-gm/2;b=2-gm/2;z=-st^2/4;
res = - gamma(-gm/2)*gm/2^gm/gamma(gm/2) + ...
    st^(2-gm)* hg1f2(a,1,b,z) /(gm-2);
%%======================================
function res=coefmat(t,bcoef)
% calc ppform coef of from bspline form
% t = [-3 -2 -1 0:l l+1 l+2 l+3]
% 
ll = length(t) - 7;
res = mexbsplpp(t,bcoef);
res = res./repmat(cumprod([1 1:3]'),1,ll);
res = res(4:-1:1,:)';
                                                                                                                                                                                                                                                                                                                                                                                                                                                    oldener.m                                                                                           0000644 0002715 0000074 00000002777 10252130221 011341  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y=ener(theta)
%% y=ener(theta,data,dista,M)
%% - likelihood as fn of spectral density
%% defined by param theta
%% theta = [smoo,sig2,dfstart,bcoef] 
%% fknots is log of spectral density values at knots (for siegman and anderson)
%% for bsplines fknots is in the original scale, not log

global Za dista M cnots indupper rr nt knots Hpol t

  %------------------------------------
  %% values of param
  %------------------------------------
  smoo = theta(1);
%% FIX this with approx for integer smoothness
%  if(smoo-floor(smoo)<0.01)
%    smoo = floor(smoo)+.01;
%  end
%  if(ceil(smoo)-smoo<0.01)
%    smoo = ceil(smoo)-.01;
%  end
   
  gm = 2*(smoo+1);
%  irango = 2*sqrt(smoo)/rango;
  sig2 = theta(2);
  bcoef = theta(3:end);
  ll = length(bcoef) - 3;

  covacnots = bhankelx(cnots,smoo,bcoef); % bhankel forces covacnots(0)=1

  % define matrix cova
  cova = zeros(size(dista));

  % get spline coeff
  cc = spline(cnots,covacnots);
  % interpolate values of cova(rr) using spline coeff
  cova(indupper) = ppval(cc,rr);

  % complete diag and lower triangular part of cova matrix
  cova = cova + cova' + eye(size(cova));
  cova = cova * sig2;
  
  %%DEBUG
%  disp('cond number');disp(cond(cova))
%  plot(cnots,covacnots)

  
if( min(eig(cova)) > 0 )
  %------------------------------------
  %% calculate likelihood or restr. likelihood
  %------------------------------------
  y = 0;
  for ii=1:nt
    y = y - gausslikeli(Za(:,ii),cova,0);
  end
  %%==========================================
else
  y = 1e300;
end
 polmatern.m                                                                                         0000664 0002715 0000074 00000001565 10256335440 011724  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res = polmatern(rvect,u,v,a,nu)
%% polmatern(r,u,v,a,nu)
%% cov correp to ((w+u)^2+v^2) ((w-u)^2+v^2) / (a^2+w^2)^(nu+1+2)

p0=(2*a^4 + 2*a^2*nu*(-u^2 + v^2) + nu*(1 + nu)*(u^2 + v^2)^2) / ...
    (2.*a^(2*(2 + nu))*nu*(1 + nu)*(2 + nu));

indnonz = find(rvect>0);
indzeros = find(rvect==0);
r = rvect(indnonz);
res = rvect;

res(indnonz) = ...
besselk(nu,a*r) .*  ...
 (  (r/a).^nu / (2^nu*gamma(1 + nu)) + ...
    (2^(-2 - nu)*(r/a).^(2 + nu)*(a^4 + 2*a^2*u^2 + u^4 - 2*a^2*v^2 + ...
                                 2*u^2*v^2 + v^4))/ gamma(3 + nu)  ...  
  ) + ...
besselk(1 + nu,a*r) .* ...
    ( ((r/a).^(1 + nu)*(-a^2 - u^2 + v^2)) / (2^nu*gamma(2 + nu)) + ...
      ( 2^(-1 - nu)*(1 + nu)*(r/a).^(2 + nu)* ...
       (a^4 + 2*a^2*u^2 + u^4 - 2*a^2*v^2 + 2*u^2*v^2 + v^4)  )./ ...
       (a*r*gamma(3 + nu)) ...
    );

res(indnonz) = res(indnonz)/p0;
res(indzeros)=1;
                                                                                                                                           rainfit.m                                                                                           0000664 0002715 0000074 00000004236 10373211076 011353  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   clear all

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



                                                                                                                                                                                                                                                                                                                                                                  read100sim_2.m                                                                                      0000664 0002715 0000074 00000013073 10314050656 012004  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   clear all

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

                                                                                                                                                                                                                                                                                                                                                                                                                                                                     read100sim.m                                                                                        0000664 0002715 0000074 00000007732 10311340732 011562  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   clear all

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

                                      readall2.m                                                                                          0000664 0002715 0000074 00000012551 10366503007 011404  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function readall()
%% of process and corresponding estimates
addpath('./mfunctions')

%dire = 'remlnu3_ok_100sim/'; %% ok ordinary kriging = use est mean
%dire = 'remlnu3_20nt_100sim/';
%dire ='nu0p5_100sim/';
dire='bic/L10000/';
for nmodelo = 1:4
   read100(nmodelo,dire);
end
%%=================
function read100(nmodelo,dire)

if nmodelo == 1 %'splinetail'
   load([dire '100sim_conv_splinetail.mat'])
elseif nmodelo == 2 %'polmatern'
   load([dire '100sim_conv_polmatern.mat'])
elseif nmodelo == 3 %'matern'
   load([dire '100sim_conv_matern.mat'])
elseif nmodelo == 4; %'expon'
   load([dire '100sim_conv_expon.mat'])
end


indi = ~isnan(mlenuvect);


%% disp summary
disp([dire ' - ' modelo])
disp(strcat('# of simu.= ',num2str(sum(indi))))
%disp('-----------')

disp(strcat('nu = ',num2str([mean(mlenuvect(indi)) std(mlenuvect(indi))])))
disp(strcat('matnu = ',num2str([mean(matnuvect(indi)) std(matnuvect(indi))])))
disp('-----------')
disp(strcat('sig2 = ',num2str([mean(mlesig2vect(indi)) std(mlesig2vect(indi))])))
disp(strcat('sig2k = ',num2str([mean(kersig2vect(indi)) std(kersig2vect(indi))])))
disp(strcat('matsig2 = ',num2str([mean(matsig2vect(indi)) std(matsig2vect(indi))])))
disp('-----------')
disp(strcat('wt = ',num2str([mean(mlewtvect(indi)) std(mlewtvect(indi))])))
disp(strcat('matwt = ',num2str([mean(matwtvect(indi)) std(matwtvect(indi))])))
disp('-----------')
llvect = sum(~isnan(mlebcoefvect),2) - 3;
llmin = min(llvect);llmax=max(llvect);
for jj=llmin:llmax
  disp(strcat('ll = ',num2str(jj), '; # = ', num2str(sum(llvect==jj))))
end
disp('-----------')
disp(strcat('rmsecovm = ',num2str([mean(rmsecovavect(indi,1)) std(rmsecovavect(indi,1))])))
disp(strcat('rmsecovk = ',num2str([mean(rmsecovavect(indi,2)) std(rmsecovavect(indi,2))])))
disp(strcat('rmsecovmat = ',num2str([mean(rmsecovavect(indi,3)) std(rmsecovavect(indi,3))])))
disp('-----------')
disp(strcat('rmsespm = ', num2str([mean(rmsespectvect(indi,1)) std(rmsespectvect(indi,1))])))
disp(strcat('rmsespk = ', num2str([mean(rmsespectvect(indi,2)) std(rmsespectvect(indi,2))])))
disp(strcat('rmsespmat =',num2str([mean(rmsespectvect(indi,3)) std(rmsespectvect(indi,1))])))
disp('-----------')
etrue=likelivect(indi,1);
emle = likelivect(indi,2)-likelivect(indi,2);
eker = likelivect(indi,3)-likelivect(indi,2);
emat = likelivect(indi,4)-likelivect(indi,2);
disp(strcat('liketrue = ',num2str([-mean(etrue) std(etrue)])))
disp(strcat('likemle-mle = ',num2str([mean(emle) std(emle)])))
disp(strcat('likeker-mle = ',num2str([mean(eker) std(eker)])))
disp(strcat('likemat-mle = ',num2str([mean(emat) std(emat)])))
disp('-----------')
msemle = mlemspevect(indi,indi); msemle = mean(msemle);
mseker = kermspevect(indi,indi); mseker = mean(mseker);
msemat = matmspevect(indi,indi); msemat = mean(msemat);

simkrigvar = simkrigvarvect(1,indi);
mlekrigvar = mlekrigvarvect(indi,indi); mlekrigvar = mean(mlekrigvar);
kerkrigvar = kerkrigvarvect(indi,indi); kerkrigvar = mean(kerkrigvar);
matkrigvar = matkrigvarvect(indi,indi); matkrigvar = mean(matkrigvar);

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
%plot(lon,lat,'.')
%hold on;plot(intpoints(:,1),intpoints(:,2),'r+');hold off
%for ii=1:ninterp
%  text(intpoints(ii,1),intpoints(ii,2),num2str(round(mseker(ii)/msemle(ii))))
%  text(intpoints(ii,1),intpoints(ii,2),num2str(round(1000*mlekrigvar(ii)) ...                                               ))
%%end

E0e0 = simkrigvar;
E0e1 = 1+msemle./E0e0;
E0e2 = 1+mseker./E0e0;
E0e3 = 1+msemat./E0e0;
E1e1 = mlekrigvar./(E0e0+msemle);
E2e2 = kerkrigvar./(E0e0+mseker);
E3e3 = matkrigvar./(E0e0+msemat);

disp('  ')
disp('Pred. error - E0e1^2/E0e0^2 - 1')
disp('  ')
disp(strcat('mle : ',num2str([quantile(E0e1',.5)-1 quantile(E0e1',.75)-quantile(E0e1',.25)])))
disp(strcat('ker : ',num2str([quantile(E0e2',.5)-1 quantile(E0e2',.75)-quantile(E0e2',.25)])))
disp(strcat('mat : ',num2str([quantile(E0e3',.5)-1 quantile(E0e3',.75)-quantile(E0e3',.25)])))
disp('  ')


%disp('Pred. variance error - E1e1^2/E0e1^2')
%disp('  ')
%disp(strcat('mle : ',num2str([ quantile(abs(E1e1'-1), .5) quantile(abs(E1e1' -1), .75)- quantile(abs(E1e1'-1),.25)] )))
%disp(strcat('ker : ',num2str([ quantile(abs(E2e2'-1), .5) quantile(abs(E2e2' -1), .75)- quantile(abs(E2e2'-1),.25)] )))
%disp(strcat('mat : ',num2str([ quantile(abs(E3e3'-1), .5) quantile(abs(E3e3' -1), .75)- quantile(abs(E3e3'-1),.25)] )))
%disp('  ')

disp('Pred. variance error log(E1e1^2/E0e1^2)')
disp('  ')
E1e1 = log(E1e1).^2;
E2e2 = log(E2e2).^2;
E3e3 = log(E3e3).^2;
disp(strcat('mle : ',num2str(sqrt([ quantile(E1e1, .5) quantile(E1e1' , .75)- quantile(E1e1',.25)]) ) ))
disp(strcat('ker : ',num2str(sqrt([ quantile(E2e2, .5) quantile(E2e2' , .75)- quantile(E2e2',.25)]) ) ))
disp(strcat('mat : ',num2str(sqrt([ quantile(E3e3, .5) quantile(E3e3' , .75)- quantile(E3e3',.25)]) ) ))
disp('  ')
disp('-----------')
disp('-----------')
disp('  ')
                                                                                                                                                       readall.m                                                                                           0000664 0002715 0000074 00000014156 10343472625 011333  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function readall()
%% of process and corresponding estimates
addpath('./mfunctions')

%dire = 'remlnu3_ok_100sim/'; %% ok ordinary kriging = use est mean
%dire = 'remlnu3_20nt_100sim/';
%dire ='nu0p5_100sim/';
dire='';
for nmodelo = 3
   read100(nmodelo,dire);
end
%%=================
function read100(nmodelo,dire)

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
%indnan = find(~isnan(estimate1(:,end-8-3))); % Kernel positive definite
indnan = find(~isnan(estkervect(:,1)) & ~isnan(estimate1(:,end-8-3)));
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
disp([dire ' - ' modelo])
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
disp(strcat('likemle-mle = ',num2str([-mean(emle)+mean(emle) std(emle-emle)])))
disp(strcat('likemat-mle = ',num2str([-mean(emat)+mean(emle) std(emat-emle)])))
disp(strcat('likeker-mle = ',num2str([-mean(eker)+mean(emle) std(eker-emle)])))

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

E1e1 = log(estvarmle./(E0e0+msemle)).^2;
E2e2 = log(estvarker./(E0e0+mseker)).^2;
E3e3 = log(estvarmat./(E0e0+msemat)).^2;
disp(sqrt([ quantile(E1e1, .5) quantile(E1e1' , .75)- quantile(E1e1',.25)]) )
disp(sqrt([ quantile(E3e3, .5) quantile(E3e3' , .75)- quantile(E3e3',.25)]) )
disp(sqrt([ quantile(E2e2, .5) quantile(E2e2' , .75)- quantile(E2e2',.25)]) )
                                                                                                                                                                                                                                                                                                                                                                                                                  rmse.m                                                                                              0000644 0002715 0000074 00000000153 10247656221 010662  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function res=rmse(x,y)
  if nargin==1
    y=x;
  end
  n=length(x);
  res=sqrt( sum( (x(:)-y(:)).^2 )/n );
                                                                                                                                                                                                                                                                                                                                                                                                                     save_pathdef.m                                                                                      0000644 0002715 0000074 00000005715 10247673415 012362  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function p = pathdef
%PATHDEF Search path defaults.
%   PATHDEF returns a string that can be used as input to MATLABPATH
%   in order to set the path.

  
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2000/06/01 16:19:21 $

% PATH DEFINED HERE -- Don't change this line.

p = [...
'./mfunctions:',...
'.:',...
matlabroot,'/toolbox/matlab/general:',...
matlabroot,'/toolbox/matlab/ops:',...
matlabroot,'/toolbox/matlab/lang:',...
matlabroot,'/toolbox/matlab/elmat:',...
matlabroot,'/toolbox/matlab/elfun:',...
matlabroot,'/toolbox/matlab/specfun:',...
matlabroot,'/toolbox/matlab/matfun:',...
matlabroot,'/toolbox/matlab/datafun:',...
matlabroot,'/toolbox/matlab/audio:',...
matlabroot,'/toolbox/matlab/polyfun:',...
matlabroot,'/toolbox/matlab/funfun:',...
matlabroot,'/toolbox/matlab/sparfun:',...
matlabroot,'/toolbox/matlab/graph2d:',...
matlabroot,'/toolbox/matlab/graph3d:',...
matlabroot,'/toolbox/matlab/specgraph:',...
matlabroot,'/toolbox/matlab/graphics:',...
matlabroot,'/toolbox/matlab/uitools:',...
matlabroot,'/toolbox/matlab/strfun:',...
matlabroot,'/toolbox/matlab/iofun:',...
matlabroot,'/toolbox/matlab/timefun:',...
matlabroot,'/toolbox/matlab/datatypes:',...
matlabroot,'/toolbox/matlab/verctrl:',...
matlabroot,'/toolbox/matlab/demos:',...
matlabroot,'/toolbox/local:',...
matlabroot,'/toolbox/simulink/simulink:',...
matlabroot,'/toolbox/simulink/blocks:',...
matlabroot,'/toolbox/simulink/simdemos:',...
matlabroot,'/toolbox/simulink/simdemos/aerospace:',...
matlabroot,'/toolbox/simulink/simdemos/automotive:',...
matlabroot,'/toolbox/simulink/simdemos/simfeatures:',...
matlabroot,'/toolbox/simulink/simdemos/simgeneral:',...
matlabroot,'/toolbox/simulink/simdemos/simnew:',...
matlabroot,'/toolbox/simulink/dee:',...
matlabroot,'/toolbox/stateflow/stateflow:',...
matlabroot,'/toolbox/stateflow/sfdemos:',...
matlabroot,'/toolbox/stateflow/coder:',...
matlabroot,'/toolbox/compiler:',...
matlabroot,'/toolbox/finance/finance:',...
matlabroot,'/toolbox/finance/calendar:',...
matlabroot,'/toolbox/finance/findemos:',...
matlabroot,'/toolbox/finance/finsupport:',...
matlabroot,'/toolbox/ident/ident:',...
matlabroot,'/toolbox/ident/idobsolete:',...
matlabroot,'/toolbox/ident/idguis:',...
matlabroot,'/toolbox/ident/idutils:',...
matlabroot,'/toolbox/ident/iddemos:',...
matlabroot,'/toolbox/ident/idhelp:',...
matlabroot,'/toolbox/nnet/nnet:',...
matlabroot,'/toolbox/nnet/nnutils:',...
matlabroot,'/toolbox/nnet/nncontrol:',...
matlabroot,'/toolbox/nnet/nndemos:',...
matlabroot,'/toolbox/nnet/nnobsolete:',...
matlabroot,'/toolbox/optim:',...
matlabroot,'/toolbox/pde:',...
matlabroot,'/toolbox/sb2sl:',...
matlabroot,'/toolbox/signal/signal:',...
matlabroot,'/toolbox/signal/fdatoolgui:',...
matlabroot,'/toolbox/signal/sptoolgui:',...
matlabroot,'/toolbox/signal/sigdemos:',...
matlabroot,'/toolbox/stats:',...
matlabroot,'/toolbox/symbolic:',...
matlabroot,'/toolbox/wavelet/wavelet:',...
matlabroot,'/toolbox/wavelet/wavedemo:',...
     ...
];

p = [userpath,p];
                                                   sieghankel.m                                                                                        0000664 0002715 0000074 00000002346 10344213650 012027  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function [rhon gg] = sieghankel(h_fun,R,varargin)
%% [rhon gg] = sieghankel(@fun,R,varargin)
%%
%% calculates hankel tr of fun
%% g(rho) = int_0^inf r fun(r) J_nu(r rho) dr
%% siegman uses g(rho) = 2pi  int_0^inf r fun(r) J_nu(2pi r rho) dr
%%
%% h_fun is handle of fun
%%
%% R is range of the function. The sampling is exponential on
%% (sqrt(rrho),sqrt(bbeta))*R on spatial domain
%% (sqrt(rrho),sqrt(bbeta))/R on spectral domain
%%
%% varargin has the additional arguments to call h_fun
%%
%% rhon and gg(rhon) is returned
  
%% parameters
  global K1 K2 N
  K1 = 2; %% r0<=1/(K1 b)
  K2 = 2; %% dr_N ~= a b
  N = 2^12; %% number of sample points
  %% b = r0 exp(a N); beta = rho0 exp(a N)
  
  % calc a, bbeta, rrho
  a = fzero(@afun,log(K1/K2));
  bbeta = fzero(@bbetafun,N/K2);
  rrho = a*K2/K1^2;
  
  b = sqrt(bbeta);
  beta = b;
  r0 = sqrt(rrho);
  
  rn = r0*R * exp(a*(0:(N-1)));
  fn = [rn .* h_fun(rn,varargin{:})  zeros(1,N)];
  jm = a*rrho*exp(a*(0:(2*N-1))).*besselj(0,rrho*exp(a*(0:(2*N-1))));

  gm = real(fft(fft(fn).*ifft(jm)));

  rhon = rn/R^2;
  gg = gm(1:N)./rhon;

%% functions
  
function res = afun(x)
  global N K1 K2
  res = x*exp(x*N)-K1/K2;
  
function res = bbetafun(x)
  global N K1 K2
  res = K2*x*log(K1*x) - N;
  
  
                                                                                                                                                                                                                                                                                          simul100_2.m                                                                                        0000664 0002715 0000074 00000013050 10350173604 011503  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   % calls nogibbs to get Nsim simulations
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




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        simul100.m                                                                                          0000664 0002715 0000074 00000013251 10307347526 011275  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% calls nogibbs to get Nsim simulations
%% of process and corresponding estimates
addpath('./mfunctions')

type simul100
type nogibbs

global Za dista M cnots indupper rr nt reml
global Hpol tailmat smoovect wtvect % tabulated values
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

ll = 3 % number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

nsimini=1;Nsim = 200;
%testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
testsplinetail = false; testmatern = true; testpolmatern = false;
nparam = 2+ll+3+1;
fixwt = false;
gibbs = false;
reml = false;

%% simulated annealing parameters
L = 20000; %%max number of iterations
Ti = 1000; %% scale temperature
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
elseif testpolmatern 
  modelo = 'polmatern';
elseif testmatern
  modelo = 'matern';
end

if nsimini ==1
  inivect1   = nan(Nsim,nparam);
  inivect50  = nan(Nsim,nparam);
  inivect100 = nan(Nsim,nparam);
  
  estimate1  = nan(Nsim,nparam+8); %% thetaf
  estimate50 = nan(Nsim,nparam+8); %% thetaf
  estimate100= nan(Nsim,nparam+8); %% thetaf
  
  msemlevect = nan(Nsim,ninterp); %% mse true kriged and mle kriged
  msekervect = nan(Nsim,ninterp); %% mse true kriged and ker kriged
  msematvect = nan(Nsim,ninterp); %% mse true kriged and mat kriged
  
  estsimvect = nan(1,ninterp);   %% estimated variance with sim param
  estmlevect = nan(Nsim,ninterp);%% estimated variance with mle param
  estkervect = nan(Nsim,ninterp);%% estimated variance with ker param
  estmatvect = nan(Nsim,ninterp);%% estimated variance with mat param
  
  save(strcat('100sim_',modelo,'.mat'),'inivect1','inivect50',...
       'inivect100','estimate1','estimate50',...
       'estimate100','msemlevect','msekervect','msematvect',...
       'estsimvect','estmlevect','estkervect','estmatvect',...
       'lon','lat','intpoints')

else
  load(strcat('100sim_',modelo,'.mat'))
end

for nsim=nsimini:Nsim

%if length(find(nsim==[3 19 31]))==0
  testwt = wtvect(1);
  %  rand('state',sum(100*clock)); randn('state',sum(100*clock));
  rand('state',1000*nsim); randn('state',11000*nsim);
  newsim = true;
  pref100sim = strcat('100sim_wt001_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect1(nsim,:)=thetavect(1,:);
  estimate1(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  msemlevect(nsim,:) = msekrigmle;
  msekervect(nsim,:) = msekrigker;
  msematvect(nsim,:) = msekrigmat;
  estsimvect(1,:) = estvarsim;
  estmlevect(nsim,:) = estvarmle;
  estkervect(nsim,:) = estvarker;
  estmatvect(nsim,:) = estvarmat;
  dEmin = dE;
  disp(strcat('---------',num2str(nsim),'--','wtvect(1)---------'))
  
  testwt = wtvect(50);
  newsim = false;
  pref100sim = strcat('100sim_wt050_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect50(nsim,:)=thetavect(1,:);
  estimate50(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  if dE<dEmin
    dEmin = dE;    
    msemlevect(nsim,:) = msekrigmle;
    msekervect(nsim,:) = msekrigker;
    msematvect(nsim,:) = msekrigmat;
    estsimvect(1,:) = estvarsim;
    estmlevect(nsim,:) = estvarmle;
    estkervect(nsim,:) = estvarker;
    estmatvect(nsim,:) = estvarmat;
  end
  disp(strcat('---------',num2str(nsim),'--','wtvect(50)---------'))

  testwt = wtvect(100);
  newsim = false;
  pref100sim = strcat('100sim_wt100_',sprintf('%03d',nsim),'_'); 
  nogibbs
  inivect100(nsim,:)=thetavect(1,:);
  estimate100(nsim,:)=[thetavect(end,:) truene Ef klike ...
                     rmsecovam rmsecovak rmsespecm rmsespeck ...
                    sig2k];
  if dE<dEmin
    dEmin = dE;    
    msemlevect(nsim,:) = msekrigmle;
    msekervect(nsim,:) = msekrigker;
    msematvect(nsim,:) = msekrigmat;
    estsimvect(1,:) = estvarsim;
    estmlevect(nsim,:) = estvarmle;
    estkervect(nsim,:) = estvarker;
    estmatvect(nsim,:) = estvarmat;
  end
  disp(strcat('---------',num2str(nsim),'--','wtvect(100)---------'))

  %% save results for later analysis
  save(strcat('100sim_conv_',modelo,'.mat'),'inivect1','inivect50',...
       'inivect100','estimate1','estimate50',...
       'estimate100','msemlevect','msekervect','msematvect',...
       'estsimvect','estmlevect','estkervect','estmatvect',...
       'lon','lat','intpoints')

%end
end


                                                                                                                                                                                                                                                                                                                                                       simulate2.m                                                                                         0000664 0002715 0000074 00000005423 10343635570 011631  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt % observations
global Hpol tailmat smoovect wtvect % tabulated values
global covario rsort rmax;

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

ll = 6; %number of polynomial pieces
load(strcat('tablahpol_',num2str(ll)))
load tablatailmat %% reads tailmat matrices

%% calculate or read tabulated values
%% hankel transform of polynomials
%tic
%Hpol = calcHpol(ll,cnots,wtvect);
%disp('hPol calc time')
%toc
%save(strcat('tablahpol_',num2str(ll)),'Hpol')

%%% hankel tranform of tail
%tic
%tailmat = calcTailmat(cnots,smoovect,wtvect);
%disp('hPol calc time')
%toc
%save tablatailmat Hpol tailmat


testsplinetail = true; testmatern = false; testpolmatern = false;
%testsplinetail = false; testmatern = false; testpolmatern = true;
%testsplinetail = false; testmatern = true; testpolmatern = false;
nt = 200; %number of replicates in time
Za = zeros(nobs,nt);
%%randn('state',0);
%------------------------------------
% generate observations
%------------------------------------
indi = find(smoovect>=3); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simwt = wtvect(find(wtvect <= 0.010)); simwt = simwt(end);
simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
simsig2=1;
simll = 4;

if testsplinetail
  simt = nodos(simwt,simll);
  simknots = simt(4:end-3);
  %  simbcoef = [0.8 1.0 .8 .2  1.2 1 .8 .6 .4]; %% ll=6 should have ll+3 elements
  simbcoef = [0.8 1 .8  .2  .8 .6 .2]; %% ll=4
  simbcoef(end)=contderiv(simbcoef,simgm);
  simtheta = [simsmoo simsig2 simbcoef simwt];
  load(strcat('tablahpol_',num2str(simll)))
  simcova = calcCova(simtheta);
elseif testpolmatern 
  % test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
  u = 0.5*simirango; v=.1*simirango;
  simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);
elseif testmatern
  % test matern cova
  simcova = simsig2*hmatern(dista,simrango,simsmoo);
end

for ii=1:nt
  Za(:,ii) = simulachol(simcova);
end

if testsplinetail
  save(strcat('splinetail_',num2str(nt),'_',num2str(simll),'.mat'),...
       'Za','simbcoef','simll','simsig2','simsmoo','simgm','simwt',...
       'simt','simknots', 'simtheta','simcova')
elseif testpolmatern 
  save(strcat('polmatern_',num2str(nt),'_','.mat'),...
       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango','u','v',...
       'simcova')
elseif testmatern
  save(strcat('matern_',num2str(nt),'.mat'),...
       'Za','simsig2','simsmoo','simgm','simwt','simrango','simirango',...
       'simcova')
end                                                                                                                                                                                                                                             spectralfn.m                                                                                        0000644 0002715 0000074 00000000710 10261315670 012050  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   function y = spectralfn(w,smoo,t,bcoef)
% arbitrary spectral density function
  y = w;
%  bcoef = [.2 1 .2 .8 .4 .2];
  sp = spmak(t,bcoef);
  wt = t(end-3);
  
  indismall = find(w<=wt);
  if length(indismall)>0
    y(indismall) = spval(sp,w(indismall));
  end
  
  indilarge = find(w>wt);
  if length(indilarge)>0  
    ft = spval(sp,wt);
    gm = smoo*2 + 2;
    consta = ft*wt^(gm);
    y(indilarge) = consta./w(indilarge).^(gm);
  end
  
  y = y/2/pi;                                                        testktilde.m                                                                                        0000664 0002715 0000074 00000004431 10341162510 012061  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% function res = mplebsplines()
%% finds mple spectral density values at knots 
%% and interpolating with cubic spline + algebraic tail
%% ------------------------------------
addpath('./mfunctions')

global Za dista M cnots indupper rr nt
global smoovect wtvect % tabulated values
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

%% smoothness and wt values for tabulation
ngm = 100;
smoovect = linspace(.05,5,ngm)+.0001;
gmvect = smoovect*2 + 2;
nwt = 100;
wtvect = linspace(1/rmax,1/rmin,nwt);

nt = 200; %number of replicates in time
Za = zeros(nobs,nt);

% generate observations

simsmoo = 3;
indi = find(smoovect>=simsmoo); simsmoo = smoovect(indi(1)); %% the first smoo>=3
simgm=2*(simsmoo+1);
simwt = 0.0100; %% simwt=0.0038
simwt = wtvect(find(wtvect <= simwt)); simwt = simwt(end);
simirango = simwt; simrango = 2*sqrt(simsmoo)/simirango; 
simsig2=1;
simll = 4;

% test polynomial matern ((w-u)^2+v^2)((w+u)^2+v^2)/(w^2+a^2)^nu+1+k
u = 0.5*simirango; v=.1*simirango;
simcova = simsig2*polmatern(dista,u,v,simirango,simsmoo);

for ii=1:nt
  Za(:,ii) = simulachol(simcova);
end  

%% calc covariogram cloud and nadwat estimate of cov fn

dindupper = find(triu(ones(size(dista))));

ncov = length(dindupper);
covario = nan(ncov,nt);
[rsort,indsort] = sort(dista(dindupper));

zm = 0;

for ii=1:nt
  z = Za(:,ii); 
  zdif = repmat(z,1,nobs);
  zdif = (zdif -zm).*(zdif -zm)';
  zdiflin = zdif(dindupper);
  covario(:,ii) = zdiflin(indsort);
end
  
covario = sum(covario,2)/nt;

rohat = kcovario(cnots);

tic
%% hankel transform of rohat
wtf = 1/50;
xxx = linspace(0,2*wtf,100);
fcovario = xxx(2:end); 
for ii = 2:length(xxx)
  fcovario(ii-1) = myhankel('kcovario',xxx(ii));
end
toc

R = 500;

tic
[ww,gg] = sieghankel(@kcovario, R);
pp = spline(ww,gg);
fcovario2 = ppval(pp,xxx);
toc

plot(ww,gg);hold on;plot(xxx(2:end),fcovario,'r.','markersize',5);hold off

tic
  R = 300;
%% get pos def covariogram
kk = ktilde(cnots,R);
plot(cnots,kk)
hold on; plot(cnots,kcovario(cnots),'r.');hold off
toc


                                                                                                                                                                                                                                       test.m                                                                                              0000644 0002715 0000074 00000014262 10252413726 010676  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% test bhankelx(cnots,gm,bcoef)~hmatern(cnots,rango,smoo)
clear all
global Za dista M cnots indupper rr nt knots Hpol t tailmat smoovect
ngm = 50;
smoovect = linspace(.05,5,ngm)+.0001;
smoo=2.0001;
gm=2*smoo+2;
rango=50;
irango=2*sqrt(smoo)/rango;
wt=2*irango;
cnots = (.1:.1:4)*rango;
%% hankel transform of polynomials
Hpol = calcHpol(knots,cnots);
%% hankel tranform of tail
tailmat = calcTailmat(cnots*wt);

ll=100;
t = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*wt;
knots = (0:ll)/ll*wt;
bcoef = fmatern(abs(t(3:end-2)),irango,smoo);
sp = spmak(t,bcoef)
xx = linspace(0,wt,20);
hold off
plot(xx,spval(sp,xx),'r.')
hold on
plot(xx,fmatern(xx,irango,smoo))
hold off

Hpol = calcHpol(knots,cnots);
trucova = hmatern(cnots,rango,smoo);
mycova = bhankelx(cnots,smoo,bcoef);

plot(cnots(1:end),trucova(1:end))
hold on
plot(cnots(1:end),mycova(1:end),'r.')
hold off

rmse(trucova,mycova)

[thetaff,fval]=fminunc(ener,thetaf)

hold off
sp = spmak(t,thetaf(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmf./xx.^gmf,'r.')
grid
print -depsc plotf.eps

hold off
sp = spmak(t,thetaff(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmff./xx.^gmff,'r.')
grid
print -depsc plotff.eps

hold off
sp = spmak(t,thetafff(3:end));
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gmfff./xx.^gmfff,'r.')
grid
print -depsc plotfff.eps



kksmoof = 2.4251;
kksig2f = 1.0129;
kkbcoeff=[ 
    1.5273
    8.1425
    1.5273
    0.2079
    0.0135
    0.2214
];

hold off
sp = spmak(t,bcoeff);
fnplt(sp,[0,.01])
% plot the tail
xx=.01:.001:.02;
hold on
plot(xx,spval(sp,.01)*wt^gm./xx.^gm,'r.')




ener([smoof sig2f [bcoeff(1:1) 100 bcoeff(3:end)]])



smoof = 2.5216;
sig2f = 1.1226;
kkbcoeff=[    0.4566
    0.3851
    0.4566
    0.0581
  100.4991
    0.6671
    0.7671
   20.4019
    7.9151
    6.6429
    0.3399
    0.8463
    3.6148
    0.8002
    0.7406
    0.1371
    0.6588
    0.2072];





xx =linspace(1,100,10) ;
yy = linspace(1,100,10);
[xxx,yyy]=meshgrid(xx,yy);
for ii=1:length(xxx);
   for jj=1:length(yyy);
     zlik(ii,jj)=-ener([smoof sig2f [bcoeff(1:3) xxx(ii,jj) yyy(ii,jj) bcoeff(6:end)]]);
   end;
 end
meshc(xxx,yyy,zlik)
xlabel('f(4/3 1e-3'),ylabel('f(2e-3)'),zlabel('likeli')


rango=2000;smoo=6.01;gm=2*(smoo+1);sig2=3;
irango = 2*sqrt(smoo)/rango;
knots = (0:ll)/ll*irango*alpha; %% interval of interest (0:wt)
t= [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*irango*alpha; %% deBoor knots
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); %% coef of bsplines approx
bcoef(end) = contderiv(bcoef,gm);  %% continuity of deriv at wt
cnots = 



%% TEST: bhankel ~= true covariance
rango=500;smoo=2.01;gm=2*(smoo+1);
ll=200;cnots = rango*(.1:.1:2);
irango = 2*sqrt(smoo)/rango;
wt=3*irango;
knots = (0:ll)/ll*wt; %% interval of interest (0:wt)
t= [-3 -2 -1 0:ll ll+1 ll+2 ll+3]/ll*wt; %% deBoor knots
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); %% coef of bsplines approx
bcoef(end) = contderiv(bcoef,gm);  
figure(1)
xx = linspace(0,wt,20);
sp=spmak(t,bcoef); 
hold off
plot(xx,spval(sp,xx),'r.')
hold on
plot(xx,fmatern(xx,irango,smoo))
hold off

Hpol = calcHpol(knots,cnots);

trucova=hmatern(cnots,rango,smoo);
mycova=bhankelx(cnots,gm,bcoef);
disp(rmse(trucova,mycova))
hold off
plot(cnots,trucova,'b-')
hold on
plot(cnots,mycova,'.r')

print -depsc plot.eps



%% TEST: profile likelihood
simrango=2000;simsmoo=6;simsig2=1;simgm=2*(simsmoo+1);
z = simulachol(hmatern(dista,simrango,simsmoo)*simsig2);

xx=.2:.4:6;
pll = zeros(size(xx));
for ii=1:length(xx)
  pll(ii)=enersmoo(xx(ii),simrango,sig2,z,dista)
end
plot(xx,pll,'.')
pll
hold on; plot([0 6],(min(pll)+2)*[1 1]);hold off





indupper = find(triu(dista,1));
rr=dista(indupper);
kkcova = zeros(size(dista));kksimcova=kkcova;

rmin=min(rr); rmax=max(rr);
kkcnots = linspace(rmin,rmax,100);
alpha = 5;

simrango = 1000;
simsmoo = 1.1; simgm = (simsmoo+1)*2;
simirango = 2*sqrt(simsmoo)/simrango;
simknots = (0:ll)*simirango/alpha;
simt = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]*simirango/alpha;
simbcoef = fmatern(abs(simt(3:(end-2))),simirango,simsmoo); 
simbcoef(end)=contderiv(simbcoef,simgm);
kksimcovacnots = bhankel(kkcnots,simgm,simknots,simbcoef);

rango = 1000;
smoo = 5.1; gm = (smoo+1)*2;
irango = 2*sqrt(smoo)/rango;
knots = (0:ll)*irango/alpha; 
t = [-3 -2 -1 0:ll ll+1 ll+2 ll+3]*irango/alpha;
bcoef = fmatern(abs(t(3:(end-2))),irango,smoo); 
bcoef(end)=contderiv(bcoef,gm);
kkcovacnots = bhankel(kkcnots,gm,knots,bcoef);

cc = spline(kkcnots,kksimcovacnots);
kksimcova(indupper)=ppval(cc,rr);
cc = spline(kkcnots,kkcovacnots);
kkcova(indupper)=ppval(cc,rr);

kksimcova = kksimcova + kksimcova' + eye(size(kkcova));
kkcova = kkcova + kkcova' + eye(size(kkcova));

cond(kksimcova)
cond(kkcova)

plot(kkcnots,kksimcovacnots,'b.');hold on
plot(kkcnots,kkcovacnots,'r.');hold off




xx=t(3:end-2);
plot(xx,simbcoef,'b.');hold on
plot(xx,bcoef,'r.');hold off








%%here
rango = 500;
gm=3.1;smoo =gm/2 - 1;
cnots = (.1:.1:5)*rango;
protknots=0:20;
knots = protknots/rango/5;
lend = protknots(end)
t = [-3 -2 -1 protknots lend+1 lend+2 lend+3]/rango/5;
xint = t(3:(end-2));
bcoef = fmatern(xint,1/rango,smoo)*2*pi;
covacnots = bhankel(cnots,gm,knots,bcoef);
notailcovacnots = kkbhankel(cnots,gm,knots,bcoef)
plot(cnots,matern(cnots,rango,smoo));hold on
plot(cnots,covacnots,'r.');
plot(cnots,notailcovacnots,'g.');hold off















tic
for ii=1:NN
  st = cnots(ii)/rango*2;
  covatail(ii)=hg1f2(1-gm/2,1,2-gm/2,-st^2/4)*st^(2-gm)/(gm-2) ...
               - gm*gamma(-gm/2)/gamma(gm/2)/2^gm ;
  covatailapp(ii)=kktail(st,gm);
end
toc

plot(cnots(100:end)*2/rango,covatail(100:end),'.');hold on;grid
plot(cnots(100:end)*2/rango,covatailapp(100:end),'r.');hold off


gm=15.1;
knots=nodos(1/30);
fknots=repmat(1,1,length(cnots));
kkcova=bhankel(cnots,gm,knots,fknots);
plot(cnots*2/rango,kkcova)

%=============================================
l=10;
knots = 0:10;
xx= [-3 -2 -1 knots 11 12 13];
xint = xx(3:(end-2));
n = length(xint)
%yy=(xint-3).*(xint-6).*(xint-9);
%yy = repmat(0,1,n);yy(6)=1;
yy= xint.^3;
kk=mexbsplpp(xx,yy)

ll = n - 3;
kk = kk./repmat(cumprod([1 1:3]'),1,ll)
kk=kk(4:-1:1,:)';
pp=mkpp(xx(4:n+1),kk);


plot(xint,yy,'r');
xtest=xint(1)+(1:100)/100*(xint(end)-xint(1));
hold on;
plot(xtest,ppval(pp,xtest),'b-')
hold off
                                                                                                                                                                                                                                                                                                                                              toyexample.m                                                                                        0000664 0002715 0000074 00000004355 10310127432 012101  0                                                                                                    ustar   haky                            phds                                                                                                                                                                                                                   %% calls nogibbs to get Nsim simulations
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
rmse(c0,c2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   