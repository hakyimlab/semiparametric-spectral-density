function covacnots = bhankelx(cnots,gm,bcoef)
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
