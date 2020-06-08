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
