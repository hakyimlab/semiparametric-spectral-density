function thetaf=bproposal(theta,fixwt,gibbs)
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
  y = vect(ceil( rand(1)*nvect ) );