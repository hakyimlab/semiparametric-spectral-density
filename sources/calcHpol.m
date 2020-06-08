function Hpol = calcHpol(ll,cnots,wtvect)
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
  