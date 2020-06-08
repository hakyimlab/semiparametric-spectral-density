function tailmat = calcTailmat(cnots,smoovect,wtvect)
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
    st^(2-gm)* hg1f2(a,1,b,z) /(gm-2);