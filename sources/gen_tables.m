%% generate tables for hankel transform of S+T
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


