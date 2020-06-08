function readall()
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
