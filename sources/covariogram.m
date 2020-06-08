ncov = length(indupper);
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
