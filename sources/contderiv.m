function flp1=contderiv(bcoef,gm)
%% fnp1 = contderiv
%% returns f_(n+1) such that deriv is cont at wt
%% called by bhankel.m and bproposal.m
    l=length(bcoef)-3;
    fl0 = bcoef(end-1);
    flm1 = bcoef(end-2);
    flp1 = ( flm1*(3*l - gm) - fl0*4*gm )/(3 * l + gm);
