function y=fmatern(x,irange,smoo)
%% calculates the matern spectral density 1/(a^2+w^2)^(nu+1/2)
    d=2;
    y = gamma(smoo+d/2)/gamma(smoo) * irange^(2*smoo) / pi^(d/2) ...
        ./ (irange^2 + x.^2 ).^(smoo+d/2);
