function res = fpolmatern(w,u,v,a,nu)
p0=(2*a^4 + 2*a^2*nu*(-u^2 + v^2) + nu*(1 + nu)*(u^2 + v^2)^2) / ...
    (2.*a^(2*(2 + nu))*nu*(1 + nu)*(2 + nu));

res = ( (w-u).^2 + v^2 ).*( (w+u).^2 + v^2) ./ ...
        ( w.^2 + a^2).^(nu+1+2);

res = res/p0/2/pi;