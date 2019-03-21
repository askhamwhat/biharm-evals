
function [fvalsp] = crescent(t)
%CRESCENT polar parametric description of a crescent shape
%r(s) = 0.2/(1+exp (4(s−3π/2)(s−π/2))) + 0.4, θ(s) = −49/50 π sin s

et = exp(4*(t-3*pi/2.0).*(t-pi/2.0));
r = 0.2./(1+et) + 0.4;
dr = 8.0*et.*(pi-t)./(5.0*(1+et).^2);
d2r = (8*et.*( -8*t.^2 + et.*(8*t.^2 - 16*pi*t + 8*pi^2 - 1) + 16*pi*t ...
    - 8*pi^2 - 1))./(5*(et + 1).^3);


th = -49.0/50.0*pi*sin(t);
dth = -49.0/50.0*pi*cos(t);
d2th = 49.0/50.0*pi*sin(t);

fvalsp(:,1) = r;
fvalsp(:,2) = th;
fvalsp(:,3) = dr;
fvalsp(:,4) = dth;
fvalsp(:,5) = d2r;
fvalsp(:,6) = d2th;

end

