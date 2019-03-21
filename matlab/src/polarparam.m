
function [fvals] = polarparam(t,polarfun)

fvalspolar = polarfun(t);
r = fvalspolar(:,1); th = fvalspolar(:,2); 
dr = fvalspolar(:,3); dth = fvalspolar(:,4);
d2r = fvalspolar(:,5); d2th = fvalspolar(:,6);

cth = cos(th);
sth = sin(th);

xs = r.*cth;
ys = r.*sth;
dxs = -r.*sth.*dth+dr.*cth;
dys = r.*cth.*dth+dr.*sth;
d2xs = -r.*cth.*dth.^2-r.*sth.*d2th-dr.*sth.*dth-dr.*sth.*dth+d2r.*cth;
d2ys = -r.*sth.*dth.^2+r.*cth.*d2th+dr.*cth.*dth+dr.*cth.*dth+d2r.*sth;

fvals(:,1) = xs;
fvals(:,2) = ys;
fvals(:,3) = dxs;
fvals(:,4) = dys;
fvals(:,5) = d2xs;
fvals(:,6) = d2ys;

end

