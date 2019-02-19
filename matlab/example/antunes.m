
function [fvals] = antunes(t)

ct = cos(t);
st = sin(t);
c2t = cos(2.0*t);
s2t = sin(2.0*t);

xs = 2.0*ct;
ys = st + 0.25*s2t;
dxs = -2.0*st;
dys = ct + 0.5*c2t;
d2xs = -2.0*ct;
d2ys = -st - s2t;

fvals(:,1) = xs;
fvals(:,2) = ys;
fvals(:,3) = dxs;
fvals(:,4) = dys;
fvals(:,5) = d2xs;
fvals(:,6) = d2ys;

end

