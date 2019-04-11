
function [fvals] = smoothbox(t,W,H,p)
%SMOOTHBOX discretize a smooth box over [0,2*pi]
%g(t) = (W cos t, H sin t)./ ((cos t).^p + (sin t).^p).^(1/p)

ct = cos(t);
st = sin(t);
ctp = ct.^p;
stp = st.^p;

xs = W*ct./((ctp + stp).^(1.0/p));
ys = H*st./((ctp + stp).^(1.0/p));
dxs = -W*st.^(p-1)./((ctp + stp).^(1.0/p+1.0));
dys = H*ct.^(p-1)./((ctp + stp).^(1.0/p+1.0));
d2xs = st.^(p-2).*(2*ct.^(p-1).*cos(2*t) + st.^(p-1).*sin(2*t) ...
    - p*ct.^(p-1))./(( ctp + stp ).^(1.0/p+2.0));
d2ys = ct.^(p-2).*(ct.^(p-1).*sin(2*t) - st.^(p-1).*cos(2*t) ...
    - p*st.^(p-1))./(( ctp + stp ).^(1.0/p+2.0));

fvals(:,1) = xs;
fvals(:,2) = ys;
fvals(:,3) = dxs;
fvals(:,4) = dys;
fvals(:,5) = d2xs;
fvals(:,6) = d2ys;

end

