
function [fvals] = circle(t,varargin)

r = 1.0;
x0 = 0.0;
y0 = 0.0;
if nargin > 1
    r = varargin{1};
end
if nargin > 2
    ctr = varargin{2};
    x0 = ctr(1); y0 = ctr(2);
end


ct = cos(t);
st = sin(t);

xs = x0+r*ct;
ys = y0+r*st;
dxs = -r*st;
dys = r*ct;
d2xs = -r*ct;
d2ys = -r*st;

fvals(:,1) = xs;
fvals(:,2) = ys;
fvals(:,3) = dxs;
fvals(:,4) = dys;
fvals(:,5) = d2xs;
fvals(:,6) = d2ys;

end

