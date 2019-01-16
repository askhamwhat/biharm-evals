
function [fvals] = fcurve(t,varargin)

narms = 5;
amp = 0.3;
if nargin > 1 && ~isempty(varargin{1})
    narms = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    amp = varargin{2};
end

fvals = zeros(length(t),6);
ct = cos(t);
st = sin(t);
cnt = cos(narms*t);
snt = sin(narms*t);

xs = (1+amp*cnt).*ct;
ys = (1+amp*cnt).*st;
dxs = -(1+amp*cnt).*st-narms*amp*snt.*ct;
dys = (1+amp*cnt).*ct-narms*amp*snt.*st;
d2xs = -dys-narms*amp*(narms*cnt.*ct-snt.*st);
d2ys = dxs-narms*amp*(narms*cnt.*st+snt.*ct);

fvals(:,1) = xs;
fvals(:,2) = ys;
fvals(:,3) = dxs;
fvals(:,4) = dys;
fvals(:,5) = d2xs;
fvals(:,6) = d2ys;

end

