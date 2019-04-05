function [chunker] = circle_chunks(nch,varargin)

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

k = 16;
chunker.k = 16;
chunker.nch = nch;

chunks = zeros(2,k*nch);
ders = zeros(2,k*nch);
ders2 = zeros(2,k*nch);
adjs = zeros(2,nch);
hs = zeros(nch,1);

ts = legpts(k);

h = 2*pi/nch;
for i=1:nch
    tstart = (i-1)*h;
    hs(i) = h/2;
    adjs(1,i) = i-1;
    adjs(2,i) = i+1;
    for j=1:k
        t = tstart + (ts(j)+1)/2*h;
        ct = cos(t);
        st = sin(t);
        chunks(1,(i-1)*k+j) = x0 + r*ct;
        chunks(2,(i-1)*k+j) = y0+r*st;
        ders(1,(i-1)*k+j) = -r*st;
        ders(2,(i-1)*k+j) = r*ct;
        ders2(1,(i-1)*k+j) = -r*ct;
        ders2(2,(i-1)*k+j) = -r*st;
   end
end

adjs(1,1) = nch;
adjs(2,nch) = 1;

chunker.chunks = chunks;
chunker.ders = ders;
chunker.ders2 = ders2;
chunker.adjs = adjs;
chunker.hs = hs;


end