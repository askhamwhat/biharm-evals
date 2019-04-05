
function chnkr = chunksplit1(chnkr,ich,opts)
%CHUNKSPLIT1 this routine takes the list of all chunks and splits one in
%      half with respect to arclength.
%   Input
%        ich - the chunk number to split

thresh=1.0e-10;
nitermax = 1000;

if nargin < 3
    opts = [];
end

if isfield(opts,'thresh')
    thresh = opts.thresh;
end
if isfield(opts,'nitermax')
    nitermax = opts.nitermax;
end


k = chnkr.k;
nch = chnkr.nch;

[x, w, u, ~] = legeexps(k);

%  first construct dsdt

dsdt = sqrt(chnkr.ders(1,:,ich).^2+chnkr.ders(2,:,ich).^2)*chnkr.hs(ich);
dsdt = dsdt(:);
rltot = dot(dsdt,w);

cdsdt = u*dsdt;

t1=0;
rlhalf=rltot/2;

% use Newton to find t such that length(-1,t) = length(t,1) = rl/2

ifdone = 0;
for ijk = 1:nitermax
    ts = -1+(t1+1)*(x + 1)/2.0;
    ws = (t1+1)*w/2.0;

    vals = legeexevvec(ts,cdsdt);
    rl1 = dot(vals,ws);
    val = legeexevvec(t1,cdsdt);
    err=rl1-rlhalf;
    if (abs(err) < thresh); ifdone=ifdone+1; end
    if (ifdone >= 3) 
        break;
    end
    t1=t1-(rl1-rlhalf)/val;
end

if (ifdone < 3); warning('did not converge'); end

% new points in parameter space

ts1 = -1+(t1+1)*(x+1)/2.0;
ts2 = t1+(1-t1)*(x+1)/2.0;

% evaluate the new values of chunks, ders, ders2 and 
% update nch, adjs, hs

xy = chnkr.chunks(:,:,ich);
dxy = chnkr.ders(:,:,ich);
d2xy = chnkr.ders2(:,:,ich);

cxy = u*(xy.');
cdxy = u*(dxy.');
cd2xy = u*(d2xy.');

i1=chnkr.adjs(1,ich);
i2=chnkr.adjs(2,ich);

ndim = 2;

chunks_1 = zeros(ndim,k,1);
ders_1 = zeros(ndim,k,1);
ders2_1 = zeros(ndim,k,1);
chunks_2 = zeros(ndim,k,1);
ders_2 = zeros(ndim,k,1);
ders2_2 = zeros(ndim,k,1);

for i = 1:ndim
    x = legeexevvec(ts1,cxy(:,i));
    dx = legeexevvec(ts1,cdxy(:,i));
    d2x = legeexevvec(ts1,cd2xy(:,i));
    chunks_1(i,:,1) = x;
    ders_1(i,:,1) = dx;
    ders2_1(i,:,1) = d2x;
    x = legeexevvec(ts2,cxy(:,i));
    dx = legeexevvec(ts2,cdxy(:,i));
    d2x = legeexevvec(ts2,cd2xy(:,i));
    chunks_2(i,:,1) = x;
    ders_2(i,:,1) = dx;
    ders2_2(i,:,1) = d2x;
end

hsold=chnkr.hs(ich);

% update chnkr

chnkr.nch=nch+1;

chnkr.hs(ich) = hsold*(t1+1)/2;
chnkr.hs(nch+1) = hsold*(1-t1)/2;

chnkr.chunks(:,:,ich) = chunks_1;
chnkr.chunks(:,:,nch+1) = chunks_2;
chnkr.ders(:,:,ich) = ders_1;
chnkr.ders(:,:,nch+1) = ders_2;
chnkr.ders2(:,:,ich) = ders2_1;
chnkr.ders2(:,:,nch+1) = ders2_2;

chnkr.adjs(2,ich)=nch+1;
chnkr.adjs(1,nch+1)=ich;
chnkr.adjs(2,nch+1)=i2;
chnkr.adjs(1,i2)=nch+1;

end
