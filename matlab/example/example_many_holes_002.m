
%EXAMPLE_MANY_HOLES_002
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D for a domain on the interval [0.1,10]
% save chebfun, domain, and some settings to a file

filebase = 'example_many_holes_002'; 
timeref = datestr(now,'_yyyymmdd_HHMMSS');
fileout = [filebase, timeref, '.mat'];

a = 0.1; b = 10.0; chebab = [a b];

max_rzk = b; lam = 2*pi/max_rzk;

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../example')
addpath('../../mwrap')

cparams.eps = 1.0e-3;
cparams.nchmax = 100000;
cparams.maxchunklen = lam;
cparams.nover = 2;
narms = 3;
amp = 0.25;
W = 3.0; H = 2.0; p = 4.0;
chunkers = {};
chunkert = chunkfunc(@(t) smoothbox(t,W,H,p),cparams);
chunkers{1} = chunkert;

nw = 4;
nh = 3;

ind = 1;
for ii = 1:nw
    for jj = 1:nh
        phi = rand()*2*pi;
        scale = 0.20;
        ctr = [-(nw-1)/2.0 + ii-1; -(nh-1)/2.0 + jj-1];
        chunkert = chunkfunc(@(t) starfish(t,narms,amp,ctr,phi,scale),cparams);
        chunkert = chunkreverse(chunkert);
        chunkert = chunksort(chunkert);
        ind = ind+1;
        chunkers{ind} = chunkert;
        
    end
end
        
[chunker,nchs] = chunkermerge(chunkers);

%

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

rnorms = chunknormals(chunker);

intparams.intorder = chunker.k;

hold off
quiver(xs,ys,rnorms(1,:).',rnorms(2,:).')
hold on
scatter(xs,ys)

%%

% build matrix

cd = -2.0 + 1i*0.0;
cs = 0.0 - 1i*2.0;
%csfun = @(zk) -2.0*1i*zk;
ncomp = length(nchs);

%

opts = []; opts.FLAM = 1; opts.verb = true;
detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
    cd,opts);   

p = chebfunpref; p.chebfuneps = 1.0e-12; p.domain = chebab; 
p.splitting=0; p.maxLength = 2049;
start = tic; detcheb = chebfun(detfun,p); t1 = toc(start);
fprintf('%5.2e time for chebfun build\n',t1)

save(fileout,'detcheb','chunker','nchs','cs','cd','opts','p','t1');