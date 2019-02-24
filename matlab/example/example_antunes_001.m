
%EXAMPLE_ANTUNES_001
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D for the Antunes domain on the interval [0.1,10]
% save chebfun, domain, and some settings to a file

fileout = 'example_antunes_001.mat';

a = 0.1; b = 10.0; chebab = [a b];

max_rzk = b; lam = 2*pi/max_rzk;

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../example')
addpath('../../mwrap')

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.maxchunklen = lam;
cparams.nover = 2;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) antunes(t),cparams);

%%

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

rnorms = chunknormals(chunker);

intparams.intorder = chunker.k;

% build matrix

cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ncomp = length(nchs);

%

opts = []; opts.FLAM = 1; opts.verb = true;
detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs,cd,opts);

p = chebfunpref; p.chebfuneps = 1.0e-12; p.domain = chebab; 
p.splitting=0;
start = tic; detcheb = chebfun(detfun,p); t1 = toc(start);
fprintf('%5.2e time for chebfun build\n',t1)

save(fileout,'detcheb','chunker','nchs','cs','cd','opts','p','t1');