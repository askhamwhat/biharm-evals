
%EXAMPLE_ANTUNES_001
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D for the Antunes domain on the interval [0.1,10]
% save chebfun, domain, and some settings to a file

fileout = 'example_antunes_001.mat';

chebabs = cell(12,1);
for j = 1:length(chebabs)
    chebabs{j} = [1.0*j 1.0*(j+1)];
end

max_rzk = chebabs{end}(end); lam = 2*pi/max_rzk;

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../example')
addpath('../../mwrap')
addpath('~/Dropbox/MATLAB/chebfun')
addpath(genpath('~/Dropbox/MATLAB/FLAM'))
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

scatter(xs,ys)

%%

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
detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
    cd,opts);   

p = chebfunpref; p.chebfuneps = 1.0e-13;
p.splitting=0; p.maxLength = 257;

t1s = cell(size(chebabs));
detchebs = cell(size(chebabs));
    
parfor j = 1:length(chebabs)
    fprintf('running on interval [ %5.2e %5.2e ]\n',chebabs{j}(1), ...
        chebabs{j}(2))
    start = tic; detchebs{j} = ...
        chebfun(detfun,chebabs{j},p); 
    t1s{j} = toc(start);
    fprintf('%5.2e time for chebfun build\n',t1s{j})
end

save(fileout,'detchebs','chunker','nchs','cs','cd','opts','p','t1s',...
    'chebabs');