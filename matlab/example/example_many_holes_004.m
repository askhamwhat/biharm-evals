
%EXAMPLE_MANY_HOLES_004
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D-2iS for a domain on finer intervals
% save chebfun, domain, and some settings to a file

filebase = 'example_many_holes_004'; 
timeref = datestr(now,'_yyyymmdd_HHMMSS');
fileout = [filebase, timeref, '.mat'];

chebabs = cell(16,1);
for j = 1:8
    chebabs{j} = 3.0 + 0.5*[1.0*(j-1) 1.0*(j)];
end

for j = 9:16
    chebabs{j} = 7 + 0.25*[1.0*(j-9) 1.0*(j-8)];
end

max_rzk = chebabs{end}(end); lam = 2*pi/max_rzk;

%%

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
amp = 0.125;
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