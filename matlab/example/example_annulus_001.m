
%EXAMPLE_ANNULUS_001
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D-2iS for a domain on the intervals [j,j+1] j=1,...,8
% save chebfun, domain, and some settings to a file

filebase = 'example_annulus_001'; 
timeref = datestr(now,'_yyyymmdd_HHMMSS');
fileout = [filebase, '_5_convplots.mat'];

% annulus params

r1 = 1.0;
r2 = 1.7;


ncells = 2;

chebabs = cell(ncells,1);
chebabs{1} = [13 14];
chebabs{2} = [14 15];

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../example')
addpath('../../mwrap')

cparams.eps = 1.0;
cparams.nchmax = 100000;
cparams.nover = 1;
narms = 3;
amp = 0.25;
W = 3.0; H = 2.0; p = 4.0;

chunkers = {};
nchi = 8;
ncho = ceil(nchi/r1*r2)+1;

chunkert = circle_chunks(ncho,r2);
chunkers{1} = chunkert;

nw = 1;
nh = 1;

ind = 2;
chunkert = circle_chunks(nchi,r1);
chunkert = chunkreverse(chunkert);
chunkert = chunksort(chunkert);
chunkers{ind} = chunkert;
        
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
    
for j = 1:length(chebabs)
    fprintf('running on interval [ %5.2e %5.2e ]\n',chebabs{j}(1), ...
        chebabs{j}(2))
    start = tic; detchebs{j} = ...
        chebfun(detfun,chebabs{j},p); 
    t1s{j} = toc(start);
    fprintf('%5.2e time for chebfun build\n',t1s{j})
end



%% Now compute everything using just double layer
cd = -2.0 + 1i*0.0;
cs = 0;
ncomp = length(nchs);

opts = []; opts.FLAM = 1; opts.verb = true;
detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
    cd,opts);   

p = chebfunpref; p.chebfuneps = 1.0e-13;
p.splitting=0; p.maxLength = 257;

t1s = cell(size(chebabs));
detchebs2 = cell(size(chebabs));
    
for j = 1:length(chebabs)
    fprintf('running on interval [ %5.2e %5.2e ]\n',chebabs{j}(1), ...
        chebabs{j}(2))
    start = tic; detchebs2{j} = ...
        chebfun(detfun,chebabs{j},p); 
    t1s{j} = toc(start);
    fprintf('%5.2e time for chebfun build\n',t1s{j})
end


save(fileout,'detchebs','detchebs2','chunker','nchs','cs','cd','opts','p','t1s',...
    'chebabs','r1','r2');