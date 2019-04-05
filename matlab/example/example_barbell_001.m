
%EXAMPLE_BARBELL_001
%
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D for a domain on the interval [0.5,6.5]
% save chebfun, domain, and some settings to a file

filebase = 'example_barbell_001'; 
timeref = datestr(now,'_yyyymmdd_HHMMSS');
fileout = [filebase, timeref, '.mat'];


chebabs = cell(12,1);
for j = 1:length(chebabs)
    chebabs{j} = 0.5*[1.0*j 1.0*(j+1)];
end

max_rzk = chebabs{end}(end); lam = 2*pi/max_rzk;

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../example')
addpath('../../mwrap')
addpath('~/Dropbox/MATLAB/chebfun/')
addpath(genpath('~/Dropbox/MATLAB/FLAM'))


opts = [];
opts.autowidths = false;
opts.autowidthsfac = 0.1;


verts = barbell_example_001(6,3);

opts.widths = 0.1*ones(size(verts,2),1);
opts.eps = 1e-12;

chunker1 = chunkcorner(verts,opts);


optsref = [];
optsref.nchmax = 10000;
optsref.maxchunklen = lam;

checkadjinfo(chunker1);

chunkert = chunkrefine(chunker1,optsref);

checkadjinfo(chunkert);

chunkers{1} = chunkert;

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
cs = 0;
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