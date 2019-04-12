%EXAMPLE_BARBELL_001_BIG_PLOTS
%
% batch make images of vorticity (this will take a while, an adaptive
% routine is used which ensures accuracy).

% addpath('../src','../test','../../mwrap');

filein = 'example_barbell_001.mat';
fileout = 'example_barbell_001_big_plots.mat';
% addpath('~/Dropbox/MATLAB/chebfun')
% addpath(genpath('~/Dropbox/MATLAB/FLAM'))

addpaths_loc();

load(filein);

%% 

nints = length(detchebs);

rtsc = cell(nints,1);

for i = 1:nints
    rtsc{i} = roots(detchebs{i},'complex');
end

nrts = 0;
for i = 1:nints
    nrts = nrts + length(rtsc{i});
end

rts = zeros(nrts,1);
ind = 1;
for i = 1:nints
    rts(ind:(ind+length(rtsc{i})-1)) = rtsc{i};
    ind = ind + length(rtsc{i});
end

%%

nxdir = 100;

xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
offsetx = 0.5*(xmax-xmin)/(nxdir-1);
offsety = 0.5*(ymax-ymin)/(nxdir-1);
xgrid = linspace(xmin+offsetx,xmax-offsetx,nxdir);
ygrid = linspace(ymin+offsety,ymax-offsety,nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

%

start = tic; in = chunkerin(chunker,targs); toc(start)

hold off
scatter(targs(1,in),targs(2,in),'go');
hold on
scatter(targs(1,~in),targs(2,~in),'rx');
scatter(xs(:),ys(:),'bo');

%%

opts.usesmooth=2;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-8,'AbsTol',1.0e-8};
opts.gausseps = 1e-7;
targsin=targs(:,in);
intparams = []; intparams.intorder = 16;
ss = cell(length(rts),1);
vorts = cell(length(rts),1);
xnulls = cell(length(rts),1);

opts.FLAM=1;

parfor i = 1:length(rts)
    fprintf('iter %d\n',i);
    zk = real(rts(i));
    
    
    start = tic; [d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts);
    toc(start)
    nsys = 2*chunker.nch*chunker.k;
    start = tic; xnull = rskelf_nullvec(F,nsys,[],1e-5,40,4); toc(start)

    ss{i} = norm(rskelf_mv(F,xnull));
    
    ss{i}

    xnulls{i} = xnull;
    
    vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
    ndimsv = [1 2];
    start = tic; 
    vortin = chunkerintkern(chunker,vortkern,ndimsv,xnull,targsin,opts); 
    toc(start)

    vorts{i} = nan(size(xx));
    vorts{i}(in) = vortin;

end

save(fileout,'vorts','in','xx','yy','ss','xnulls','rts');
