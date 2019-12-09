%EXAMPLE_MANY_HOLES_005_BIG_PLOTS
%
% batch make images of vorticity 

% addpath('../src','../test','../../mwrap');

filein = 'example_many_holes_005.mat';
fileout = 'example_many_holes_005_big_plots.mat';
% addpath('~/Dropbox/MATLAB/chebfun')
% addpath(genpath('~/Dropbox/MATLAB/FLAM'))

addpaths_loc();

load(filein);

%% 

nints = length(detchebs);

rtsc = cell(nints,1);
frtsc = cell(nints,1);
dfrtsc = cell(nints,1);
dfrtsc2 = cell(nints,1);

nremove_imag = 0;
nremove_real = 0;
nremove = 0;
for i = 1:nints
    a = chebabs{i}(1);
    b = chebabs{i}(2);
    f = detchebs{i};
    df = diff(f);
    fmax = norm(f,'inf');
    rtstemp = roots(f,'complex','norecursion');
    rtstemp_trans = (rtstemp-a)*2/(b-a)-1;
    rho = 1+10*sqrt(p.chebfuneps); aell = (rho+1/rho)/2; bell = (rho-1/rho)/2;
    indkeep = and(real(rtstemp_trans).^2/aell^2 + imag(rtstemp_trans).^2/bell^2 <= 1, ...
        and(real(rtstemp) >= chebabs{i}(1),real(rtstemp) <= chebabs{i}(2)));
    rtstemp(~indkeep);
    
%     indkeep = and(imag(rtstemp_trans) > -sqrt(p.chebfuneps), ...
%         and(real(rtstemp) >= chebabs{i}(1),real(rtstemp) <= chebabs{i}(2)));
%     rtstemp(~indkeep);
%     
%     ind_badimag = and(imag(rtstemp_trans) <= -sqrt(p.chebfuneps),...
%         and(real(rtstemp) >= a,real(rtstemp) <= b));
%     rtstemp(ind_badimag)
%     nremove_imag = nnz(and(imag(rtstemp_trans) <= -sqrt(p.chebfuneps),...
%         and(real(rtstemp) >= a,real(rtstemp) <= b)))+nremove_imag;
%     nremove_real = nnz(or(real(rtstemp) < chebabs{i}(1),real(rtstemp) ...
%         > chebabs{i}(2))) + nremove_real;
%     
    nremove = nnz(~indkeep)+nremove;
    rtstemp = rtstemp(indkeep);
    [~,inds] = sort(real(rtstemp));
    rtsc{i} = rtstemp(inds);
    
    frtsc{i} = f(rtsc{i});
    dfrtsc{i} = df(rtsc{i});
    dfrtsc2{i} = fmax./df(rtsc{i});
end

nrts = 0;
for i = 1:nints
    nrts = nrts + length(rtsc{i});
end

rts = zeros(nrts,1);
frts = zeros(nrts,1);
dfrts = zeros(nrts,1);
dfrts2 = zeros(nrts,1);
ind = 1;
for i = 1:nints
    rts(ind:(ind+length(rtsc{i})-1)) = rtsc{i};
    frts(ind:(ind+length(rtsc{i})-1)) = frtsc{i};
    dfrts(ind:(ind+length(rtsc{i})-1)) = dfrtsc{i};
    dfrts2(ind:(ind+length(rtsc{i})-1)) = dfrtsc2{i};
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
    start = tic; xnull = rskelf_nullvec(F,nsys,4,40); toc(start)

    ss{i} = sqrt(sum(rskelf_mv(F,xnull).^2,1));
    
    ss{i}

    xnulls{i} = xnull;
    
    vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
    ndimsv = [1 2];
    start = tic; 
    vortin = chunkerintkern(chunker,vortkern,ndimsv,xnull(:,1),targsin,opts); 
    toc(start)

    vorts{i} = nan(size(xx));
    vorts{i}(in) = vortin;

end

save(fileout,'vorts','in','xx','yy','ss','xnulls','rts');

