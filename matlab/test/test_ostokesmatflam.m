
%TEST_OSTOKESMATFLAM
%
% This file tests using the matrix routines with
% FLAM (Ken's Fast Linear Algebra in Matlab)

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

rnorms = chunknormals(chunker);

% build matrix

zk = 1.0;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ncomp = length(nchs);

intparams.intorder = chunker.k;

%


%% full dense

fprintf('\nUsing dense routines...\n\n')
fprintf('forming full matrix ...\n')
start = tic; sysmat2 = ostokesmat(chunker,zk,cs,cd,intparams); 
t1 = toc(start);
fprintf('time %5.2e\n',t1)

fprintf('computing determinant...\n')
start = tic; d = det(sysmat2); t2 = toc(start);
fprintf('time %5.2e dense determinant\n',t2)

fprintf('\ntime %6.3e total dense\n',t1+t2)

%% FLAM
npts = length(xs);
xflam = zeros(2,2*npts);
xflam(:,1:2:1+2*(npts-1)) = chunker.chunks(1:2,:);
xflam(:,2:2:2+2*(npts-1)) = chunker.chunks(1:2,:);

xnorm = zeros(2,2*npts);
xnorm(:,1:2:1+2*(npts-1)) = rnorms(1:2,:);
xnorm(:,2:2:2+2*(npts-1)) = rnorms(1:2,:);

rnorm1 = rnorms(:);

whts = chunkwhts(chunker);
whtsflam = repmat((whts(:)).',2,1);

p = 64;
[proxy,pnorm,pw] = proxy_square_pts(p);

k = chunker.k;

kern = @(s,t,sn,tn,sw,tw,slf) ostokes_pxy_kern(s,t,sn,tn,sw,tw,slf, ...
    zk,cs,cd);
pxyfun = @(x,slf,nbr,l,ctr)pxyfun_sq(kern,proxy,pnorm,pw, ...
    x,xnorm,whtsflam,slf,nbr,l,ctr);
matfun = @(i,j) ostokes_matfun_wpre(xflam,xnorm,whtsflam,i,j,zk,cs,cd,...
    k,nchs,ncomp,stokestd,stokesinds);

fprintf('\n Using FLAM ...\n\n')
fprintf('precomputing tri-diagonal part ...\n')
start = tic; [stokestd,stokesinds] = ostokesmat_pre(chunker,nchs, ...
    zk,cs,cd,intparams);
t3 = toc(start);
fprintf('time %5.2e\n',t3)

fprintf('computing skeletonization w/ precomp tri-diagonal ...\n')
start = tic; F_wpre = rskelf(matfun,xflam,200,1e-14,pxyfun); 
t4 = toc(start);
fprintf('time %5.2e\n',t4)

fprintf('computing determinant ...\n')
start = tic; drskelf = exp(rskelf_logdet(F_wpre)); t5 = toc(start);
fprintf('time %5.2e\n\n',t5)

fprintf('time %6.3e total rskelf\n',t3+t4+t5)

%% TEST

ntest = 5;
[nsys,~] = size(sysmat2);
mutest = randn(nsys,ntest);

ytest = sysmat2*mutest;
y_wpre = rskelf_mv(F_wpre,mutest);
%
err_mv = norm(ytest-y_wpre,'fro')/norm(ytest,'fro');

fprintf('\nerrors...\n\n')
fprintf('error %5.2e matvec\n',err_mv)

err_det = abs(d-drskelf)/abs(d);

fprintf('error %5.2e determinant\n',err_det)




