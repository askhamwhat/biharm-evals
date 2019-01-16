
%TESTHBHDIRMAT
%
% This file tests whether or not chunks 
% lib is wrapped properly

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 0;
narms = 5;
amp = 0.5;
tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc

chunker.nch
chunker.k*chunker.nch

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

tic; mat1 = chunknormonesmat(chunker); toc
tic; mat2 = normalonesmatwrap(chunker); toc

rnorms = chunknormals(chunker);
rnx = rnorms(1,:,:); rnx = rnx(:);
rny = rnorms(2,:,:); rny = rny(:);

figure(1)
clf
scatter(xs,ys)
axis equal
hold on
quiver(xs,ys,rnx,rny)

%%

zk = 10.0;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ccs = zeros(2,1);

intparams.intorder = chunker.k;

start = tic; sysmat1 = hbhdirmat(chunker,nchs,ccs,zk,cs,cd,intparams); toc(start)

%%

wgeo = chunkpack(chunker);
start = tic; sysmat2 = zhbhstokesmatbuild(zk,wgeo,nchs,ccs,cs,cd); toc(start)

%%
hold off
spy(abs(sysmat1-sysmat2)>1.0e-12)   