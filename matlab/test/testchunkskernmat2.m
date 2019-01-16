
%TESTCHUNKSKERNMAT2
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

zk = 1.0;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
ndims(1) = 2; ndims(2) = 2;
intparams.intorder = chunker.k;
tic; matul = chunkskernelmat(chunker,fkern,ndims,intparams); toc

normonesmat = chunknormonesmat(chunker);
matul = eye(size(matul)) + matul + normonesmat;

%%

wgeo = chunkpack(chunker);

nchs(1) = chunker.nch;
ccs = zeros(2,1);

q1 = cs; q2 = cd;
tic; sysmat = zhbhstokesmatbuild(zk,wgeo,nchs,ccs,q1,q2); toc
