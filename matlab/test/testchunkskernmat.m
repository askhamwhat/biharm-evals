
%TESTCHUNKSKERNMAT
%
% This file tests whether or not chunks 
% lib is wrapped properly

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc

chunker.nch

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

figure(1)
clf
scatter(xs,ys)
axis equal

fkern = @(s,t,sn,tn) glapkern(s,t,sn,tn,'D');
ndims(1) = 1; ndims(2) = 1;
intparams.intorder = chunker.k;
tic; D = chunkskernelmat(chunker,fkern,ndims,intparams); toc

%
sysmat = eye(size(D))-2*D;

tic; d1 = det(sysmat); toc

% hifie
npts = length(xs);
xhifie = chunker.chunks(1:2,:);

%
tic; F = hifie2(sysmat,xhifie,10,1e-10); toc
tic; d2 = exp(hifie_logdet(F)); toc

d1
d2

%scatter(real(e),imag(e))