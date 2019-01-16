
%TESTWITHFLAM
%
% This file tests using the matrix routines with
% FLAM (slow version of matrix build for now).

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc

chunker.nch
chunker.k*chunker.nch

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

%tic; mat1 = chunknormonesmat(chunker); toc
%tic; mat2 = normalonesmatwrap(chunker); toc

rnorms = chunknormals(chunker);
rnx = rnorms(1,:,:); rnx = rnx(:);
rny = rnorms(2,:,:); rny = rny(:);

%figure(1)
%clf
%scatter(xs,ys)
%axis equal
%hold on
%quiver(xs,ys,rnx,rny)

% build matrix

zk = 0.1;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ccs = zeros(2,1);

intparams.intorder = chunker.k;

%

start = tic; sysmat1 = hbhdirmat(chunker,nchs,ccs,zk,cs,cd,intparams); toc(start)

%%

%fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
%ndims(1) = 2; ndims(2) = 2;
%tic; [stokestd,inds] = chunkskernelmattd(chunker,fkern,ndims,intparams); toc

%%

start = tic; yylap = hbhdirmat_lap_pre(chunker,nchs,intparams); toc(start)

%% 

start = tic; stokestd = hbhdirmat_pre(chunker,nchs,ccs,zk,cs,cd,intparams,yylap); toc(start)

%% hifie
[ntot,~] = size(sysmat1);
npts = length(xs);
xhifie = zeros(2,ntot);
xhifie(:,1:2:1+2*(npts-1)) = chunker.chunks(1:2,:);
xhifie(:,2:2:2+2*(npts-1)) = chunker.chunks(1:2,:);
xhifie(:,2*npts+1:ntot) = ccs;

xnorm = zeros(2,ntot);
xnorm(:,1:2:1+2*(npts-1)) = rnorms(1:2,:);
xnorm(:,2:2:2+2*(npts-1)) = rnorms(1:2,:);
xnorm(:,2*npts+1:ntot) = ccs;

p = 64;
[proxy,pnorm] = proxy_square_pts(p);
[proxyc,pnormc] = proxy_circ_pts(p);

kern = @(s,t,sn,tn,slf) kern_hs(s,t,sn,tn,slf,npts,zk,cs,cd);
pxyfun = @(x,slf,nbr,l,ctr)pxyfun_sq(kern,proxy,pnorm,x,xnorm, ... 
    slf,nbr,l,ctr);
%pxyfun = @(x,slf,nbr,l,ctr)pxyfun_circ(kern,proxyc,pnormc,x,xnorm, ... 
%    slf,nbr,l,ctr);

%%

start = tic; F = hifie2(sysmat1,xhifie,80,1e-12,pxyfun); toc(start)
start = tic; F = rskelf(sysmat1,xhifie,80,1e-12,pxyfun); toc(start)

%%

start = tic; detF = rskelf_logdet(F); toc(start)


%%



function [Kpxy] = matfun(s,sn,slf,npts,inds,zk,cs,cd)

Kpxy = helmstokessubmat(zk,s(:,slf),t,sn(:,slf),cs,cd);
indcs = 2*(1:length(slf))-mod(slf,2);

Kpxy = Kpxy(:,indcs);
[m,~] = size(Kpxy);
iout = (slf > 2*npts);
Kpxy(:, iout) = randn(m,sum(iout));

end


function [Kpxy] = kern_hs(s,t,sn,tn,slf,npts,zk,cs,cd)

Kpxy = helmstokessubmat(zk,s(:,slf),t,sn(:,slf),cs,cd);
indcs = 2*(1:length(slf))-mod(slf,2);

Kpxy = Kpxy(:,indcs);
[m,~] = size(Kpxy);
iout = (slf > 2*npts);
Kpxy(:, iout) = randn(m,sum(iout));

end
