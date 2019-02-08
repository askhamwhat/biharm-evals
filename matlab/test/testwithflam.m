
%TESTWITHFLAM
%
% This file tests using the matrix routines with
% FLAM (slow version of matrix build for now).

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

rnorms = chunknormals(chunker);

% build matrix

zk = 1;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ncomp = length(nchs);
ccs = zeros(2,1);

intparams.intorder = chunker.k;

%

start = tic; sysmat1 = hbhdirmat(chunker,nchs,ccs,zk,cs,cd,intparams); toc(start)

%

start = tic; yylap = hbhdirmat_lap_pre(chunker,nchs,intparams); toc(start)

% 

start = tic; [stokestd,stokesinds,sysmatbl,sysmatbr,sysmattr] = ...
    hbhdirmat_pre(chunker,nchs,ccs,zk,cs,cd,intparams,yylap); toc(start)

%% FLAM
[ntot,~] = size(sysmat1);
npts = length(xs);
xflam = zeros(2,ntot);
xflam(:,1:2:1+2*(npts-1)) = chunker.chunks(1:2,:);
xflam(:,2:2:2+2*(npts-1)) = chunker.chunks(1:2,:);
xflam(:,2*npts+1:ntot) = ccs;

xnorm = zeros(2,ntot);
xnorm(:,1:2:1+2*(npts-1)) = rnorms(1:2,:);
xnorm(:,2:2:2+2*(npts-1)) = rnorms(1:2,:);
xnorm(:,2*npts+1:ntot) = ccs;

rnorm1 = rnorms(:);

whts = chunkwhts(chunker);
whtsflam = repmat((whts(:)).',2,1);
whtsflam = [whtsflam(:); zeros(ncomp,1)];
whtsflam = whtsflam(:);

p = 64;
[proxy,pnorm] = proxy_square_pts(p);

k = chunker.k;

kern = @(s,t,sn,tn,slf) kern_hs(s,t,sn,tn,slf,npts,zk,cs,cd,ccs,whtsflam);
pxyfun = @(x,slf,nbr,l,ctr)pxyfun_sq(kern,proxy,pnorm,x,xnorm, ... 
    slf,nbr,l,ctr);

kernr = @(s,t,sn,tn,slf) pxykerndirr(s,t,sn,tn,slf,npts,zk,cs,cd,...
    whtsflam) ;
kernc = @(s,t,sn,tn,slf) pxykerndirc(s,t,sn,tn,slf,npts,zk,cs,cd,...
    whtsflam) ;
pxyfunr = @(rc,rx,cx,slf,nbr,l,ctr) pxyfun_sqr(rc,rx,cx,slf,nbr, ...
    l,ctr,kernr,kernc,xnorm,xnorm,proxy,pnorm);

matfun = @(i,j) hbhdir_matfun_wpre(xflam,xnorm,i,j,zk,cs,cd,whtsflam,k,...
    nchs,ncomp,stokestd,stokesinds,sysmatbl,sysmatbr,sysmattr,rnorm1);

%%

opts.rdpiv='Q';
start = tic; F_nopre = rskelfr(sysmat1,xflam,xflam,200,1e-14,[],opts); toc(start)


%%

opts.nc = [2*chunker.k*chunker.nch+1];
opts.nc = [];
start = tic; F_wpre = rskelf(matfun,xflam,200,1e-14,[],opts); 
toc(start)

%%

% test

mutest = randn(ntot,5);
ytrue = sysmat1*mutest;
yflam_nopre = rskelfr_mv(F_nopre,mutest,'N');
yflam_wpre = rskelf_mv(F_wpre,mutest,'N');

norm(ytrue-yflam_nopre,'fro')/norm(ytrue,'fro')
norm(ytrue-yflam_wpre,'fro')/norm(ytrue,'fro')

%% 

sysmat2 = rskelfr_mv(F_wpre,eye(ntot));

%%

norm(sysmat2-sysmat1,'fro')
spy(abs(sysmat2-sysmat1) > 1e-5)

%%

iall = 1:ntot; jall = 1:ntot;
start = tic; sysmat4 = matfun(iall,jall); toc(start)

%%

norm(sysmat4-sysmat1,'fro')
nnz(abs(sysmat4-sysmat1) > 1e-10)

function Kpxy = pxykerndirc(s,t,sn,tn,slf,npts,zk,cs,cd,swhts) 
%TODO handle multiply connected (log terms)
iin = (slf <= 2*npts);
slfin = slf(iin);
Ktemp = helmstokessubmat(zk,s(:,iin),t,sn(:,iin),cs,cd);
indcs = 2*(1:length(slfin))-mod(slfin,2);

Kpxy = zeros(2*size(t,2),length(slf));
Kpxy(:,iin) = Ktemp(:,indcs);

[m,~] = size(Kpxy);
iout = (slf > 2*npts);
Kpxy(:, iout) = zeros(m,sum(iout));

end

function Kpxy = pxykerndirr(s,t,sn,tn,slf,npts,zk,cs,cd,swhts) 
    iin = (slf <= 2*npts);
    slfin = slf(iin);
    Ktemp = helmstokessubmat(zk,s,t(:,iin),sn,cs,cd);
    indcs = 2*(1:length(slfin))-mod(slfin,2);
    
    Kpxy = zeros(length(slf),size(s,2)*2);
    Kpxy(iin,:) = Ktemp(indcs,:);
    
    iout = (slf > 2*npts);
    Kpxy(iout,:) = randn(sum(iout),size(s,2)*2);
    
end

function [Kpxy] = kern_hs(s,t,sn,tn,slf,npts,zk,cs,cd,ccs,whtsflam)
% TODO HANDLE MULTIPLY CONNECTED
iin = (slf <= 2*npts);
slfin = slf(iin);
Ktemp = helmstokessubmat(zk,s(:,slfin),t,sn(:,slfin),cs,cd);
indcs = 2*(1:length(slfin))-mod(slfin,2);

Kpxy = zeros(2*size(t,2),length(slf));
Kpxy(:,iin) = Ktemp(:,indcs);

[m,~] = size(Kpxy);
iout = (slf > 2*npts);
Kpxy(:, iout) = zeros(m,sum(iout));

whtslf = whtsflam(slf);
Kpxy = bsxfun(@times,Kpxy,(whtslf(:)).');

end

function K = hbhdir_matfun_wpre(s,sn,i,j,zk,cs,cd,whtsflam,k,nchs,ncomp, ...
    stokestd,stokesinds,sysmatbl,sysmatbr,sysmattr,rnorm1)
%% evaluate system matrix by entry index utility function
% 
% this function primarily calls the FORTRAN based submat evaluator
% but makes some attempts at efficiency%
%

nch = sum(nchs);
ngeo = nch*k;

% find rows and columns which are not part of stokes top left block

iicomp = i > 2*ngeo;
icomp = i(iicomp);
iistokes = i <= 2*ngeo;
ijcomp = j > 2*ngeo;
jcomp = j(ijcomp);
ijstokes = j <= 2*ngeo;

% this first part checks for column pairs in stokes part which come
% from the same source point (j(ii) = j(ii)+1) for such indices)

ijodd = and(mod(j,2) == 1,j <= 2*ngeo); jodd = j(ijodd);
[jeveninc,ijeveninc,ijoddinc] = intersect(j(:),jodd(:)+1); 
[jevenout,ijout] = setdiff(j,[jodd(:);jeveninc(:);jcomp(:)]);

% ditto for row pairs

iiodd = and(mod(i,2) == 1,i <= 2*ngeo); iodd = i(iiodd);
[ieveninc,iieveninc,iioddinc] = intersect(i(:),iodd(:)+1); 
[ievenout,iiout] = setdiff(i,[iodd(:);ieveninc(:);icomp(:)]);

Kroddcodd = helmstokessubmat(zk,s(:,jodd),s(:,iodd),sn(:,jodd), ...
    cs,cd);
Kroddcout = helmstokessubmat(zk,s(:,jevenout),s(:,iodd), ...
    sn(:,jevenout),cs,cd);
Kroutcodd = helmstokessubmat(zk,s(:,jodd),s(:,ievenout), ...
    sn(:,jodd), cs,cd);
Kroutcout= helmstokessubmat(zk,s(:,jevenout),s(:,ievenout), ...
    sn(:,jevenout), cs,cd);

K = zeros(length(i),length(j));

% grab everything from odd - odd matrix (have to account for the even
% pairs of these indices, if present)

K(iiodd,ijodd) = Kroddcodd(1:2:end,1:2:end);
K(iieveninc,ijodd) = Kroddcodd(2*iioddinc,1:2:end);
K(iiodd,ijeveninc) = Kroddcodd(1:2:end,2*ijoddinc);
K(iieveninc,ijeveninc) = Kroddcodd(2*iioddinc,2*ijoddinc);

% others are a little simpler

K(iiodd,ijout) = Kroddcout(1:2:end,2:2:end);
K(iieveninc,ijout) = Kroddcout(2*iioddinc,2:2:end);

K(iiout,ijodd) = Kroutcodd(2:2:end,1:2:end);
K(iiout,ijeveninc) = Kroutcodd(2:2:end,2*ijoddinc);

K(iiout,ijout) = Kroutcout(2:2:end,2:2:end);

wj = whtsflam(j);
K = bsxfun(@times,K,wj(:).');

% relatively efficient replacement of entries which require special
% quadrature (such entries should have been precomputed)

nrow = size(s,2);
assert(nrow == 2*ngeo + ncomp);
lininds = bsxfun(@plus,i(:),2*ngeo*(j(:)-1).'); % pretend it's just 
                        % stokes matrix for reference purposes
lininds(iicomp,:) = -1; % ignores non-stokes entries(won't be matched)
lininds(:,ijcomp) = -1;
lininds = lininds(:);

istokes = i(iistokes);
stokesindsi= stokesinds(istokes,:); % iicomp rows would throw error
stokeslinindsi = bsxfun(@plus,istokes(:),(stokesindsi-1)*2*ngeo);
stokeslinindsi = stokeslinindsi(:);
stokestdi = stokestd(istokes,:);

[~,ik,itdi] = intersect(lininds,stokeslinindsi);

K(ik) = stokestdi(itdi);

% replace other non-stokes parts with precomputes

K(iistokes,ijcomp) = sysmattr(i(iistokes),j(ijcomp)-2*ngeo);
K(iicomp,ijstokes) = sysmatbl(i(iicomp)-2*ngeo,j(ijstokes));
K(iicomp,ijcomp) = sysmatbr(i(iicomp)-2*ngeo,j(ijcomp)-2*ngeo);

% finally add onemat correction for stokes type entries

istokes = i(iistokes); jstokes = j(ijstokes);
ri = rnorm1(istokes); rj = rnorm1(jstokes);
onefix = bsxfun(@times,ri(:),rj(:).');
wj = whtsflam(jstokes);
onefix = bsxfun(@times,onefix,(wj(:)).');
K(iistokes,ijstokes) = K(iistokes,ijstokes)+onefix;

end
