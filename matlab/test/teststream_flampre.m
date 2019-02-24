%% build a recursive skeletonization of the stream matrix 
% using only a precomputed tri-diagonal (everything else on 
% the fly)

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

% geometry parameters and construction

cparams.eps = 1.0e-3;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

ngeo = chunker.nch*chunker.k;

whts = chunkwhts(chunker); whts = whts(:);
whts2 = repmat(whts.',2,1); whts2 = reshape(whts2,2*ngeo,1);

rnorms = chunknormals(chunker); rnorms = reshape(rnorms,2,ngeo);
rnorms2 = repmat(rnorms,2,1); rnorms2 = reshape(rnorms2,2,2*ngeo);
x = chunker.chunks; x = reshape(x,2,ngeo);
x2 = repmat(x,2,1); x2 = reshape(x2,2,2*ngeo);

% problem/representation specification

zk = 0.1;
cd = -2.0 + 1i*0.0;
cs = 0.0 + 1i*0.0;
nchs = ones(1,1);
nchs(1) = chunker.nch;
ccs = zeros(2,1);

% integration parameters

intparams.intorder = chunker.k;

% kernel settings

fkernstream = @(s,t,sn,tn) helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
ndims(1) = 1; ndims(2) = 2;

% proxy fun stuff (rectangular matrix type proxy fun evaluator(s))

np = 64;
[p,pn] = proxy_square_pts(np);

pkernstreamr = @(s,t,sn,tn,slf) pxykernstreamr(s,t,sn,tn,slf,zk,cs,cd,...
    whts) ;
pkernstreamc = @(s,t,sn,tn,slf) pxykernstreamc(s,t,sn,tn,slf,zk,cs,cd,...
    whts2) ;
pxyfunstream = @(rc,rx,cx,slf,nbr,l,ctr) pxyfun_sqr(rc,rx,cx,slf,nbr, ...
    l,ctr,pkernstreamr,pkernstreamc,rnorms,rnorms2,p,pn);

% precompute tridiag for matrix entry evaluator

start = tic; [streamtd,streaminds] = chunkskernelmattd(chunker, ...
    fkernstream,ndims,intparams); t = toc(start);

fprintf('%5.2e : time for tridiag comp\n',t)
% define matrix entry evaluator

matfunstream = @(i,j) stream_matfun_wd(x2,x,rnorms2,[],i,j,zk,cs,cd, ...
    whts2,streamtd,streaminds);

% compute compressed rep of matrix with FLAM's rskelfr

start = tic; F = rskelfr(matfunstream,x,x2,200,1e-12,pxyfunstream); 
t= toc(start);

fprintf('%5.2e : time for rskelfr\n',t)

% compute compressed FMM format of matrix 

matfunstream = @(i,j) stream_matfun_wd(x2,x,rnorms2,[],i,j,zk,cs,cd, ...
    whts2,streamtd,streaminds);

opts.store = 'A';
opts.near = 0;
start = tic; F_fmm = ifmm(matfunstream,x,x2,200,1e-14,pxyfunstream,opts); 
t=toc(start);

fprintf('%5.2e : time for ifmm\n',t)

% compute matrix, straightforward way

start = tic; 
streamfull = chunkskernelmat(chunker,fkernstream,ndims,intparams); 
t = toc(start);
fprintf('%5.2e : time for full eval\n',t)

% test

start = tic; streamfull2 = rskelfr_mv(F,eye(2*ngeo),'n'); t = toc(start);
fprintf('%5.2e : time for mat-mat, rskelfr\n',t)
start = tic; streamfull3 = ifmm_mv(F_fmm,eye(2*ngeo)); t= toc(start);
fprintf('%5.2e : time for mat-mat, ifmm\n',t)


norm(streamfull-streamfull2,'fro')/norm(streamfull,'fro')
norm(streamfull-streamfull3,'fro')/norm(streamfull,'fro')

% test 2

ncomp = length(nchs);

whtsbycomp = zeros(length(whts),ncomp);

k = chunker.k;
ind = 1;
for i = 1:ncomp
    whtsbycomp(ind:ind+nchs(i)*k-1,i) = whts(ind:ind+nchs(i)*k-1);
    ind = ind+nchs(i)*k;
end

yy2 = whtsbycomp.'*streamfull;
yy3 = integrate_stream(chunker,zk,cs,cd,intparams,whtsbycomp);

yy4 = (ifmm_mv(F_fmm,whtsbycomp,[],'T')).';

norm(yy2-yy3,'fro')/norm(yy2,'fro')
norm(yy2-yy4,'fro')/norm(yy2,'fro')

%%

function K = pxykernstreamc(s,t,sn,tn,slf,zk,cs,cd,swhts) 
    K = helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
    indcs = 2*(1:length(slf))-mod(slf,2);
    K = K(:,indcs)*diag(swhts(slf));
end

function K = pxykernstreamr(s,t,sn,tn,slf,zk,cs,cd,swhts) 
    K = helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
end

function K = stream_matfun(s,t,sn,tn,i,j,zk,cs,cd,swhts)

K= helmstokesstreamsubmat(zk,s(:,j),t(:,i),sn(:,j),cs,cd);
indcs = 2*(1:length(j))-mod(j,2);

K= K(:,indcs)*diag(swhts(j));

end

function K = stream_matfun_wd(s,t,sn,tn,i,j,zk,cs,cd,swhts,streamtd,streaminds)

ijodd = find(mod(j,2) == 1); jodd = j(ijodd);
[jeveninc,ijeveninc,ijoddinc] = intersect(j(:),jodd(:)+1); 
[jevenout,ijout] = setdiff(j,[jodd(:);jeveninc(:)]);

Kodd= helmstokesstreamsubmat(zk,s(:,jodd),t(:,i),sn(:,jodd),cs,cd);
Kout = helmstokesstreamsubmat(zk,s(:,jevenout),t(:,i),sn(:,jevenout),cs,cd);
K = zeros(length(i),length(j));
K(:,ijodd) = Kodd(:,1:2:end);
K(:,ijeveninc) = Kodd(:,ijoddinc*2);
K(:,ijout) = Kout(:,2:2:end);

%Ktemp = helmstokesstreamsubmat(zk,s(:,j),t(:,i),sn(:,j),cs,cd);
%indices = 2*(1:length(j)) - mod(j,2);
%K = Ktemp(:,indices);
wj = swhts(j);
K = bsxfun(@times,K,wj(:).');

nrow = size(t,2);
lininds = bsxfun(@plus,i(:),nrow*(j(:)-1).');
lininds = lininds(:);
streamindsi= streaminds(i,:);
streamlinindsi = bsxfun(@plus,(streamindsi-1)*nrow,i(:));
streamlinindsi = streamlinindsi(:);
streamtdi = streamtd(i,:);

[~,ik,itdi] = intersect(lininds,streamlinindsi);

K(ik) = streamtdi(itdi);


%for ii = 1:length(i)
%    [~,ij,id] = intersect(j(:),streaminds(i(ii),:));
%    K(ii,ij) = streamtd(i(ii),id);
%end

end