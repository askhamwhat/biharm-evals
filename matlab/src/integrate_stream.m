
function yyis = integrate_stream(chunker,zk,cs,cd,intparams,whtsbycomp,opt)
  %% compute the integral of the stream function
  % matrix along each component

if nargin < 7
    opt = [];
end
if ~isfield(opt,'fast')
    opt.fast = 1;
end
if ~isfield(opt,'fstream_tol')
    opt.fstream_tol = 1.0e-14;
end
if ~isfield(opt,'fstream_np')
    opt.fstream_np = 64;
end
if ~isfield(opt,'fstream_occ')
    opt.fstream_occ = 200;
end
if ~isfield(opt,'verb')
    opt.verb = 0;
end

tol = opt.fstream_tol;
occ = opt.fstream_occ;
np = opt.fstream_np;

% geometry info (doubling for 2d kernel parts)

ngeo = chunker.nch*chunker.k;
whts = chunkwhts(chunker); whts = whts(:);
whts2 = repmat(whts.',2,1); whts2 = reshape(whts2,2*ngeo,1);
rnorms = chunknormals(chunker); rnorms = reshape(rnorms,2,ngeo);
rnorms2 = repmat(rnorms,2,1); rnorms2 = reshape(rnorms2,2,2*ngeo);
x = chunker.chunks; x = reshape(x,2,ngeo);
x2 = repmat(x,2,1); x2 = reshape(x2,2,2*ngeo);

% kernel settings

fkernstream = @(s,t,sn,tn) helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
ndims(1) = 1; ndims(2) = 2;

% proxy settings

[p,pn] = proxy_square_pts(np);

pkernstreamr = @(s,t,sn,tn,slf) pxykernstreamr(s,t,sn,tn,slf,zk,cs,cd,...
    whts) ;
pkernstreamc = @(s,t,sn,tn,slf) pxykernstreamc(s,t,sn,tn,slf,zk,cs,cd,...
    whts2) ;
pxyfunstream = @(rc,rx,cx,slf,nbr,l,ctr) pxyfun_sqr(rc,rx,cx,slf,nbr, ...
    l,ctr,pkernstreamr,pkernstreamc,rnorms,rnorms2,p,pn);

% compute special quadrature on tridiagonal part

start = tic; [streamtd,streaminds] = chunkskernelmattd(chunker, ...
    fkernstream,ndims,intparams); t = toc(start);

if opt.verb; fprintf('%5.2e : time for tridiag comp\n',t); end

% use flam to get fast mv

matfunstream = @(i,j) stream_matfun_wd(x2,x,rnorms2,[],i,j,zk,cs,cd, ...
    whts2,streamtd,streaminds);

optsfmm.store = 'A';
optsfmm.near = 0;
start = tic; F_fmm = ifmm(matfunstream,x,x2,occ,tol,pxyfunstream,optsfmm); 
t=toc(start);

if opt.verb; fprintf('%5.2e : time for ifmm\n',t); end

yyis = (ifmm_mv(F_fmm,whtsbycomp,[],'T')).';

end

function K = pxykernstreamc(s,t,sn,tn,slf,zk,cs,cd,swhts) 
    K = helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
    indcs = 2*(1:length(slf))-mod(slf,2);
    K = K(:,indcs)*diag(swhts(slf));
end

function K = pxykernstreamr(s,t,sn,tn,slf,zk,cs,cd,swhts) 
    K = helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
end

function K = stream_matfun_wd(s,t,sn,tn,i,j,zk,cs,cd,swhts,streamtd,...
    streaminds)
%% evaluate stream matrix by entry index utility function
% 
% this function primarily calls the FORTRAN based submat evaluator
% but makes some attempts at efficiency%
%
% in particular, this first part checks for column pairs which come
% from the same source point (j(ii) = j(ii)+1 for such indices)

ijodd = find(mod(j,2) == 1); jodd = j(ijodd);
[jeveninc,ijeveninc,ijoddinc] = intersect(j(:),jodd(:)+1); 
[jevenout,ijout] = setdiff(j,[jodd(:);jeveninc(:)]);

Kodd= helmstokesstreamsubmat(zk,s(:,jodd),t(:,i),sn(:,jodd),cs,cd);
Kout = helmstokesstreamsubmat(zk,s(:,jevenout),t(:,i),sn(:,jevenout),cs,cd);
K = zeros(length(i),length(j));
K(:,ijodd) = Kodd(:,1:2:end);
K(:,ijeveninc) = Kodd(:,ijoddinc*2);
K(:,ijout) = Kout(:,2:2:end);

wj = swhts(j);
K = bsxfun(@times,K,wj(:).');

% relatively efficient replacement of entries which require special
% quadrature (such entries should have been precomputed)

nrow = size(t,2);
lininds = bsxfun(@plus,i(:),nrow*(j(:)-1).');
lininds = lininds(:);
streamindsi= streaminds(i,:);
streamlinindsi = bsxfun(@plus,(streamindsi-1)*nrow,i(:));
streamlinindsi = streamlinindsi(:);
streamtdi = streamtd(i,:);

[~,ik,itdi] = intersect(lininds,streamlinindsi);

K(ik) = streamtdi(itdi);


end
