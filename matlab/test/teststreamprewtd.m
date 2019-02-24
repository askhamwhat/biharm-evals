%% matrix entry function w/ precomputed tri-diagonal

addpath('../src')
addpath('../../mwrap')

% geometry parameters and construction

cparams.eps = 1.0e-3;
cparams.nchmax = 100000;
cparams.nover = 2;
narms = 5;
amp = 0.5;
tic; chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams); toc

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

%% compute matrix, straightforward way

start = tic; stream1 = chunkskernelmat(chunker,fkernstream,ndims,intparams); toc(start)

%% compute matrix, using indices (smooth rule) with tridiagonal replaced

matfun2 = @(i,j) stream_matfun(x2,x,rnorms2,[],i,j,zk,cs,cd,whts2);
alli = 1:ngeo; allj = 1:2*ngeo;
start = tic; stream2 =  matfun2(alli,allj); toc(start)

% get tridiagonal and replace (singular rule)

start = tic; [streamtd,streaminds] = chunkskernelmattd(chunker,fkernstream,ndims,intparams); toc(start)
lininds = (streaminds-1)*ngeo+(1:ngeo).'; lininds = lininds(:);
stream2(lininds) = streamtd(:);

%% compute matrices using indices (tridiag fixed by routine)

matfun3 = @(i,j) stream_matfun_wd(x2,x,rnorms2,[],i,j,zk,cs,cd,whts2,streamtd,streaminds);
alli = 1:ngeo; allj = 1:2*ngeo;
start = tic; stream3 =  matfun3(alli,allj); toc(start)

%% test

norm(stream1-stream2,'fro')
norm(stream1-stream3,'fro')


function K = stream_matfun(s,t,sn,tn,i,j,zk,cs,cd,swhts)

K= helmstokesstreamsubmat(zk,s(:,j),t(:,i),sn(:,j),cs,cd);
indcs = 2*(1:length(j))-mod(j,2);

K= K(:,indcs)*diag(swhts(j));

end

function K = stream_matfun_wd(s,t,sn,tn,i,j,zk,cs,cd,swhts,streamtd,streaminds)

K= helmstokesstreamsubmat(zk,s(:,j),t(:,i),sn(:,j),cs,cd);
indcs = 2*(1:length(j))-mod(j,2);

K= K(:,indcs)*diag(swhts(j));

for ii = 1:length(i)
    [~,ij,id] = intersect(j(:),streaminds(i(ii),:).');
    K(ii,ij) = streamtd(i(ii),id);
end

end