
function [stokestd,stokesinds,sysmatbl,sysmatbr,sysmattr] = hbhdirmat_pre(chunker,nchs,ccs,zk,cs,cd,intparams,varargin)
%HBHDIRMAT
%
% This function builds the system matrix for the Helmholtz biharmonic
% Dirichlet problem on the given domain. This is the 'slow' version...
%



% geometry stuff

k = chunker.k;
nch = chunker.nch;

whts = chunkwhts(chunker); whts = whts(:);
rnorms = chunknormals(chunker); rnorms = reshape(rnorms,2,k*nch);

ncomp = length(nchs);

whtsbycomp = zeros(length(whts),ncomp);

ind = 1;
for i = 1:ncomp
    whtsbycomp(ind:ind+nchs(i)*k-1,i) = whts(ind:ind+nchs(i)*k-1);
    ind = ind+nchs(i)*k;
end

ngeo = 2*k*nch;
nsys = 2*k*nch+ncomp;

%sysmat = zeros(nsys,nsys) + 1i*zeros(nsys,nsys);


%% effect of sources in holes

chunks = reshape(chunker.chunks,2,k*nch);
logs = zeros(k*nch,ncomp);
logperps = zeros(2*k*nch,ncomp);

% outer boundary 'source' is just a constant

logs(:,1) = 1;

% get other sources if they exist

if ncomp > 1
    [val,grad] = glapfun(ccs(:,2:end),chunks);
    logs(:,2:end) = val;
    logperps(1:2:end,2:end) = -grad(:,:,2);
    logperps(2:2:end,2:end) = grad(:,:,1);
end

% bottom right is integral of log charges

sysmatbr = whtsbycomp.'*logs + 1i*0.0;

% top right is velocity due to log charges (grad perp of charges)

sysmattr = logperps + 1i*0.0;

%% effect of stokeslet-type layer potentials

% form stokes densities to velocity field (standard stokes dirichlet mat)
% near and self only

fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
ndims(1) = 2; ndims(2) = 2;
[stokestd,stokesinds] = chunkskernelmattd(chunker,fkern,ndims,intparams);

stokestd(:,1+ndims(2)*k:ndims(2)*2*k) = stokestd(:,1+ndims(2)*k:ndims(2)*2*k) ...
-0.5*cd*kron(ones(nch,1),eye(ndims(2)*k,ndims(2)*k));

% the tricky part: evaluate corresponding stream function on
% boundary components and integrate

% it has two parts: one is straightforward, just an evaluation
% the other is not given as a stream function, so the corresponding
% stream must be computed via Dirichlet to Neumann map...

% this is now performed by the hbhdirmat_lap_pre routine...

if nargin > 7
    yy = varargin{1};
    size(yy)
else
    yy = hbhdirmat_lap_pre(chunker,nchs,intparams);
end

% stream function part

fkernstream = @(s,t,sn,tn) helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
ndims(1) = 1; ndims(2) = 2;
streammat = chunkskernelmat(chunker,fkernstream,ndims,intparams);
size(streammat)
% integrate stream part

yy2 = whtsbycomp.'*streammat;
size(repmat(rnorms(:).',ncomp))
size(yy2)
size(kron(yy.',ones(1,2)))
size(-cd*kron(yy.',ones(1,2)).*(repmat(rnorms(:).',ncomp)) ...
    +yy2);

sysmatbl = -cd*kron(yy.',ones(1,2)).*(repmat(rnorms(:).',ncomp)) + yy2;




