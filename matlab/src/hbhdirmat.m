
function sysmat = hbhdirmat(chunker,nchs,ccs,zk,cs,cd,intparams)
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

sysmat = zeros(nsys,nsys) + 1i*zeros(nsys,nsys);


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

sysmat(ngeo+1:end,ngeo+1:end) = whtsbycomp.'*logs;

% top right is velocity due to log charges (grad perp of charges)

sysmat(1:ngeo,ngeo+1:end) = logperps;

%% effect of stokeslet-type layer potentials

% form stokes densities to velocity field (standard stokes dirichlet mat)

fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
ndims(1) = 2; ndims(2) = 2;
sysmat(1:ngeo,1:ngeo) = chunkskernelmat(chunker,fkern,ndims,intparams);
normonesmat = chunknormonesmat(chunker);
sysmat(1:ngeo,1:ngeo) = -0.5*cd*eye(ngeo,ngeo) ...
    + sysmat(1:ngeo,1:ngeo) + normonesmat;

% the tricky part: evaluate corresponding stream function on
% boundary components and integrate

% it has two parts: one is straightforward, just an evaluation
% the other is not given as a stream function, so the corresponding
% stream must be computed via Dirichlet to Neumann map...

% want whtsbycomp.'*S*(1/2+sprime+ones)^(-1)*Stau, use transpose trick

% evaluate S^T*whtsbycomp 

fkernS = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');
ndims(1) = 1; ndims(2) = 1;
slap = chunkskernelmat(chunker,fkernS,ndims,intparams);

slapintstrans = slap.'*whtsbycomp;

% then solve (1/2+sprime+ones)^(-T)*(S^T*whtsbycomp)

fkernSprime = @(s,t,sn,tn) glapkern(s,t,sn,tn,'sprime');
ndims(1) = 1; ndims(2) = 1;
sprime = chunkskernelmat(chunker,fkernSprime,ndims,intparams);

sprime = sprime + 0.5*eye(ngeo/2,ngeo/2) + chunkonesmat(chunker);

afun = @(x) sprime.'*x;

yy = zeros(size(slapintstrans));

for i = 1:ncomp
    yy(:,i) = gmres(afun,slapintstrans(:,i),[],1.0e-12,50);
end

% finally evaluate stau^T*((1/2+sprime+ones)^(-T)*(S^T*whtsbycomp))

ndim1 = 1;
stau = zeros(size(slap));
for i = 1:ngeo/2
    stau(:,i) = chunkderf(chunker,slap(:,i),ndim1);
end

yy = stau.'*yy;

% stream function part

fkernstream = @(s,t,sn,tn) helmstokesstreamsubmat(zk,s,t,sn,cs,cd);
ndims(1) = 1; ndims(2) = 2;
streammat = chunkskernelmat(chunker,fkernstream,ndims,intparams);

% integrate stream part

yy2 = whtsbycomp.'*streammat;

sysmat(ngeo+1:end,1:ngeo) = -cd*kron(yy.',ones(1,2)).*(repmat(rnorms(:).',ncomp)) ...
    +yy2;




