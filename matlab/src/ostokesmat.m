function sysmat = ostokesmat(chunker,zk,cs,cd,intparams)

k = chunker.k;
nch = chunker.nch;
ngeo = 2*k*nch;

%% effect of stokeslet-type layer potentials

% form stokes densities to velocity field (standard stokes dirichlet mat)

fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
ndims(1) = 2; ndims(2) = 2;
sysmat(1:ngeo,1:ngeo) = chunkskernelmat(chunker,fkern,ndims,intparams);
normonesmat = chunknormonesmat(chunker);
sysmat(1:ngeo,1:ngeo) = -0.5*cd*eye(ngeo,ngeo) ...
    + sysmat(1:ngeo,1:ngeo) + normonesmat;
