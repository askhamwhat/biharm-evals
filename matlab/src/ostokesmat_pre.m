
function [stokestd,stokesinds] = ...
    ostokesmat_pre(chunker,nchs,zk,cs,cd,intparams)
%STOKESMAT_PRE
%
% This function computes the tridiagonal part of the oscillatory Stokes
% system matrix for the Dirichlet problem on the given domain.
% 

% geometry stuff

k = chunker.k;
nch = chunker.nch;

% form stokes densities to velocity field (standard stokes dirichlet mat)
% near and self only

fkern = @(s,t,sn,tn) helmstokessubmat(zk,s,t,sn,cs,cd);
ndims(1) = 2; ndims(2) = 2;
[stokestd,stokesinds] = chunkskernelmattd(chunker,fkern,ndims,intparams);

stokestd(:,1+ndims(2)*k:ndims(2)*2*k) = ...
    stokestd(:,1+ndims(2)*k:ndims(2)*2*k) ...
    - 0.5*cd*kron(ones(nch,1),eye(ndims(2)*k,ndims(2)*k));

end
