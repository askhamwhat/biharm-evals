function mat = helmstokesstreamsubmat_wpre(zk,src,targ,src_norm,cs,cd,streammattd,streammatinds)

[~,ns] = size(src);
[~,nt] = size(targ);
mat2 = zeros(1,2,nt,ns);

mex_id_ = 'zhelmstokesallstreammatmany(i dcomplex[x], i int[x], i double[], i int[x], i double[], i dcomplex[], i dcomplex[x], i dcomplex[x], io dcomplex[])';
[mat2] = bhdeiggateway(mex_id_, zk, ns, src, nt, targ, src_norm, cs, cd, mat2, 1, 1, 1, 1, 1);

mat2 = reshape(mat2,1,2,nt,ns);

mat = zeros(1*nt,2*ns);
mat(:,1:2:2*ns) = reshape(mat2(:,1,:,:),nt,ns);
mat(:,2:2:2*ns) = reshape(mat2(:,2,:,:),nt,ns);

end

