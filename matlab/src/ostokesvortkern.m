function mat = ostokesvortkern(zk,cs,cd,src,targ,srcn,~)

[~,nt] = size(targ);
[~,ns] = size(src);

srct = zeros(size(srcn));
srct(1,:) = -srcn(2,:);
srct(2,:) = srcn(1,:);

[~,g,h] = helmfun(zk,src,targ);

srcnxmat = repmat(srcn(1,:),nt,1);
srcnymat = repmat(srcn(2,:),nt,1);
srctxmat = repmat(srct(1,:),nt,1);
srctymat = repmat(srct(2,:),nt,1);

g2 = g(:,:,2); g1 = g(:,:,1);
spart = [g2; -g1];

spart = reshape(spart,nt,2*ns);

dpartn = 2.0*(srcnxmat.*srctxmat.*h(:,:,1) ...
    + (srcnxmat.*srctymat+srcnymat.*srctxmat).*h(:,:,2) ...
    + srcnymat.*srctymat.*h(:,:,3));

dpartt = ( (srctxmat.*srctxmat-srcnxmat.*srcnxmat).*h(:,:,1) ...
    + 2.0*(srctxmat.*srctymat-srcnxmat.*srcnymat).*h(:,:,2) ...
    + (srctymat.*srctymat - srcnymat.*srcnymat).*h(:,:,3));

dpartn2 = repmat(dpartn,2,1); dpartn2 = reshape(dpartn2,nt,2*ns);
dpartt2 = repmat(dpartt,2,1); dpartt2 = reshape(dpartt2,nt,2*ns);

mat = cd*(bsxfun(@times,dpartn2,(srcn(:)).')+...
    bsxfun(@times,dpartt2,(srct(:)).')) + cs*spart;
