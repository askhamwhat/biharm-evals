
function yylap = hbhdirmat_lap_pre(chunker,nchs,intparams)

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

% want whtsbycomp.'*S*(1/2+sprime+ones)^(-1)*Stau, use transpose trick

% evaluate S^T*whtsbycomp 

fkernS = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');
ndims(1) = 1; ndims(2) = 1;
tic; slap = chunkskernelmat(chunker,fkernS,ndims,intparams); toc
%tic; [Std,Sinds] = chunkskernelmattd(chunker,fkernS,ndims,intparams); toc

slapintstrans = slap.'*whtsbycomp;

% then solve (1/2+sprime+ones)^(-T)*(S^T*whtsbycomp)

fkernSprime = @(s,t,sn,tn) glapkern(s,t,sn,tn,'sprime');
ndims(1) = 1; ndims(2) = 1;
tic; sprime = chunkskernelmat(chunker,fkernSprime,ndims,intparams); toc

sprime = sprime + 0.5*eye(ngeo/2,ngeo/2) + chunkonesmat(chunker);

afun = @(x) sprime.'*x;

yylap = zeros(size(slapintstrans));

for i = 1:ncomp
    yylap(:,i) = gmres(afun,slapintstrans(:,i),[],1.0e-12,50);
end

% finally evaluate stau^T*((1/2+sprime+ones)^(-T)*(S^T*whtsbycomp))

ndim1 = 1;
stau = zeros(size(slap));
for i = 1:ngeo/2
    stau(:,i) = chunkderf(chunker,slap(:,i),ndim1);
end

yylap = stau.'*yylap;
