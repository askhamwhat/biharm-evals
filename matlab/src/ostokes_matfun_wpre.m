function K = ostokes_matfun_wpre(s,sn,sw,i,j,zk,cs,cd,k,nchs, ...
    ncomp,stokestd,stokesinds)
%% evaluate system matrix by entry index utility function for 
% combined field oscillatory stokes integral equation
% 
% this function primarily calls the FORTRAN based submat evaluator
% but makes some attempts at efficiency%
%
% see also OSTOKESMAT_PRE
nch = sum(nchs);
ngeo = nch*k;


% this first part checks for column pairs in stokes part which come
% from the same source point (j(ii) = j(ii)+1) for such indices)

ijodd = mod(j,2) == 1; jodd = j(ijodd);
[jeveninc,ijeveninc,ijoddinc] = intersect(j(:),jodd(:)+1); 
[jevenout,ijout] = setdiff(j,[jodd(:);jeveninc(:)]);

% ditto for row pairs

iiodd = mod(i,2) == 1; iodd = i(iiodd);
[ieveninc,iieveninc,iioddinc] = intersect(i(:),iodd(:)+1); 
[ievenout,iiout] = setdiff(i,[iodd(:);ieveninc(:)]);

Kroddcodd = helmstokessubmat(zk,s(:,jodd),s(:,iodd),sn(:,jodd), ...
    cs,cd);
Kroddcout = helmstokessubmat(zk,s(:,jevenout),s(:,iodd), ...
    sn(:,jevenout),cs,cd);
Kroutcodd = helmstokessubmat(zk,s(:,jodd),s(:,ievenout), ...
    sn(:,jodd), cs,cd);
Kroutcout= helmstokessubmat(zk,s(:,jevenout),s(:,ievenout), ...
    sn(:,jevenout), cs,cd);

K = zeros(length(i),length(j));

% grab everything from odd - odd matrix (have to account for the even
% pairs of these indices, if present)

K(iiodd,ijodd) = Kroddcodd(1:2:end,1:2:end);
K(iieveninc,ijodd) = Kroddcodd(2*iioddinc,1:2:end);
K(iiodd,ijeveninc) = Kroddcodd(1:2:end,2*ijoddinc);
K(iieveninc,ijeveninc) = Kroddcodd(2*iioddinc,2*ijoddinc);

% others are a little simpler

K(iiodd,ijout) = Kroddcout(1:2:end,2:2:end);
K(iieveninc,ijout) = Kroddcout(2*iioddinc,2:2:end);

K(iiout,ijodd) = Kroutcodd(2:2:end,1:2:end);
K(iiout,ijeveninc) = Kroutcodd(2:2:end,2*ijoddinc);

K(iiout,ijout) = Kroutcout(2:2:end,2:2:end);

wj = sw(j);
K = bsxfun(@times,K,wj(:).');

% relatively efficient replacement of entries which require special
% quadrature (such entries should have been precomputed, see
% ostokesmat_pre.m)

lininds = bsxfun(@plus,i(:),2*ngeo*(j(:)-1).'); 
lininds = lininds(:);

i = i(:);
stokesindsi= stokesinds(i,:); 
stokeslinindsi = bsxfun(@plus,i(:),(stokesindsi-1)*2*ngeo);
stokeslinindsi = stokeslinindsi(:);
stokestdi = stokestd(i,:);

[~,ik,itdi] = intersect(lininds,stokeslinindsi);

K(ik) = stokestdi(itdi);

% finally add onemat correction for stokes type entries

ri = sn(:,i); ri=ri(2*reshape(1:length(i),size(i))-mod(i,2));
rj = sn(:,j); rj=rj(2*reshape(1:length(j),size(j))-mod(j,2));
onefix = bsxfun(@times,ri(:),rj(:).');
wj = sw(j);
onefix = bsxfun(@times,onefix,(wj(:)).');
K = K+onefix;

end
