
function [Kpxy] = ostokes_pxy_kern(s,t,sn,tn,sw,tw,slf,zk,cs,cd)
%OSTOKES_PXY_KERN this function evaluates the interactions between
% the given source points and a proxy surface for the oscillatory 
% Stokes equation (single and double layer potentials)
%
% on input: 
%    s - full list of sources (2,ns)
%    t - proxy targets (2,nt)
%    sn - full list of source normals (2,ns)
%    tn - proxy target normals (2,nt)
%    sw - smooth integration weights for the source
%    tw - smooth integration weights for the targ
%    slf - relevant sources in s,sn,sw arrays
%    zk - oscillatory stokes parameter
%    cs - coefficient of single layer potential
%    cd - coefficient of double layer potential
%
% on output:
%    Kpxy - matrix of interactions

Ktemp = helmstokessubmat(zk,s(:,slf),t,sn(:,slf),cs,cd);

% index trickery to grab appropriate part of source 
% (the sources for a given point correspond to the 1st
% or second component of the density for odd and even 
% indices, respectively. The submatrix routine always 
% computes the contribution for both)

indcs = 2*(1:length(slf))-mod(slf,2);

% scale columns by int weights and add normal integration
% ones-type-matrix for appropriate indices
wj = sw(slf(:));
Kpxy = bsxfun(@times,Ktemp(:,indcs),(wj(:)).');
rnslf = sn(:,slf); 
rn1slf = rnslf(2*(1:length(slf))-mod(slf,2));
rw = rn1slf(:).*wj;
onefix = bsxfun(@times,tn(:),rw(:).');
Kpxy = Kpxy + onefix;

% scale columns by int weights and add normal integration
% ones-type-matrix for appropriate indices
Ktemp = helmstokessubmat(zk,t,s(:,slf),tn,cs,cd);
tw2 = [(tw(:)).';(tw(:)).'];
% scale by int whts
Ktemp = bsxfun(@times,Ktemp(indcs,:),(tw2(:)).');
% add ones-mat-like effect
rw = tn(:).*tw2(:);
onefix = bsxfun(@times,rn1slf(:),rw(:).');
Ktemp = Ktemp + onefix;

Kpxy = [Kpxy; Ktemp.'];


end

