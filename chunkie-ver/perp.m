function vperp = perp(v)
%CHNK.UTIL.PERP gives perp of 2 vectors. if array input, assumes 
% first dimension is the one to perp.

szv = size(v);

vperp = v;
vperp(1,:) = v(2,:);
vperp(2,:) = -v(1,:);

end
