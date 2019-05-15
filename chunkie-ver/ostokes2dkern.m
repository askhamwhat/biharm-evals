function mat = ostokes2dkern(zk,src,targ,srctau,targtau,type,addw)

if nargin < 7
   addw = false;
end

if strcmpi(type,'stresslet')
   targn = perp(targtau);
    cd = 1.0 + 1i*0.0; cs = 0.0 + 1i*0.0;
    mat = ostokes2dsubmat(zk,targ,src,targn,cs,cd);
    mat = mat.';
elseif strcmpi(type,'double')
    srcn = perp(srctau);
    cs = 0.0 + 1i*0.0; cd = 1.0 + 1i*0.0;
    mat = ostokes2dsubmat(zk,src,targ,srcn,cs,cd);
else
    srcn = perp(srctau);
    cs = 1.0 + 1i*0.0; cd = 0.0 + 1i*0.0;
    mat = ostokes2dsubmat(zk,src,targ,srcn,cs,cd);
end

if addw
   targn = perp(targtau);
   srcn = perp(srctau);
   mat = mat + targn(:)*(srcn(:).');
end
