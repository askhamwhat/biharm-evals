function mat = ostokeskern(zk,src,targ,srcn,targn,type)

if strcmpi(type,'stresslet')
    cd = 1.0 + 1i*0.0; cs = 0.0 + 1i*0.0;
    mat = helmstokessubmat(zk,targ,src,targn,cs,cd);
    mat = mat.';
elseif strcmpi(type,'double')
    cs = 0.0 + 1i*0.0; cd = 1.0 + 1i*0.0;
    mat = helmstokessubmat(zk,src,targ,srcn,cs,cd);
else
    cs = 1.0 + 1i*0.0; cd = 0.0 + 1i*0.0;
    mat = helmstokessubmat(zk,src,targ,srcn,cs,cd);
end
