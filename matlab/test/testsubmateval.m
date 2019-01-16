%TESTSUBMATEVAL


ns = 16;
nt = 16;

src = randn(2,ns);
targ = randn(2,nt);
src_norm = randn(2,ns);

cs = 1.0;
cd = 1.0;
zk = 1.0;

tic; mat = helmstokessubmat(zk,src,targ,src_norm,cs,cd); time2 = toc

ns*nt/time2
