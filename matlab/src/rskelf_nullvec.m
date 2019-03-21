
function xnull = rskelf_nullvec(F,nsys,restart,tol,maxit,nref)
%RSKELF_NULLVEC find nullvec with cheap randomized method
%
% (F + u*v.')x = u
%

ur = randn(nsys,1);
vr = randn(nsys,1);

xnull = gmres(@(x) rskelf_plus_mv(x,F,ur,vr),ur,restart,tol,maxit);

xnull = xnull/norm(xnull);
for i = 1:nref
   xnull = rskelf_sv(F,xnull,'T');
   xnull = rskelf_sv(F,xnull);
   xnull = xnull/norm(xnull);
end

end


function y = rskelf_plus_mv(x,F,u,v)

y = rskelf_mv(F,x);
y = y+u*v.'*x;

end