
function xnull = rskelf_nullvec(F,nsys,p,q)
%RSKELF_NULLVEC find nullvec with cheap randomized method
%
% (F + u*v.')x = u
%

%ur = randn(nsys,nvec);
%vr = randn(nsys,nvec);

%xnull = zeros(nsys,nvec);
%for i = 1:nvec
%    xnull(:,i) = gmres(@(x) rskelf_plus_mv(x,F,ur(:,i),vr(:,i)), ...
%        ur(:,i),restart,tol,maxit);
%end

xnull = rskelf_sv(F,randn(nsys,p));
[xnull,~,~] = qr(xnull,0);
for i = 1:q
   xnull = rskelf_sv(F,xnull,'T');
   [xnull,~,~] = qr(xnull,0);
   xnull = rskelf_sv(F,xnull);
   [xnull,~,~] = qr(xnull,0);
   
end

end


function y = rskelf_plus_mv(x,F,u,v)

y = rskelf_mv(F,x);
y = y+u*v.'*x;

end