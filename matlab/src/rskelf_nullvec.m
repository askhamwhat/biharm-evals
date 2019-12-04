
function xnull = rskelf_nullvec(F,nsys,p,q)
%RSKELF_NULLVEC find nullvec with cheap randomized method
%
% (F + u*v.')x = u
%


xnull = rskelf_sv(F,randn(nsys,p));
[xnull,~,~] = qr(xnull,0);
for i = 1:q
   xnull = rskelf_sv(F,xnull,'C');
   [xnull,~,~] = qr(xnull,0);
   xnull = rskelf_sv(F,xnull);
   [xnull,~,~] = qr(xnull,0);
end

end


function y = rskelf_plus_mv(x,F,u,v)

y = rskelf_mv(F,x);
y = y+u*v.'*x;

end