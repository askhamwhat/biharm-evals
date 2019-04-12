function ier = checkadjinfo(chunker)

i1 = 1;

hit = zeros(chunker.nch,1);

for i = 1:chunker.nch
    i2 = chunker.adjs(1,i1);
    hit(i2) = hit(i2)+1;
    i1 = i2;
end

ier = 0;
if nnz(hit == 0) > 0
    ier = ier + 1;
end
if nnz(hit > 1) > 0
    ier = ier + 2;
end
% fprintf('nch %d\n',chunker.nch)
% fprintf('number hit %d\n',nnz(hit > 0))
% fprintf('number missed %d\n',nnz(hit == 0))
% fprintf('number doubled up %d\n',nnz(hit > 1))