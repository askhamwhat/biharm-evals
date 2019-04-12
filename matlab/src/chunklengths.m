function chunklens = chunklengths(chunker)

[~,w] = legeexps(chunker.k);

ds = reshape(sum((chunker.ders).^2,1),chunker.k,chunker.nch);
dt = reshape(kron(chunker.hs,w),chunker.k,chunker.nch);

chunklens = reshape(sum(ds.*dt,1),chunker.nch,1);
