function checkchunker(chunker)

assert(isstruct(chunker),'chunker must be a struct');
assert(isfield(chunker,'k'));
assert(isfield(chunker,'nch'));
assert(isfield(chunker,'chunks'));
assert(isfield(chunker,'ders'));
assert(isfield(chunker,'ders2'));
assert(isfield(chunker,'adjs'));
assert(isfield(chunker,'hs'));

nch = chunker.nch; k = chunker.k;
assert(numel(chunker.chunks)==2*nch*k);
assert(numel(chunker.ders)==2*nch*k);
assert(numel(chunker.ders2)==2*nch*k);
assert(numel(chunker.adjs)==2*nch);
assert(numel(chunker.hs)==nch);
%assert(checkadjinfo(chunker)==0);