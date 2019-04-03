
%TEST_CORNERS
%
% This file tests whether or not the polygon routine
% in chunks.mw is working

addpath('../src','../../mwrap')

opts = [];
opts.autowidths = false;
opts.autowidthsfac = 0.1;

opts.eps = 1e2;

verts = barbell_example_001(6,3);

opts.widths = 0.3*ones(size(verts,2),1);

chunker = chunkcorner(verts,opts);

chunker.nch

%

nch = chunker.nch;
k = chunker.k;

optsref = [];
optsref.nchmax = 10000;
optsref.maxchunklen = 0.01;

chunker2 = chunkrefine(chunker,optsref);

%%

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

figure(1)
clf
scatter(xs,ys,'rx')
axis equal
