
%PLOT_ANNULUS

addpaths_loc();

filebase = 'example_annulus_001'; 
timeref = datestr(now,'_yyyymmdd_HHMMSS');
fileout = [filebase, '_reptest.mat'];

% annulus params

r1 = 1.0;
r2 = 1.7;


chunkers = {};
nchi = 6;
ncho = 10;

chunkert = circle_chunks(ncho,r2);
chunkers{1} = chunkert;

nw = 1;
nh = 1;

ind = 2;
chunkert = circle_chunks(nchi,r1);
chunkert = chunkreverse(chunkert);
chunkert = chunksort(chunkert);
chunkers{ind} = chunkert;
        
[chunker,nchs] = chunkermerge(chunkers);

%

nn = nchi+ncho;
xx1 = chunker.chunks(1,:,1:2:nn); xx1 = xx1(:);
yy1 = chunker.chunks(2,:,1:2:nn); yy1 = yy1(:);


xx2 = chunker.chunks(1,:,2:2:nn); xx2 = xx2(:);
yy2 = chunker.chunks(2,:,2:2:nn); yy2 = yy2(:);

rnorms = chunknormals(chunker);

intparams.intorder = chunker.k;

plot(xx1,yy1,'r.','MarkerSize',15); hold on;
plot(xx2,yy2,'b.','MarkerSize',15);
saveas(gcf,'annulus_discretization.pdf')