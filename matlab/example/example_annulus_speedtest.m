
%EXAMPLE_ANNULUS_SPEEDTEST
%
% record some timing information as the frequency and number 
% of discretization points is varied


% annulus params

r1 = 1.0;
r2 = 1.7;



nzkvals = 10;

zkvals = rand(1,nzkvals) + (1:1:nzkvals);
zkvals = zkvals*10;

nchvals = 5:5:100;


tt = zeros(length(zkvals),length(nchvals));
detvals = zeros(length(zkvals),length(nchvals));

seed = 8675309;
rng(seed);
% addpath('../src')
% addpath('../example')
% addpath('../../mwrap')

addpaths_loc();

for ii=1:length(nchvals)

    chunkers = {};
    nchi = nchvals(ii);
    ncho = ceil(nchi/r1*r2)+1;

    chunkert = circle_chunks(ncho,r2);
    checkchunker(chunkert);
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


%%

% build matrix

    cd = -2.0 + 1i*0.0;
    cs = 0.0 - 1i*2.0;
    ncomp = length(nchs);

%

    opts = []; opts.FLAM = 1; opts.verb = true;
    detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
        cd,opts);   
    for jj=1:length(zkvals)
        aa=tic; detvals(jj,ii)=abs(detfun(zkvals(jj)));
        tt(jj,ii) = toc(aa);
    end
    
end

save('example_annulus_speedtest','tt','zkvals','nchvals','detvals')