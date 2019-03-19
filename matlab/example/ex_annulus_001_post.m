
% replace with location of chebfun/flam on your computer

addpath('../src','../../mwrap','~/Dropbox/MATLAB/chebfun');
addpath(genpath('~/Dropbox/MATLAB/FLAM'));


%%

load example_annulus_001.mat % I-2D - 2W -2iS

detchebs1 = detchebs;

load example_annulus_002.mat % I-2D-2W

detchebs2 = detchebs;

d45_1 = detchebs1{4};
d56_1 = detchebs1{5};

d45_2 = detchebs2{4};
d56_2 = detchebs2{5};

rts3_2 = roots(detchebs2{3},'complex');
zk = real(rts3_2(3));


%% 

% look at determinants

figure(1)
clf
for i = 1:16
    subplot(4,4,i)
    length(detchebs1{i})
    length(detchebs2{i})
    semilogy(abs(detchebs1{i}),'r')
    hold on
    semilogy(abs(detchebs2{i}),'b')
end

%%

% evaluate determinant at zk using FLAM
% d is determinant
% F is matrix in RSKELF format
opts = []; opts.verb = true; opts.FLAM=1;
[d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,...
    opts);

%%

% find approximate null vector for factorization
nsys = chunker.nch*chunker.k*2;
tol = 1e-5; maxit=40; nref=4;
xnull = rskelf_nullvec(F,nsys,[],tol,maxit,nref);

norm(xnull,'fro')
norm(rskelf_mv(F,xnull),'fro')


%%

% set up grid for plotting
nxdir = 40;
xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
xgrid = linspace(xmin,xmax,nxdir);
ygrid = linspace(ymin,ymax,nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

% find target points which lie inside domain

start = tic; in = chunkerin(chunker,targs); toc(start)

hold off
scatter(targs(1,in),targs(2,in),'go');
hold on
scatter(targs(1,~in),targs(2,~in),'rx');
scatter(xs(:),ys(:),'bo');

%%

% evaluate with adaptive quadrature 

opts.usesmooth=2;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-8,'AbsTol',1.0e-8};
opts.gausseps = 1e-7; % if gauss ID accurate to 7 digits, use smooth
targsin=targs(:,in);
intparams = []; intparams.intorder = 16;

dkern = @(s,t,sn,tn) ostokeskern(zk,s,t,sn,tn,'double');
vkern = @(s,t,sn,tn) ostokesvortkern(zk,0.0,1.0,s,t,sn);
ndimsD = [2 2];
ndimsV = [1 2];
start = tic; 
Din = chunkerintkern(chunker,dkern,ndimsD,xnull,targsin,opts); 
toc(start)
start = tic; 
Vin = chunkerintkern(chunker,vkern,ndimsV,xnull,targsin,opts); 
toc(start)

%%

Din = reshape(Din,[ndimsD(1) nnz(in)]);

D = nan([size(xx) , 2]);
D(in,1) = Din(1,:);
D(in,2) = Din(2,:);

Vin = reshape(Vin,[ndimsV(1) nnz(in)]);

V = nan([size(xx)]);
V(in) = Vin(:);


%%

figure(2)
clf
h = pcolor(xx,yy,real(V));
hold on
quiver(xx(in),yy(in),D(in,1),D(in,2))


