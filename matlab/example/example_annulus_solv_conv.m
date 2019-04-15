
% EXAMPLE_ANNULUS_EIG_CONV
% Solve an integral equation with rhs 1, and compute the vorticity
% at an interior target for zk=13.7 using both 
% I-2D-2iS and I-2D and do a convergence study
% The geometry is an annulus of size 1<r<1.7;

fileout = ['mat-files/annulus_solv_conv.mat'];

zk = 13.7; % Oscillatory stokes parameter

% annulus params
r1 = 1.0;
r2 = 1.7;

addpaths_loc();
nnchi = 10;
nchi_start = 4;

chunkers = cell(2,nnchi);

vals_comb = zeros(nnchi,1)+1i*zeros(nnchi,1);
vals_dl = zeros(nnchi,1)+1i*zeros(nnchi,1);


errs_comb = zeros(nnchi,1);
sigma = cell(2,nnchi);
rhs = cell(nnchi,1);
errs_dl = zeros(nnchi,1);
npts = zeros(nnchi,1);

% Set options for FLAM
opts2 = []; opts2.FLAM = 1; opts2.verb = true;

r = 1.2 +rand*0.3;
thet = rand*2*pi;
targs = zeros(2,1);
targs(1) = r*cos(thet); targs(2)=r*sin(thet);

opts=[];
opts.usesmooth=2;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-13,'AbsTol',1.0e-13};
opts.gausseps = 1e-13;
intparams = []; intparams.intorder = 16;


for i=1:nnchi
    % setup geometry
    nchi = nchi_start+i-1;
    ncho = ceil(nchi/r1*r2)+1;

    chunkert = circle_chunks(ncho,r2);
    chunkers{1,i} = chunkert;


    chunkert = circle_chunks(nchi,r1);
    chunkert = chunkreverse(chunkert);
    chunkert = chunksort(chunkert);
    chunkers{2,i} = chunkert;
        
    [chunker,nchs] = chunkermerge({chunkers{:,i}});
    
    
    npts(i) = (nchi+ncho)*16;
    n = npts(i)*2;
    rhs{i} = zeros(n,1) + 1i*zeros(n,1);
    n1 = ncho*16*2;
    rhs{i}(1:2:n1) = 1;
    rhs{i}(n1+2:2:n) = 1;
    

    % Compute determinant using combined field reprensetation

    cd = -2.0 + 1i*0.0;
    cs = 0.0 - 1i*2.0;
    
    [d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts2);
    sigma{1,i} = rskelf_sv(F,rhs{i});
    
    vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
    ndimsv = [1 2];
    
    vals_comb(i) = chunkerintkern(chunker,vortkern,ndimsv,sigma{1,i},targs,opts); 
    
    cd = -2.0 + 1i*0.0;
    cs = 0.0;
    
    [d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts2);
    sigma{2,i} = rskelf_sv(F,rhs{i});
    
    vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
    ndimsv = [1 2];
    
    vals_dl(i) = chunkerintkern(chunker,vortkern,ndimsv,sigma{2,i},targs,opts);
    
end

%% Compute exact solution
% setup geometry
nchi = 20;
ncho = ceil(nchi/r1*r2)+1;

chunkers_ex = {};
chunkert = circle_chunks(ncho,r2);
chunkers_ex{1} = chunkert;


chunkert = circle_chunks(nchi,r1);
chunkert = chunkreverse(chunkert);
chunkert = chunksort(chunkert);
chunkers_ex{2} = chunkert;

[chunker,nchs] = chunkermerge(chunkers_ex);


nn = (nchi+ncho)*16;
n = nn*2;
sigma_ex = zeros(n,1) + 1i*zeros(n,1);
rhs_ex = zeros(n,1) + 1i*zeros(n,1);
n1 = ncho*16*2;
rhs_ex(1:2:n1) = 1;
rhs_ex(n1+2:2:n) = 1;

sigma_ex = cell(2,1);


% Compute determinant using combined field reprensetation

cd = -2.0 + 1i*0.0;
cs = 0.0 - 1i*2.0;

[d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts2);
sigma_ex{1} = rskelf_sv(F,rhs_ex);

vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
ndimsv = [1 2];

vals_comb_true = chunkerintkern(chunker,vortkern,ndimsv,sigma_ex{1},targs,opts); 

cd = -2.0 + 1i*0.0;
cs = 0.0;

[d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts2);
sigma_ex{2} = rskelf_sv(F,rhs_ex);

vortkern = @(s,t,sn,tn) ostokesvortkern(zk,cs,cd,s,t,sn,tn);
ndimsv = [1 2];

vals_dl_true = chunkerintkern(chunker,vortkern,ndimsv,sigma_ex{2},targs,opts);

%% Compute the errors
errs_comb = abs(vals_comb-vals_comb_true)/abs(vals_comb_true);
errs_dl = abs(vals_dl-vals_dl_true)/abs(vals_dl_true);

figure(1)
loglog(npts,errs_comb,'k.','MarkerSize',20), hold on
loglog(npts,(npts/150).^-20,'k--')
grid on
xlim([150,650])
ylim([10^-13,10^-1])
xticks([150,300,450,600])
saveas(gcf,'res-files/solv-conv-comb.pdf')

figure(2)
loglog(npts,errs_dl,'k.','MarkerSize',20), hold on
loglog(npts,(npts/150).^-20,'k--')
grid on
xlim([150,650])
ylim([10^-13,10^-1])
xticks([150,300,450,600])
saveas(gcf,'res-files/solv-conv-dl.pdf')

save(fileout,'chunkers','chunkers_ex','errs_comb','errs_dl','npts','nchi',...
      'opts','opts2','targs','vals_comb','vals_comb_true','zk','vals_dl_true',...
      'rhs','rhs_ex','sigma','sigma_ex','r1','r2')
