%TEST_GREENOSTOKES test the routines for integrating over chunks against 
% the Green's id for oscillatory Stokes
%
% 

clear all; clc; clf;
seed = 8675309;
rng(seed);

addpaths_loc();

doadap = false;

% wave number

zk = 0.4;
opdims = [2 2];

% geometry parameters and construction

iseed = 8675309;
rng(iseed);

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 2;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref);
t1 = toc(start);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(opdims(2)*ns,1);

% targets

nt = 1;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

figure(1)
clf
hold off
plot(chnkr)
hold on
quiver(chnkr)
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')
axis equal 


%%

% kernel defs

kernd = @(s,t,stau,ttau) ostokes2dkern(zk,s,t,stau,ttau,'double');
kerns = @(s,t,stau,ttau) ostokes2dkern(zk,s,t,stau,ttau,'single');
kernstress = @(s,t,stau,ttau) ostokes2dkern(zk,s,t,stau,ttau,'stresslet');

% eval u and traction on boundary

targs = chnkr.r; targs = targs(:,:);
targstau = taus(chnkr); targstau = targstau(:,:);

kernmats = kerns(sources,targs,sources,targstau);
kernmatstress = kernstress(sources,targs,sources,targstau);
densu = kernmats*strengths(:);
denstrac = kernmatstress*strengths(:);

%

% eval u at targets

% for ordinary charges, tau vectors don't matter
kernmatstarg = kerns(sources,targets,sources,targets);
utarg = kernmatstarg*strengths(:); utarg = reshape(utarg,2,nt);

% test green's id

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-8,'AbsTol',1.0e-8};
start=tic; Du = chunkerintkern(chnkr,kernd,opdims,densu,targets,opts); 
toc(start)
start=tic; Strac = chunkerintkern(chnkr,kerns,opdims,denstrac,...
    targets,opts); toc(start)

%

utarg2 = Strac-Du; utarg2 = reshape(utarg2,2,nt);

%

norm(utarg-utarg2,'fro')/norm(utarg,'fro')

