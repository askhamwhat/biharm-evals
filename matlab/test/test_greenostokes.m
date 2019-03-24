%TEST_GREENOSTOKES test the routines for integrating over chunks against 
% the Green's id for oscillatory Stokes
%
% 

clear all; clc; clf;
seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

doadap = false;

% geometry parameters and construction

ndims = [2 2];

zk = 0.4;

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 1;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
fvals = starfish(ts,narms,amp);
sources = zeros(2,ns);
sources(1,:) = fvals(:,1); sources(2,:) = fvals(:,2);
sources = 1.5*sources;
strengths = randn(2,ns,1);

% targets

nt = 100;
ts = 0.0+2*pi*rand(nt,1);
fvals = starfish(ts,narms,amp);
targets = zeros(2,nt);
targets(1,:) = fvals(:,1); targets(2,:) = fvals(:,2);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));

hold off
scatter(xs(:),ys(:))
hold on
scatter(sources(1,:),sources(2,:),'o')
scatter(targets(1,:),targets(2,:),'x')

%

% kernel defs

kernd = @(s,t,sn,tn) ostokeskern(zk,s,t,sn,tn,'double');
kerns = @(s,t,sn,tn) ostokeskern(zk,s,t,s,t,'single');
kernstress = @(s,t,sn,tn) ostokeskern(zk,s,t,s,tn,'stresslet');

% eval u and traction on boundary

targs = chunker.chunks; targs = reshape(targs,2,chunker.k*chunker.nch);
targsn = chunknormals(chunker); 
targsn = reshape(targsn,2,chunker.k*chunker.nch);

kernmats = kerns(sources,targs,sources,targsn);
kernmatstress = kernstress(sources,targs,sources,targsn);
densu = kernmats*strengths(:);
denstrac = kernmatstress*strengths(:);

%%

% eval u at targets

kernmatstarg = kerns(sources,targets,sources,targets);
utarg = kernmatstarg*strengths(:); utarg = reshape(utarg,2,nt);

% test green's id

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1.0e-8,'AbsTol',1.0e-8};
start=tic; Du = chunkerintkern(chunker,kernd,ndims,densu,targets,opts); 
toc(start)
start=tic; Strac = chunkerintkern(chunker,kerns,ndims,denstrac,...
    targets,opts); toc(start)

%

utarg2 = Strac-Du; utarg2 = reshape(utarg2,2,nt);

%

norm(utarg-utarg2,'fro')/norm(utarg,'fro')

