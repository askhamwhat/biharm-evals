function [d,varargout] = ostokes_determinant(zk,chunker,nchs,cs,cd,...
    opts)
%OSTOKES_DETERMINANT
%
% A wrapper which computes the determinant of the specified 
% integral operator on the domain defined by chunker and nchs
%
% On input:
%
% zk - the frequency to check
% chunker - chunks definition of domain
% nchs - number of chunks on each part of boundary (first nchs(1) is outer)
% cs - coefficient of single layer
% cd - coefficient of double layer
% opts - options structure, available options (default):
%       opts.FLAM = 1 use FLAM routines, 0 use dense built-ins (0)
%       opts.FLAMtol = tolerance for compression (1e-14)
%       opts.FLAMocc = occupancy parametr for FLAM tree (200)
%       opts.FLAMpxy = use proxy surface for FLAM compress (1)
%       opts.FLAMpxyorder = number of points on proxy surface, should 
%                be positive multiple of 4 for square surface (64)
%       opts.FLAMpxyshape = shape of proxy surface 'circle' or ('square')
%       opts.intparams = integration parameters structure (
%       intparams.intorder = chunker.k)
%       opts.verb = print out some basic timing info (false)
% 
% On output:
%
% d - determinant of the integral operator 
%        cd(-0.5+D) + cs S
%     where D and S denote the oscillatory stokes double/single layer
% varargout{1} = rskelf factorization if opts.FLAM = 1
%                dense matrix if opts.FLAM = 0
%

if (nargin < 6)
    opts = [];
end

if (~isfield(opts,'intparams'))
    intparams = [];
    intparams.intorder = chunker.k;
    opts.intparams = intparams;
end
if (~isfield(opts,'FLAM'))
    opts.FLAM = 0;
end
if (~isfield(opts,'FLAMtol'))
    opts.FLAMtol = 1.0e-14;
end
if (~isfield(opts,'FLAMocc'))
    opts.FLAMocc = 200;
end
if (~isfield(opts,'FLAMpxy'))
    opts.FLAMpxy = 1;
end
if (~isfield(opts,'FLAMpxyorder'))
    opts.FLAMpxyorder = 64;
end
if (~isfield(opts,'FLAMpxyshape'))
    opts.FLAMpxyshape = 'square';
end
if (~isfield(opts,'verb'))
    opts.verb = false;
end

assert( opts.FLAMocc > 0 , 'opts.FLAMocc should be positive int')
assert( strcmpi(opts.FLAMpxyshape,'square') || ...
    strcmpi(opts.FLAMpxyshape,'circle'), 'opts.FLAMpxyshape invalid')
if (strcmpi(opts.FLAMpxyshape,'square'))
    assert( mod(opts.FLAMpxyorder,4)==0 && opts.FLAMpxyorder > 0, ...
        'opts.FLAMpxyorder should be a positive multiple of 4')
end

nch = chunker.nch; k = chunker.k; npts = nch*k; ncomp = length(nchs(:));

assert(sum(nchs(:))==nch,'total of nchs is inconsistent with chunker.nch')

verb = opts.verb;
cfprintf = @(msg,varargin) cfprintf1(verb,msg,varargin{:});

if opts.FLAM == 1
    
    p = opts.FLAMpxyorder;
    
    if strcmpi(opts.FLAMpxyshape,'square')
        [proxy,pnorm,pw] = proxy_square_pts(p);
    else
        [proxy,pnorm,pw] = proxy_circ_pts(p);
    end
    rnorms = chunknormals(chunker);
    
    xflam = zeros(2,2*npts);
    xflam(:,1:2:1+2*(npts-1)) = chunker.chunks(1:2,:);
    xflam(:,2:2:2+2*(npts-1)) = chunker.chunks(1:2,:);

    xnorm = zeros(2,2*npts);
    xnorm(:,1:2:1+2*(npts-1)) = rnorms(1:2,:);
    xnorm(:,2:2:2+2*(npts-1)) = rnorms(1:2,:);

    whts = chunkwhts(chunker);
    whtsflam = repmat((whts(:)).',2,1);

    k = chunker.k;

    
    cfprintf('\nComputing det w/ zk = %5.2e %5.2ei using FLAM ...\n', ...
        real(zk),imag(zk))
    cfprintf('precomputing tri-diagonal part ...\n')
    start = tic; [stokestd,stokesinds] = ostokesmat_pre(chunker,nchs, ...
        zk,cs,cd,intparams);
    t3 = toc(start);
    cfprintf('time %5.2e\n',t3)
    
    kern = @(s,t,sn,tn,sw,tw,slf) ostokes_pxy_kern(s,t,sn,tn,sw,tw,slf, ...
        zk,cs,cd);
    pxyfun = @(x,slf,nbr,l,ctr)pxyfun_sq(kern,proxy,pnorm,pw, ...
        x,xnorm,whtsflam,slf,nbr,l,ctr);
    matfun = @(i,j) ostokes_matfun_wpre(xflam,xnorm,whtsflam,i,j,zk,cs,cd,...
        k,nchs,ncomp,stokestd,stokesinds);
        
    cfprintf('computing skeletonization w/ precomp tri-diagonal ...\n')
    start = tic; F_wpre = rskelf(matfun,xflam,200,1e-14,pxyfun); 
    t4 = toc(start);
    cfprintf('time %5.2e\n',t4)

    cfprintf('computing determinant ...\n')
    start = tic; d = exp(rskelf_logdet(F_wpre)); t5 = toc(start);
    cfprintf('time %5.2e\n',t5)

    cfprintf('total time %6.3e\n',t3+t4+t5)
    cfprintf('determinant = %5.2e\n',d)
    
    varargout{1} = F_wpre;


else
    cfprintf('\nComputing det w/ zk = %5.2e %5.2ei using dense routs...\n',...
        real(zk),imag(zk))
    cfprintf('forming full matrix ...\n')
    start = tic; sysmat2 = ostokesmat(chunker,zk,cs,cd,intparams); 
    t1 = toc(start);
    cfprintf('time %5.2e\n',t1)

    cfprintf('computing determinant...\n')
    start = tic; d = det(sysmat2); t2 = toc(start);
    cfprintf('time %5.2e\n',t2)

    cfprintf('total time %6.3e\n',t1+t2)
    
    varargout{1} = sysmat2;
end


end

function cfprintf1(verb,msg,varargin)
if (verb)
    fprintf(msg,varargin{:});
end
end