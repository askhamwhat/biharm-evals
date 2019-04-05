
function chnkro = chunkrefine(chnkri,opts)
%CHUNKREFINE refine the chunks structure according to any of a number of
% rules. This routine takes the chnkr to be resolved, i.e. the underlying
% curve is given to be exactly the chnkr.
%
% Input:
%   chnkr - chunker object
%   opts - options structure
%       opts.nchmax = maximum number of chunks on refined chunker
%                   (10*chnkri.nch)
%       opts.lvlr = level restriction flag (true), if true, enforce that no
%                   two adjacent chunks should differ in length by more 
%                   than a factor of two
%       opts.maxchunklen = maximum chunk length (Inf). enforce that no
%                   chunk is larger than this maximum length
%       opts.farfac = enforce that for each chunk, any non-adjacent 
%                   chunk is closer than farfac*chunklength (Inf)
%                   NOTE: a value of 0.5 would be like ensuring that 
%                   the level restriction holds for near intersecting 
%                   pieces of boundary
%       opts.nover = oversample boundary nover times (0)
%       opts.maxiter_lvlr = number of iterations allowed when attempting
%                           level restriction (1000)
%
% TODO implement farfac option
%
%

if nargin < 1
    opts = [];
end

lvlr = true;
maxchunklen = Inf;
farfac = Inf;
nover = 0;

maxiter_lvlr=1000;
maxiter_maxlen=1000;

nchmax = 10*chnkri.nch;

if isfield(opts,'lvlr'); lvlr = opts.lvlr; end
if isfield(opts,'maxchunklen'); maxchunklen = opts.maxchunklen; end
if isfield(opts,'farfac'); farfac = opts.farfac; end
if isfield(opts,'nchmax'); nchmax = opts.nchmax; end
if isfield(opts,'nover'); nover= opts.nover; end


if farfac < Inf
    warning('far-field refinement not implemented')
end


nch = chnkri.nch;
k = chnkri.k;
[ndim,~] = size(chnkri.chunks);

% compute lengths of chunks at start and update along the way

chunklens = zeros(nchmax,1);
whts = chunkwhts(chnkri);
chunklens(1:nch) = sum(whts,1);

chnkro = [];
chnkro.chunks = zeros(ndim,k,nchmax);
chnkro.ders = zeros(ndim,k,nchmax);
chnkro.ders2 = zeros(ndim,k,nchmax);
chnkro.adjs = zeros(2,nchmax);
chnkro.hs = zeros(nchmax,1);

chnkro.chunks(:,:,1:nch) = chnkri.chunks(:,:,1:nch);
chnkro.ders(:,:,1:nch) = chnkri.ders(:,:,1:nch);
chnkro.ders2(:,:,1:nch) = chnkri.ders2(:,:,1:nch);
chnkro.adjs = chnkri.adjs(:,1:nch);
chnkro.hs(1:nch) = chnkri.hs(1:nch);

chnkro.nch = nch;
chnkro.k = k;

% 

[~,w] = legeexps(k);

% maximum chunklength

if maxchunklen < Inf
    for ijk = 1:maxiter_maxlen

        nchold=chnkro.nch;
        ifdone=1;

        for i = 1:nchold

            rlself = chunklens(i);

    %       only check if self is sufficiently small

            if rlself > maxchunklen

    %       split chunk i now, and recalculate nodes, ders, etc

                if (chnkro.nch + 1 > nchmax)
                    error('too many chunks')
                end

                chnkro = chunksplit1(chnkro,i);

                % update chunklens 

                nch = chnkro.nch;
                dersi = chnkro.ders(:,:,i);
                hi = chnkro.hs(i);
                dersnch = chnkro.ders(:,:,nch);
                hnch = chnkro.hs(nch);

                dsdti = sum(dersi.^2,1)*hi;
                dsdtnch = sum(dersnch.^2,1)*hnch;

                chunklens(i) = dot(dsdti,w);
                chunklens(nch) = dot(dsdtnch,w);

                ifdone=0;

            end
        end

        if (ifdone == 1)
            break;
        end

    end
end

% level restriction

if lvlr
    for ijk = 1:maxiter_lvlr

        nchold=chnkro.nch;
        ifdone=1;

        for i = 1:nchold
            i1=chnkro.adjs(1,i);
            i2=chnkro.adjs(2,i);

            rlself = chunklens(i);

            rl1=rlself;
            rl2=rlself;

            if (i1 > 0)
                rl1 = chunklens(i1);
            end
            if (i2 > 0)
                rl2 = chunklens(i2);
            end

    %       only check if self is larger than either of adjacent blocks,
    %       iterating a couple times will catch everything

            sc = 2.05d0;
            if (rlself > sc*rl1 || rlself > sc*rl2)

    %       split chunk i now, and recalculate nodes, ders, etc

                if (chnkro.nch + 1 > nchmax)
                    error('too many chunks')
                end

                chnkro = chunksplit1(chnkro,i);

                % update chunklens 

                nch = chnkro.nch;
                dersi = chnkro.ders(:,:,i);
                hi = chnkro.hs(i);
                dersnch = chnkro.ders(:,:,nch);
                hnch = chnkro.hs(nch);

                dsdti = sum(dersi.^2,1)*hi;
                dsdtnch = sum(dersnch.^2,1)*hnch;

                chunklens(i) = dot(dsdti,w);
                chunklens(nch) = dot(dsdtnch,w);

                ifdone=0;

            end
        end

        if (ifdone == 1)
            break;
        end

    end
end

% oversample 

for ijk = 1:nover

    nchold=chnkro.nch;

    for i = 1:nchold

%       split chunk i now, and recalculate nodes, ders, etc
        if (chnkro.nch + 1 > nchmax)
            error('too many chunks')
        end

        chnkro = chunksplit1(chnkro,i);

        % update chunklens 

        nch = chnkro.nch;
        dersi = chnkro.ders(:,:,i);
        hi = chnkro.hs(i);
        dersnch = chnkro.ders(:,:,nch);
        hnch = chnkro.hs(nch);

        dsdti = sum(dersi.^2,1)*hi;
        dsdtnch = sum(dersnch.^2,1)*hnch;

        chunklens(i) = dot(dsdti,w);
        chunklens(nch) = dot(dsdtnch,w);

    end

end

nch = chnkro.nch;
chnkro.chunks = chnkro.chunks(:,:,1:nch);
chnkro.ders = chnkro.ders(:,:,1:nch);
chnkro.ders2 = chnkro.ders2(:,:,1:nch);
chnkro.adjs = chnkro.adjs(:,1:nch);
chnkro.hs = chnkro.hs(1:nch);