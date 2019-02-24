function fints = chunkerintkern(chunker,kern,ndims,dens,targs,opts)
%CHUNKERINTKERN compute the comvolution of the integral kernel with
% the density defined on the chunk geometry. 
%
% input:
%   chunker - chunks description of curve
%   kern - integral kernel taking inputs kern(s,t,sn,tn) where 
%          s is a source, sn is the normal at the source, t is a target
%          and tn is the normal at the target
%   ndims - input and output dimensions of the kernel (ndims(1) dimension
%           output, ndims(2) dimension of input)
%   dens - density on boundary, should have size ndims(2) x k x nch
%          where k = chunker.k, nch = chunker.nch
%   targs - targ(1:2,i) gives the coords of the ith target
%   opts - structure for setting various parameters
%       opts.targsn - if provided, the normals at the targets
%       opts.usesmooth - if = 1, then just use the smooth integration
%          rule for each chunk. otherwise, adaptive integration 
%          (quadgk) is used (default 0)
%       opts.quadgkparams - if non-empty this is a cell structure
%       containing string,value pairs to be sent to quadgk (default {})
%
% output:
%   fints - ndims(1) x nt array of integral values
%
% see also QUADGK

if nargin < 6
    opts = [];
end

[k2,nt] = size(targs);
assert(k2==2,'targs array of incorrect shape, should be (2,nt)');

if ~isfield(opts,'usesmooth'); opts.usesmooth = false; end
if ~isfield(opts,'quadgkparams'); opts.quadgkparams = {}; end
if ~isfield(opts,'targsn'); opts.targsn = zeros(2,nt); end
if ~isfield(opts,'verb'); opts.verb= false; end

targsn = opts.targsn;

[k2,nt2] = size(targsn);

assert(k2==2 && nt2==nt,...
    'opts.targsn should have same dimensions as targs');

k = chunker.k;
nch = chunker.nch;

assert(numel(dens) == ndims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,ndims(2),k,nch);

[~,w,u,~] = legeexps(k);

rnorms = chunknormals(chunker);

fints = zeros(ndims(1)*nt,1);

if opts.usesmooth
    % assume smooth weights are good enough
    for i = 1:nch
        densvals = dens(:,:,i); densvals = densvals(:);
        dsdtdt = sqrt(chunker.ders(1,:,i).^2 + ...
            chunker.ders(2,:,i).^2);
        dsdtdt = dsdtdt(:).*w(:)*chunker.hs(i);
        dsdtdt = repmat( (dsdtdt(:)).',ndims(2),1);
        densvals = densvals.*(dsdtdt(:));
        kernmat = kern(chunker.chunks(:,:,i),targs, ...
            rnorms(:,:,i),targsn);

        fints = fints + kernmat*densvals;
    end

else
    % using adaptive quadrature
    [xc,yc,xpc,ypc] = chunkerexps(chunker);
    for i = 1:nch
        if opts.verb; fprintf('chunk %d integral\n',i); end
        xci = xc(:,i); yci = yc(:,i); 
        xpci = xpc(:,i); ypci = ypc(:,i);
        densvals = dens(:,:,i); densvals = densvals.';
        densc = u*densvals; % each column is set of coefficients
                        % for one dimension of density on chunk
        for j = 1:nt
            indj = (j-1)*ndims(1);
            for l = 1:ndims(1)
                ind = indj+l;
                temp = chunkerintchunk_kernfcoefs(kern,ndims,l,...
                    densc,xci,yci,xpci,ypci,targs(:,j),targsn(:,j));

                fints(ind) = fints(ind) + temp*chunker.hs(i);
            end
        end
    end
end
