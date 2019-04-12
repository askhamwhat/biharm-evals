
function [Kpxy,nbr] = pxyfun_sqr(rc,rx,cx,slf,nbr,l,ctr, ...
    kernr,kernc,rxn,cxn,p,pn)
%PXYFUN_SQR proxy surface utility for square proxy surface - rect matrix
%
%
% Input 
%
% rc - char 'R' or 'C' determines row or column compression
% rx - 2 x nrows array of row points
% cx - 2 x nrows array of col points
% slf - relevant indices in rx or cx
% nbr - indices of candidate neighbor points
% l - length of a box at this level of the tree
%
% kernr should be a kernel function of the form 
% submat = kernr(src, targ, srcn, targn, slf)
% kernc should be a kernel function of the form 
% submat = kernc(src, targ, srcn, targn, slf)
% rxn and cxn and are corresponding normals (same format as rx and cx)
% p should be a reference set of proxy points
% pn should be corresponding outward normals
%
% Output
%
% Kpxy - matrix for proxy surface compression of row or col points
% nbr - subset of input nbr indices which are within the proxy
% surface
%

  pxy = bsxfun(@plus,p*l,ctr(:));
  if strcmpi(rc,'r')
    Kpxy = kernr(pxy,rx(:,slf),pn,rxn(:,slf),slf)*l/length(p);
    dx = cx(1,nbr) - ctr(1);
    dy = cx(2,nbr) - ctr(2);
  elseif strcmpi(rc,'c')
    Kpxy = kernc(cx(:,slf),pxy,cxn(:,slf),pn,slf);
    dx = rx(1,nbr) - ctr(1);
    dy = rx(2,nbr) - ctr(2);
  end

  dist = max(abs([dx; dy]),[],1);
  nbr = nbr(dist/l < 1.5);
end
