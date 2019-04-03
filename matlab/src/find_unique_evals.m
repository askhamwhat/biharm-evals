function iplot = find_unique_evals(rts,sings,abs_eps,imag_eps)
%FIND_UNIQUE_EVALS returns indices corresponding to the roots 
% which seem like they correspond to unique eigenvalues to the
% desired precision
%

[~,isort] = sort(real(rts));
iplot1 = isort(abs(imag(rts(isort))) < imag_eps);

if length(iplot1) < 2
    iplot = iplot1;
else
    diffs = diff(rts(iplot1));
    iremove = find(abs(diffs) < abs_eps);
    iremover = iplot1( iremove(sings(iplot1(iremove))<sings(iplot1(iremove+1)))+1);
    iremovel = iplot1( iremove(sings(iplot1(iremove))>= sings(iplot1(iremove+1))));
    iplot = setdiff(iplot1(:),[iremover(:);iremovel(:)]);
    [~,isort] = sort(real(rts(iplot)));
    iplot = iplot(isort);
end



