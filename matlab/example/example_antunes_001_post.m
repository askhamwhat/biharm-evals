%EXAMPLE_ANTUNES_001_POST
%
% post-processing for example_antunes_001
%
% reads in chebfun and writes info for figures to 
% various files 

filein='example_antunes_001.mat';
load(filein);

rdetcheb = real(detcheb);
idetcheb = imag(detcheb);
detcheb_roots = roots(detcheb,'complex'); 
detcheb_roots = detcheb_roots(:).';
detcheb_coeffs = chebcoeffs(detcheb); 
detcheb_coeffs = detcheb_coeffs(:).';

dom = detcheb.domain;

nplot = 2000;
xx = linspace(dom(1),dom(2),nplot);
yyr = rdetcheb(xx);
yyi = idetcheb(xx);

fileout = 'ex_antunes_001_plotreal.tex';
fid = fopen(fileout,'w'); fprintf(fid,'%7.4e %7.4e\n',[xx;yyr]); 
fclose(fid);

fileout = 'ex_antunes_001_plotimag.tex';
fid = fopen(fileout,'w'); fprintf(fid,'%7.4e %7.4e\n',[xx;yyi]); 
fclose(fid);

fileout = 'ex_antunes_001_roots.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[detcheb_roots;zeros(size(detcheb_roots))]); 
fclose(fid);

fileout = 'ex_antunes_001_coeffs.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(detcheb_coeffs)); ...
    abs(detcheb_coeffs)/abs(detcheb_coeffs(1))]);
fclose(fid);


