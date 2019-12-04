%EXAMPLE_BARBELL_001_BIG_PLOTS_POST
%
% process saved vorticity and singular value data

% addpath('../src','../test','../../mwrap');

addpaths_loc();

filein1 = 'example_barbell_001.mat';
filein2 = 'example_barbell_001_big_plots.mat';

load(filein1); load(filein2);

chunkers = multichunksort(chunker,nchs,true);
%%

sings = zeros(length(ss),1);
for i = 1:length(ss)
    sings(i) = ss{i}(1);
end


iplot = 1:length(rts);

%

%% 

figure(1)
clf

mrow = 12; ncol = 10;
[ha, pos] = tight_subplot(mrow,ncol,[.01 .01],[.01 .01],[.01 .01]) ;

nsamp=10;

for ii = 1:min(length(iplot),mrow*ncol)
    i = iplot(ii);
    axes(ha(ii))
    vort = nan(size(vorts{i}));
    
    vort(in) = vorts{i}(in); iin = find(in); isel = randperm(length(iin),nsamp);
    vort = vort/(mean(vort(iin(isel))));
    norm(imag(vort(in)),'fro')/norm(vort(in),'fro')
    vort = real(vort);
    vmax = max(max(abs(vort))); vort = vort/vmax;
    hold off
    h = pcolor(xx,yy,vort); set(h,'EdgeColor','none');
    hold on
    for i = 1:length(chunkers)
        chnkr = chunkers{i};
        plot(chnkr.chunks(1,:),chnkr.chunks(2,:),'-k')
    end
    
    axis equal
    colormap(redblue)
    caxis([-1 1])

    set(gca,'XTick',[], 'YTick', [])
    axis equal    
    axis tight
    axis equal
end


%% PRINT OUT SINGULAR VALUES TO TEX FILES

fileout = 'ex_barbell_001_sings.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(iplot)); ...
    reshape(abs(sings(iplot)),1,length(iplot))]);
fclose(fid);

%%

for i = 1:length(chebabs)
    fileout = sprintf('ex_barbell_001_coeffs_%d.tex',i);
    fid = fopen(fileout,'w'); 
    coefft = chebcoeffs(detchebs{i});
    fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(coefft)); ...
        reshape(abs(coefft),1,length(coefft))]);
    fclose(fid);
end
for i = 1:length(chebabs)
    fileout = sprintf('ex_barbell_001_coeffs_scaled_%d.tex',i);
    fid = fopen(fileout,'w'); 
    coefft = chebcoeffs(detchebs{i});
    fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(coefft)); ...
        reshape(abs(coefft)/abs(coefft(1)),1,length(coefft))]);
    fclose(fid);
end

%% 
% 
% figure(1)
% clf
% hold on
% for i = 1:length(detchebs)
%     semilogy(abs(detchebs{i}),'b')
% end
% axis tight
% 

%% 

%%

for i = 1:length(chunkers)
    chnkr = chunkers{i};
    fileout = sprintf('ex_barbell_001_bdry_%d.tex',i);
    fid = fopen(fileout,'w'); 
    fprintf(fid,'%7.4e %7.4e\n',[chnkr.chunks(1,:); ...
        chnkr.chunks(2,:)]);
    fclose(fid);
end    


%% 

% 
% 
% nperab = 400;
% xxtot = [];
% yytot = [];
% for i = 1:length(detchebs)
%     ab = chebabs{i};
%     xx = linspace(ab(1),ab(2),nperab);
%     xxtot = [xxtot; xx(:)];
%     yytot = [yytot; detchebs{i}(xx(:))];
% end
%     
% figure(1)
% clf
% hold on
% plot(xxtot,real(yytot),'r')
% plot(xxtot,imag(yytot),'b')
% 
% figure(2)
% clf
% semilogy(xxtot,abs(yytot))