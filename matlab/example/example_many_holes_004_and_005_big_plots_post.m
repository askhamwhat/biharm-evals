%EXAMPLE_MANY_HOLES_004_and_005_BIG_PLOTS_POST
%
% process saved vorticity and singular value data

clear all

% addpath('../src','../test','../../mwrap');

addpaths_loc();

filein1 = 'example_many_holes_004.mat';
filein2 = 'example_many_holes_004_big_plots.mat';

filein0 = 'example_many_holes_005_big_plots.mat';
filein00 = 'example_many_holes_005.mat';

load(filein0); load(filein00);
ss5 = ss; rts5 = rts; detchebs5 = detchebs;
vorts5 = vorts; xnulls5 = xnulls;

load(filein1); load(filein2);

chunkers = multichunksort(chunker,nchs,true);
%%

sings = zeros(length(ss),1);
for i = 1:length(ss)
    sings(i) = ss{i}(1);
end

sings5 = zeros(length(ss5),1);
for i = 1:length(ss5)
    sings5(i) = ss5{i}(1);
end

%


%% 

figure(1)
clf

mrow = 8; ncol = 10;
[ha, pos] = tight_subplot(mrow,ncol,[.01 .01],[.01 .01],[.01 .01]);

nsamp=10;

for ii = 1:mrow*ncol
    axes(ha(ii))
    
    if ii > length(rts)
        axis off
    else
        i = ii;

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
end

%%

%%

%%


%%

fileout = 'ex_many_holes_004_sings.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(sings)); ...
    reshape(abs(sings),1,length(sings))]);
fclose(fid);

%%

for i = 1:length(chebabs)
    fileout = sprintf('ex_many_holes_004_coeffs_%d.tex',i);
    fid = fopen(fileout,'w'); 
    coefft = chebcoeffs(detchebs{i});
    fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(coefft)); ...
        reshape(abs(coefft),1,length(coefft))]);
    fclose(fid);
end
for i = 1:length(chebabs)
    fileout = sprintf('ex_many_holes_004_coeffs_scaled_%d.tex',i);
    fid = fopen(fileout,'w'); 
    coefft = chebcoeffs(detchebs{i});
    fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(coefft)); ...
        reshape(abs(coefft)/abs(coefft(1)),1,length(coefft))]);
    fclose(fid);
end


%%

for i = 1:length(chunkers)
    chnkr = chunkers{i};
    fileout = sprintf('ex_many_holes_004_bdry_%d.tex',i);
    fid = fopen(fileout,'w'); 
    fprintf(fid,'%7.4e %7.4e\n',[chnkr.chunks(1,:); ...
        chnkr.chunks(2,:)]);
    fclose(fid);
end    


%% 

figure(1)
clf
hold on
for i = 1:length(detchebs)
    semilogy(abs(detchebs{i}),'b')
end
axis tight

%%
for i = 1:length(detchebs)
    cf{i} = chebcoeffs(detchebs{i});
    chebabs{i}
    max(abs(cf{i}))/abs(cf{i}(1))
    min(abs(cf{i}))/max(abs(cf{i}))
    min(abs(cf{i}))/abs(cf{i}(1))
end

%%

for i = 1:length(detchebs)
    ab = chebabs{i};
    a = ab(1); b = ab(2);
    iiii = and(rts > a,rts < b);
    subplot(4,4,i)
    plot(real(diff(detchebs{i})),'r')
    hold on 
    plot(imag(diff(detchebs{i})),'b')
    plot(rtsp(iiii),zeros(nnz(iiii),1),'go')
end

%% error estimate


reldir = zeros(length(rts),1);

for i = 1:length(detchebs)
    ab = chebabs{i}
    a = ab(1); b = ab(2);
    iiii = and(rts > a,rts < b);
    
    ff = detchebs{i};
    dff = diff(ff);
    reldir(iiii) = norm(ff,'inf')./abs(dff(rts(iiii)));
end

%

figure(1)
clf
semilogy(abs(sings))
hold on

sing_err_estimate = max(reldir*p.chebfuneps, ...
    sqrt(chunker.nch*chunker.k)*1e-14);

semilogy(sing_err_estimate)



%%

fileout = 'ex_many_holes_004_sings_est.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(rts)); ...
    reshape(sing_err_estimate,1,length(rts))]);
fclose(fid);

fileout = 'ex_many_holes_005_sings.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*find(any(abs(bsxfun(@plus,rts.',-rts5(:)))<1e-7))
; ...
    reshape(abs(sings5),1,length(sings5))]);
fclose(fid);


%% 


nperab = 400;
xxtot = [];
yytot = [];
for i = 1:length(detchebs)
    ab = chebabs{i};
    xx = linspace(ab(1),ab(2),nperab);
    xxtot = [xxtot; xx(:)];
    yytot = [yytot; detchebs{i}(xx(:))];
end
    
figure(1)
clf
hold on
plot(xxtot,real(yytot),'r')
plot(xxtot,imag(yytot),'b')

figure(2)
clf
semilogy(xxtot,abs(yytot))
