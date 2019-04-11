%EXAMPLE_MANY_HOLES_004_BIG_PLOTS_POST
%
% process saved vorticity and singular value data

clear all

addpath('../src','../test','../../mwrap');

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

indimagbig = abs(imag(rts)) > sqrt(p.chebfuneps);
nnz(indimagbig)
find(indimagbig)
sings(indimagbig)

%

[iplot,ichkdbl] = find_unique_evals(rts,sings,1e-5, ...
    sqrt(p.chebfuneps));

%% 

figure(1)
clf

mrow = 8; ncol = 10;
[ha, pos] = tight_subplot(mrow,ncol,[.01 .01],[.01 .01],[.01 .01]);

nsamp=10;

for ii = 1:mrow*ncol
    axes(ha(ii))
    
    if ii > length(iplot)
        axis off
    else
        i = iplot(ii);

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

tempi = [97, 98];

mmrow = 1; nncol =2;

for ii = 1:2
    i =  tempi(ii);
    subplot(mmrow,nncol,ii)
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

%%

% Run through and check for double roots at suspicious sings

sings2 = zeros(length(ss),1);
for i = 1:length(ss)
    sings2(i) = ss{i}(2);
end

[sings(ichkdbl(:,1)), sings2(ichkdbl(:,1))]

%%

% compare computed vorticity at suspicious sings

projs = zeros(length(ichkdbl),1);
projsx = zeros(length(ichkdbl),1);

for ii =  1:length(ichkdbl)
    i1 = ichkdbl(ii,1); i2 = ichkdbl(ii,2);
    v1 = vorts{i1}(in); v2 = vorts{i2}(in);
    x1 = xnulls{i1}(:,1); x2 =xnulls{i2}(:,1);
    v1 = v1/norm(v1,'fro');
    projs(ii) = norm(v2-v1*v1'*v2,'fro')/norm(v2,'fro');
    x1 = x1/norm(x1,'fro');
    x2 = x2/norm(x2,'fro');
    projsx(ii) = norm(x2-x1*x1'*x2);
end


%%

fileout = 'ex_many_holes_004_sings.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(iplot)); ...
    reshape(abs(sings(iplot)),1,length(iplot))]);
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

rtsp = rts(iplot);

for i = 1:length(detchebs)
    ab = chebabs{i};
    a = ab(1); b = ab(2);
    iiii = and(rtsp > a,rtsp < b);
    subplot(4,4,i)
    plot(real(diff(detchebs{i})),'r')
    hold on 
    plot(imag(diff(detchebs{i})),'b')
    plot(rtsp(iiii),zeros(nnz(iiii),1),'go')
end

%%


rtsp = rts(iplot);
reldir = zeros(length(iplot),1);

for i = 1:length(detchebs)
    ab = chebabs{i}
    a = ab(1); b = ab(2);
    iiii = and(rtsp > a,rtsp < b);
    
    ff = detchebs{i};
    dff = diff(ff);
    reldir(iiii) = norm(ff,'inf')./abs(dff(rtsp(iiii)));
end

%

figure(1)
clf
semilogy(abs(sings(iplot)))
hold on

sing_err_estimate = max(reldir*p.chebfuneps, ...
    sqrt(chunker.nch*chunker.k)*1e-14);

semilogy(sing_err_estimate)



%%

fileout = 'ex_many_holes_004_sings_est.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*(1:length(iplot)); ...
    reshape(sing_err_estimate,1,length(iplot))]);
fclose(fid);

fileout = 'ex_many_holes_005_sings.tex';
fid = fopen(fileout,'w'); 
fprintf(fid,'%7.4e %7.4e\n',[1.0*find(any(abs(bsxfun(@plus,rts(iplot).',-rts5(:)))<1e-7))
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
