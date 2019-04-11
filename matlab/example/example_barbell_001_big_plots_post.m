%EXAMPLE_BARBELL_001_BIG_PLOTS_POST
%
% process saved vorticity and singular value data

addpath('../src','../test','../../mwrap');

filein1 = 'example_barbell_001.mat';
filein2 = 'example_barbell_001_big_plots.mat';

load(filein1); load(filein2);

chunkers = multichunksort(chunker,nchs,true);
%%

sings = zeros(length(ss),1);
for i = 1:length(ss)
    sings(i) = ss{i}(end);
end

%%

indimagbig = abs(imag(rts)) > sqrt(p.chebfuneps);
nnz(indimagbig)
find(indimagbig)
sings(indimagbig)

%%
%

[iplot,ichkdbl] = find_unique_evals(rts,sings,1e-5,...
    sqrt(p.chebfuneps));

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

pord = 10;
q = 40;
ss2 = zeros(p,length(ichkdbl));

for ii = 1:length(ichkdbl)
    i = ichkdbl(ii,1);
    zk = real(rts(i));
    start = tic; [d,F] = ostokes_determinant(zk,chunker,nchs,cs,cd,opts);
    toc(start)
    nsys = 2*chunker.nch*chunker.k;
    start = tic; xnull = rskelf_nullvec(F,nsys,pord,q); toc(start)
    ynull = rskelf_mv(F,xnull);
    ss2(:,ii) = sqrt(sum(abs(ynull).^2,1));
end

%%

% compare computed vorticity at suspicious sings

projs = zeros(length(ichkdbl),1);

for ii =  1:length(ichkdbl)
    i1 = ichkdbl(ii,1); i2 = ichkdbl(ii,2);
    v1 = vorts{i1}(in); v2 = vorts{i2}(in);
    v1 = v1/norm(v1,'fro');
    projs(ii) = norm(v2-v1*v1'*v2);
end


%%

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

figure(1)
clf
hold on
for i = 1:length(detchebs)
    semilogy(abs(detchebs{i}),'b')
end
axis tight


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