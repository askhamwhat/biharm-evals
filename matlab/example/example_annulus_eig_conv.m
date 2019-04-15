
% EXAMPLE_ANNULUS_EIG_CONV
% 
% compute a chebfun representation of the fredholm determinant 
% of the system I-2D-2iS and I-2D for the annulus 1<r<1.7
% on the interval [13,14] (for testing convergence to the true
% eigenvalue, rtrue below) 
% and on the interval [14,15] (for testing convergence
% of I-2D representation to spurious eigenvalue rspur below)
% save chebfun, domain, and some settings to a file

fileout = ['mat-files/annulus_eig_conv.mat'];

rtrue = 13.48025717955054967695879; %true eigenvalue
rspur = 14.79595178235126074666147; %spurious eigenvalue (Interior Neumann eigenvalue for unit disk)


% annulus params
r1 = 1.0;
r2 = 1.7;


ncells = 2;
chebabs = cell(ncells,1);
chebabs{1} = [13 14];
chebabs{2} = [14 15];

addpaths_loc();
nnchi = 5;
nchi_start = 4;

detchebs_comb = cell(ncells,nnchi); 
detchebs_dl = cell(ncells,nnchi);
t1s_comb = cell(ncells,nnchi);
t1s_dl = cell(ncells,nnchi);
chunkers = cell(2,nnchi);

errs_true_comb = zeros(nnchi,1);
errs_true_dl = zeros(nnchi,1);
errs_spur_dl = zeros(nnchi,1);
npts = zeros(nnchi,1);


% Set options for FLAM
opts = []; opts.FLAM = 1; opts.verb = true;

% Set chebfun prefs
p = chebfunpref; p.chebfuneps = 1.0e-13;
p.splitting=0; p.maxLength = 257;
t_all = tic;
for i=1:nnchi
    % setup geometry
    nchi = nchi_start+i-1;
    ncho = ceil(nchi/r1*r2)+1;

    chunkert = circle_chunks(ncho,r2);
    chunkers{1,i} = chunkert;


    chunkert = circle_chunks(nchi,r1);
    chunkert = chunkreverse(chunkert);
    chunkert = chunksort(chunkert);
    chunkers{2,i} = chunkert;
        
    [chunker,nchs] = chunkermerge({chunkers{:,i}});
    
    npts(i) = (nchi+ncho)*16;

    % Compute determinant using combined field reprensetation

    cd = -2.0 + 1i*0.0;
    cs = 0.0 - 1i*2.0;
    
    opts = []; opts.FLAM = 1; opts.verb = true;
    detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
        cd,opts);   

    
    for j = 1:length(chebabs)
        fprintf('running on interval [ %5.2e %5.2e ]\n',chebabs{j}(1), ...
            chebabs{j}(2))
        start = tic; detchebs_comb{j,i} = ...
            chebfun(detfun,chebabs{j},p); 
        t1s_comb{j,i} = toc(start);
        fprintf('%5.2e time for chebfun build\n',t1s_comb{j,i})
    end
    
    rts = roots(detchebs_comb{1,i},'complex');
    errs_true_comb(i) = min(abs(rts-rtrue))/abs(rtrue);



    % Now compute determinant using double layer representation
    cd = -2.0 + 1i*0.0;
    cs = 0;
    ncomp = length(nchs);

    detfun = @(zk) ostokes_determinant(zk,chunker,nchs,cs, ...
        cd,opts);   


    
    for j = 1:length(chebabs)
        fprintf('running on interval [ %5.2e %5.2e ]\n',chebabs{j}(1), ...
        chebabs{j}(2))
        start = tic; detchebs_dl{j,i} = ...
            chebfun(detfun,chebabs{j},p); 
        t1s_dl{j,i} = toc(start);
        fprintf('%5.2e time for chebfun build\n',t1s_dl{j,i})
    end
    rts = roots(detchebs_dl{1,i},'complex');
    errs_true_dl(i) = min(abs(rts-rtrue))/abs(rtrue);
    
    rts = roots(detchebs_dl{2,i},'complex');
    errs_spur_dl(i) = min(abs(rts-rspur))/abs(rspur);
end
total_time = toc(t_all);
fprintf('%5.2e time for computing all determinants\n',total_time)

% Post processing

figure(1)
loglog(npts,errs_true_comb,'k.','MarkerSize',20), hold on,
loglog(npts,(npts/90).^(-20),'k--')
grid on
xlim([150,400])
ylim([10^-13,10^-4])
saveas(gcf,'res-files/eig-conv-true-comb.pdf')

figure(2)
loglog(npts,errs_true_dl,'k.','MarkerSize',20), hold  on,
loglog(npts,(npts/90).^(-20),'k--')
grid on
xlim([150,400])
ylim([10^-13,10^-4])
saveas(gcf,'res-files/eig-conv-true-dl.pdf')

figure(3)
loglog(npts,errs_spur_dl,'k.','MarkerSize',20), hold on,
loglog(npts,(npts/105).^(-20),'k--')
grid on
xlim([150,400])
ylim([10^-13,10^-4])
saveas(gcf,'res-files/eig-conv-spur-dl.pdf')



figure(4)
rmax = max(abs(detchebs_comb{2,nnchi}));
plot(real(detchebs_comb{2,nnchi})/rmax,'k-','LineWidth',1.5);
hold on
plot(imag(detchebs_comb{2,nnchi})/rmax,'k--','LineWidth',1.5);
zval = detchebs_comb{2,nnchi}(rspur)/rmax;
plot(rspur,real(zval),'k.','MarkerSize',15);
plot(rspur,imag(zval),'k.','MarkerSize',15);
grid on;
saveas(gcf,'res-files/comb_det_14-15.pdf')

figure(5)
rmax = max(abs(detchebs_dl{2,nnchi}));
plot(real(detchebs_dl{2,nnchi})/rmax,'k-','LineWidth',1.5);
hold on
plot(imag(detchebs_dl{2,nnchi})/rmax,'k--','LineWidth',1.5);

zval = detchebs_dl{2,nnchi}(rspur);
plot(rspur,real(zval),'k.','MarkerSize',15);
plot(rspur,imag(zval),'k.','MarkerSize',15);
grid on;
saveas(gcf,'res-files/dl_det_14-15.pdf')

% Save all results in a mat-file
save(fileout,'chebabs','chunkers','detchebs_comb','detchebs_dl',...
    'errs_spur_dl','errs_true_comb','errs_true_dl','npts','r1','r2',...
    'rtrue','rspur','opts','nnchi','nchi_start','t1s_comb','t1s_dl','p');

