%%
detchebs_all = cell(5,1);
detchebs2_all = cell(5,1);

rtrue = 13.48025717955054967695879;
rspur = 14.795951782351260746661471;
errs = zeros(5,1);
errs2 = zeros(5,1);
npts = zeros(5,1);
for i=1:5
    fileout = ['example_annulus_001_' int2str(i) '_convplots.mat'];
    load(fileout);
    detchebs_all{i} = detchebs{1};
    detchebs2_all{i} = detchebs2{2};
    rts = roots(detchebs_all{i},'complex');
    
    rts2 = roots(detchebs2_all{i},'complex');
    npts(i) = sum(nchs)*chunker.k;
    errs(i) = min(abs(rts-rtrue));
    errs2(i) = min(abs(rts2-rspur));
end


figure(1)
semilogy(npts,errs,'k.','MarkerSize',20)
grid on
xlim([150,400])
ylim([10^-12,10^-4])
saveas(gcf,'trueconv.pdf')

figure(2)
semilogy(npts,errs2,'k.','MarkerSize',20)
grid on
xlim([150,400])
ylim([10^-12,10^-4])
saveas(gcf,'spurconv.pdf')

figure(3)
rmax = 4330;
plot(real(detchebs{2})/rmax,'k-','LineWidth',1.5);
hold on
plot(imag(detchebs{2})/rmax,'k--','LineWidth',1.5);
zval = detchebs{2}(rspur)/rmax;
plot(rspur,real(zval),'k.','MarkerSize',15);
plot(rspur,imag(zval),'k.','MarkerSize',15);
grid on;
saveas(gcf,'comb_det_14-15.pdf')

figure(4)
rmax = 4.98E04;
plot(real(detchebs2{2})/rmax,'k-','LineWidth',1.5);
hold on
plot(imag(detchebs2{2})/rmax,'k--','LineWidth',1.5);

zval = detchebs2{2}(rspur);
plot(rspur,real(zval),'k.','MarkerSize',15);
plot(rspur,imag(zval),'k.','MarkerSize',15);
grid on;
saveas(gcf,'dl_det_14-15.pdf')


