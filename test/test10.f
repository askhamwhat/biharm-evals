      implicit real *8 (a-h,o-z)
      real *8 w(1000000),xy(2,100),ck(1000),cms(2,1000),targ(2,100)

      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n

      done = 1
      pi = atan(done)*4

      rin = 1.2d0
      rout = 1.7d0

      ncomp = 2
      cms(1,1) = 3.1d0
      cms(2,1) = 3.2d0

      cms(1,2) = 0.0d0
      cms(2,2) = 0.0d0

      nsrc = 2
      pert = 0.1d0
      xy(1,1) = cms(1,1) + pert*hkrand(0)
      xy(2,1) = cms(2,1) + pert*hkrand(0)

      xy(2,1) = cms(1,2) + pert*hkrand(0)
      xy(2,2) = cms(2,2) + pert*hkrand(0)

      ck(1) = hkrand(0)
      ck(2) = hkrand(0)

      nchci = 20

      iw2 = 24

      call testhbhdir2(iw2,rin,rout,nchci,nsrc,
     1   ncomp,cms,xy,ck,ntot,erra,rcond)

      call prinf('ntot=*',ntot,1)
      call prin2('erra=*',erra,1)

      stop
      end

c------------------------------------------------------     
c
      subroutine testhbhdir2(ifile,rin,rout,nchci,
     1   nsrc,ncomp,cms,xy,ck,ntot,errl2,rcond)

      implicit real *8 (a-h,o-z)

      real *8 :: dist, xy(2,*), ck(*), xy_norm(2),cms(2,*)
      real *8 xsrc(100000),ysrc(100000)
      real *8, allocatable :: targ(:,:)
      real *8, allocatable :: targx(:),targy(:)
      real *8, allocatable :: errt(:),pvals(:),pvals2(:)
      real *8, allocatable :: chunks(:,:,:)
      real *8, allocatable :: ders(:,:,:)
      real *8, allocatable :: rn(:,:,:)
      real *8 sc
      real *8, allocatable :: tails(:,:)
      real *8, allocatable :: ders2(:,:,:),dsdta(:,:)
      real *8, allocatable :: hs(:)
      integer, allocatable :: adjs(:,:),ca(:,:)
      real *8, allocatable :: cc(:,:,:),cd(:,:,:),cd2(:,:,:),ch(:)

      real *8, allocatable :: wgeos(:)
      real *8, allocatable :: wgeos3(:)

      complex *16, allocatable :: work(:),workspace(:)
      complex *16, allocatable :: rhs(:),soln(:),soln2(:)
      complex *16, allocatable :: zdens(:,:)
      integer, allocatable :: nchs(:)

      complex *16 zk

c     matrix formation, eigenvalue calculation

      real *8, allocatable :: whts(:,:)
      complex *16, allocatable :: sysmat(:,:)
      complex *16 p1(10),p2,p3,p4
      
      real *8, allocatable :: xs(:),ys(:),cxs(:),cys(:),qwts(:)
      real *8, allocatable :: dxdt(:),dydt(:),dsdt(:),rnx(:),rny(:)
      real *8, allocatable :: rkappa(:)

      integer novers(1000), iparsindeces(1000)
      integer ifcloseds(1000)
      real *8 chsmalls(1000), tas(1000), tbs(1000)
      real *8 epss(1000), pars(10 000), pars1(100),pars2(100)

      real *8 src(2),verts(2,10000),widths(10000)
      real *8 errmax1(100),errmax2(100)
      
      complex *16, allocatable :: pot(:,:)
      complex *16, allocatable :: pottau(:,:)
      complex *16, allocatable :: grad(:,:,:)
      complex *16, allocatable :: potn(:,:)
      complex *16, allocatable :: wbdry(:,:)
      complex *16 gradtmp(2), wintex
      
      complex *16 sing(100000), q1, q2

      integer ncomp, nwiggles
      real *8 tmp, errp, rp
      real *8 errs(1000)

      complex *16 eye
      integer nd(10000)
      integer, allocatable :: isint(:)
      complex *16, allocatable :: wex(:),wval(:),wgrad(:)

      data eye/(0.0d0,1.0d0)/

      external multa, zhelmstokes_kern, fgreensdummy


      done = 1
      pi = atan(done)*4
      nchunkmax = 10000

      zk = 12.57124227636001d0

      zk = 14.0d0

      k = 16

      allocate(chunks(2,k,nchunkmax))
      allocate(tails(nchunkmax,4))
      allocate(ders(2,k,nchunkmax),dsdta(k,nchunkmax))
      allocate(rn(2,k,nchunkmax))
      allocate(ders2(2,k,nchunkmax))
      allocate(hs(nchunkmax))
      allocate(adjs(2,nchunkmax))


      allocate(cc(2,k,nchci),cd(2,k,nchci),cd2(2,k,nchci),ca(2,nchci))
      allocate(ch(nchci))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Setup geometry
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc      


      lwgeos = 1000000


      allocate(wgeos3(lwgeos),wgeos(lwgeos))
c
cc     update chunks array
c 

      iort = 1
      call getchunks(iort,k,nchci,chunks,adjs,ders,ders2,hs)

      iort = -1
      call getchunks(iort,k,nchci,cc,ca,cd,cd2,ch)

      

      do ich=1,nchci
      do j=1,k

      chunks(1,j,ich) = chunks(1,j,ich)*rout
      chunks(2,j,ich) = chunks(2,j,ich)*rout

      ders(1,j,ich) = rout*ders(1,j,ich)
      ders(2,j,ich) = rout*ders(2,j,ich)

      ders2(1,j,ich) = rout*ders2(1,j,ich)
      ders2(2,j,ich) = rout*ders2(2,j,ich)

      enddo
      enddo

      nchstart = nchci


      call prin2('rin=*',rin,1)
      call prin2('rout=*',rout,1)


      do ich=1,nchci
      i = nchstart + ich
      do j=1,k

      chunks(1,j,i) = rin*cc(1,j,ich)
      chunks(2,j,i) = rin*cc(2,j,ich)

      ders(1,j,i) = -rin*cd(1,j,ich)
      ders(2,j,i) = -rin*cd(2,j,ich)

      ders2(1,j,i) = rin*cd2(1,j,ich)
      ders2(2,j,i) = rin*cd2(2,j,ich)

      enddo

      hs(i) = ch(ich)
      adjs(1,i) = ca(1,ich) + nchstart
      adjs(2,i) = ca(2,ich) + nchstart

      enddo


c
cc      end of update chunks array
c

      nch = 2*nchci
      call prinf('adjs=*',adjs,2*nch)
      call prinf('nch=*',nch,1)
      call prinf('ncomp=*',ncomp,1)


      
c     Allocate potential, gradient and potn arrays

      allocate(nchs(ncomp))

      nchs(1) = nchci
      nchs(2) = nchci

      call prinf('nch=*',nch,1)

      call prinf('nchs=*',nchs,ncomp)


      call chunkpack(k,nch,chunks,adjs,ders,ders2,hs,wgeos,lused)
      call prinf('lused=*',lused,1)

      n = k*nch
      allocate(xs(n),ys(n))
      allocate(dsdt(n),rkappa(n),dxdt(n),dydt(n),rnx(n),rny(n))
      allocate(whts(k,nch),qwts(n))

      call chunkunpack1(wgeos,k,nch,ichunks,iadjs,iders,iders2,ihs)
      
      call chunkwhts(k,nch,wgeos(ichunks),wgeos(iders),wgeos(ihs),whts)
      call prin2('whts=*',whts,24)

      call prin2('hs=*',hs,24)

      
      i = 0
      do ich=1,nch
      do j=1,k

      i = i+1

      xs(i) = chunks(1,j,ich)
      ys(i) = chunks(2,j,ich)
      dxdt(i) = ders(1,j,ich)
      dydt(i) = ders(2,j,ich)
      dsdt(i) = sqrt(ders(1,j,ich)**2 + ders(2,j,ich)**2)
      dsdta(j,ich) = dsdt(i)
      rnx(i) = dydt(i)/dsdt(i)
      rny(i) = -dxdt(i)/dsdt(i)
      rn(1,j,ich) = rnx(i)
      rn(2,j,ich) = rny(i)
      qwts(i) = whts(j,ich)
      rkappa(i) = (dxdt(i)*ders2(2,j,ich) - dydt(i)*ders2(1,j,ich))/
     1              dsdt(i)**3

      enddo
      enddo

      ifile2 = 23
      call pyplot(ifile2,xs,ys,n,2,'a*')

      ifile2 = ifile2 + 1

      ntot = 2*k*nch + ncomp
      allocate(sysmat(ntot,ntot))

      call prinf('k=*',k,1)
      call prinf('ntot=*',ntot,1)
      call prinf('n=*',n,1)
      q1 = (1.0d0,0.0d0)
      q2 = (1.0d0,0.0d0)

      call zhbhstokesmatbuild(zk,wgeos,ncomp,nchs,cms,
     1     q1,q2,ntot,sysmat,ier)

c
cc     square root scale the matrices
c
      do ich = 1,nch
      do inode =1,k

      ipt = (ich-1)*k + inode

      do jch = 1,nch
      do jnode = 1,k

      jpt = (jch-1)*k + jnode

      rr = sqrt(whts(inode,ich)/whts(jnode,jch))
      sysmat(2*ipt-1,2*jpt-1) = sysmat(2*ipt-1,2*jpt-1)*rr
      sysmat(2*ipt-1,2*jpt) = sysmat(2*ipt-1,2*jpt)*rr
      sysmat(2*ipt,2*jpt-1) = sysmat(2*ipt,2*jpt-1)*rr
      sysmat(2*ipt,2*jpt) = sysmat(2*ipt,2*jpt)*rr

      enddo
      enddo
      enddo
      enddo

      do ipt=2*k*nch+1,ntot
      do jch=1,nch
      do jnode=1,k

      jpt = (jch-1)*k+jnode

      rr = 1.0d0/sqrt(whts(jnode,jch))
      sysmat(ipt,2*jpt-1) = sysmat(ipt,2*jpt-1)*rr
      sysmat(ipt,2*jpt) = sysmat(ipt,2*jpt)*rr

      enddo
      enddo
      enddo

      do ich=1,nch
      do inode=1,k
      ipt = (ich-1)*k+inode

      rr = sqrt(whts(inode,ich))

      do jpt = 2*k*nch+1,ntot

      sysmat(2*ipt-1,jpt) = sysmat(2*ipt-1,jpt)*rr
      sysmat(2*ipt,jpt) = sysmat(2*ipt,jpt)*rr

      enddo
      enddo
      enddo

     
      
c
c     Allocate potential, gradient and potn arrays
c
      allocate(pot(k,nch),potn(k,nch),grad(2,k,nch),pottau(k,nch))
      allocate(rhs(ntot),soln(ntot),soln2(ntot2))

      nn = k*nch
      allocate(zdens(nn,2))
      do i=1,ntot
         rhs(i) = 0.0d0
         soln(i) = 0.0d0
      enddo

c     Generate boundary data     
c
      ii = 0
      call prin2('cms=*',cms,2*ncomp)
      do i=1,nch
         do j=1,k
             pot(j,i) = 0.0d0
             pottau(j,i) = 0.0d0
             gradtmp(1) = 0
             gradtmp(2) = 0
             ntt = 1
             call gethboundarydata(zk,nsrc,xy,ck,ntt,chunks(1,j,i),
     1                            pot(j,i),gradtmp)
             potn(j,i) = gradtmp(1)*rn(1,j,i)+gradtmp(2)*rn(2,j,i)
             ii = ii + 1
             rhs(2*ii-1) = -gradtmp(2)
             rhs(2*ii) = gradtmp(1)
         enddo
      enddo

      do ich=1,nch
      do inode=1,k

      ipt = (ich-1)*k+inode
      zdens(ipt,1) = rhs(2*ipt-1)
      zdens(ipt,2) = rhs(2*ipt)

      enddo
      enddo

      do j=1,4

      errmax1(j) = 0
      errmax2(j) = 0

      do i=1,nch

      tails(i,j) = 0

      enddo
      enddo

      call chunkzres(k,nch,chunks,ders,hs,zdens(1,1),errmax1(1),
     1    errmax2(1),tails(1,1))

      call chunkzres(k,nch,chunks,ders,hs,zdens(1,2),errmax1(2),
     1    errmax2(2),tails(1,2))

c     get integral of exact solution on each
c     boundary component
      
      ich = 1
      do nbod = 1,ncomp
         wintex = 0
         do i = 1,nchs(nbod)
            do j = 1,k
               wintex = wintex + whts(j,ich)*pot(j,ich)
            enddo
            ich = ich+1
         enddo
         rhs(2*k*nch+nbod) = wintex
      enddo

      call prin2('rhs=*',rhs,24)

      call prin2('rhs end = *',rhs(2*k*nch+1),2*ncomp)
c
cc      square root scale rhs
c
      do ich=1,nch
      do inode = 1,k

      ipt = (ich-1)*k+inode
      rr= sqrt(whts(inode,ich))

      rhs(2*ipt-1) = rhs(2*ipt-1)*rr
      rhs(2*ipt) = rhs(2*ipt)*rr

      enddo
      enddo
      

      ra = 0.0d0

      do i=1,ntot

      ra = ra + abs(rhs(i))**2

      enddo

      ra = sqrt(ra)
      call prin2('l2 norm of rhs=*',ra,1)
      
c     End of generating boundary data

      do i=1,ntot

      soln(i) = 0

      enddo


      ngmrec = 200
      numit = 200

      lw = (ngmrec*2 + 4)*ntot

      allocate(work(lw))
      job = 0
      
      ier = 0
      epsg = 1.0d-14

      call cgmres(ier,ntot,sysmat,multa,p1,p2,p3,p4,rhs,epsg,numit,soln,
     1        niter,errs,ngmrec,work)


      do ich=1,nch
      do inode=1,k

      ipt = (ich-1)*k+inode

      rr = sqrt(1.0d0/whts(inode,ich))
      soln(2*ipt-1) = soln(2*ipt-1)*rr
      soln(2*ipt) = soln(2*ipt)*rr

      enddo
      enddo

cc      call prin2('soln final=*',soln,2*ntot)

      do ich=1,nch
      do inode=1,k

      ipt = (ich-1)*k+inode
      zdens(ipt,1) = soln(2*ipt-1)
      zdens(ipt,2) = soln(2*ipt)

      enddo
      enddo


      call chunkzres(k,nch,chunks,ders,hs,zdens(1,1),errmax1(3),
     1    errmax2(3),tails(1,3))

      call chunkzres(k,nch,chunks,ders,hs,zdens(1,2),errmax1(4),
     1    errmax2(4),tails(1,4))


      call prin2('errmax1=*',errmax1,4)
      call prin2('errmax2=*',errmax2,4)

c
cc      compute expansion tails

      call prinf('ier=*',ier,1)
      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)


      call prin2('soln=*',soln,24)

      call prin2('computed correction=*',soln(2*k*nch+1),2)


      nlat = 300
      nt = nlat*nlat
      allocate(targ(2,nt),wval(nt),wgrad(2*nt),wex(nt),isint(nt))
      allocate(errt(nt),pvals(nt),pvals2(nt))
      allocate(targx(nt),targy(nt))



      xmin = -rout*1.1d0
      xmax = rout*1.1d0

      ymin = -rout*1.1d0
      ymax = rout*1.1d0
      

      ii = 1

      do i=1,nlat
      do j=1,nlat

      xt = xmin + (i-1.0d0)/(nlat-1.0d0)*(xmax-xmin)
      yt = ymax + (j-1.0d0)/(nlat-1.0d0)*(ymin-ymax)



      targ(1,ii) = xt
      targ(2,ii) = yt
      wex(ii) = 0
      wval(ii) = 0
      wgrad(ii) = 0
      errt(ii) = -8

      ii = ii+1

      enddo
      enddo

      call prinf('ii=*',ii,1)

c
cc       find targets in the interior of the domain
c
      
      nnn = k*nchs(1)

      call prinf('nnn=*',nnn,1)
      call prinf('ncomp=*',ncomp,1)
      call prin2('cms=*',cms,2*ncomp)




      ntin = 0
      do i=1,nt

      isint(i) = 0
      xt = targ(1,i)
      yt = targ(2,i)

      rr1 = xt**2 + yt**2
      if(rr1.le.rout**2.and.rr1.gt.rin**2) isint(i) = 1

      if(isint(i).eq.1) then

      ntin = ntin + 1
      targx(ntin) = xt
      targy(ntin) = yt

      endif

      enddo

      call prinf('ntin=*',ntin,1)

      itypes = 1
      itypet = 2
      call pyplot2(ifile2,xs,ys,n,itypes,targx,targy,ntin,itypet,'a*')


      allocate(wbdry(k,nch))

      t1 = second()

      call wevaltargstoke2(zk,nt,targ,isint,wgeos,ncomp,nchs,cms,
     1     q1,q2,soln,wval,wbdry,wintex)

      t2 = second()

      call prin2('total time=*',t2-t1,1)

      call prin2('wbdry*',wbdry,12)
      call prin2('pot*',pot,12)

      call prin2('wval=*',wval,2)
      call prin2('cms=*',cms,2*ncomp)
      call prin2('ck=*',ck,nsrc)


      call gethboundarydata(zk,nsrc,xy,ck,nt,targ(1,1),
     1     wex,wgrad)

      do i=1,nt

      if(isint(i).ne.1) then

      wval(i) = wex(i)

      endif

 
      pvals(i) = abs(wval(i))
      pvals2(i) = abs(wex(i))

      if(isint(i).eq.1) then


      errt(i) = log(abs(wval(i)-wex(i)))/log(10.0d0)
      if(errt(i).gt.-6) then

      call prinf('i=*',i,1)
      call prin2('targ=*',targ(1,i),2)
      call prin2('wval=*',wval(i),2)
      call prin2('wex=*',wex(i),2)

      endif
      
      endif

      enddo


      errl2 = 0
      rl2 = 0
      do i = 1,nt
         rl2 = rl2 + abs(wex(i))**2
         errl2 = errl2 + abs(wex(i)-wval(i))**2
      enddo

      errl2 = sqrt(errl2/rl2)

      call prin2('errl2 = *',errl2,1)


      call prinf('ntot=*',ntot,1)

      call prin2('ck=*',ck(1),nsrc)
      call prin2('ck like terms=*',soln(2*k*nch+2),2*(ncomp-1))

      iw = ifile+ 1

      do ib=2,ncomp

      nd(ib-1) = k*nchs(ib)

      enddo

      call prinf('nd=*',nd,ncomp-1)

      is = k*nchs(1)+1

      call pyimage3(iw,nlat,nlat,pvals,ncomp-1,xs(is),ys(is),nd,xmin,
     1       xmax,ymin,ymax)


      iw = ifile+2

      call pyimage3(iw,nlat,nlat,pvals2,ncomp-1,xs(is),ys(is),nd,xmin,
     1       xmax,ymin,ymax)


      iw = ifile+3

      call pyimage3(iw,nlat,nlat,errt,ncomp-1,xs(is),ys(is),nd,xmin,
     1       xmax,ymin,ymax)


      stop

      return
      end
c----------------------------------------------------------------------

c-----------------------------------------------------------------------      


      subroutine gethboundarydata(zk,k,zout,cs,nt,targ,pot,grad)
c     This subroutine computes the biharmic potential given by
c     \phi = \sum_{j=0}^{k} cs(j) (log (r_j) + h0(r_{j}) where
c     r_j = sqrt((xt - zout(1,j)**2 + (yt-zout(j))**2). The
c     subroutine returns the potential and the corresponding
c     gradient

      implicit real *8 (a-h,o-z)
      integer k,i,nt
      real *8 zout(2,1:k), cs(1:k), xt,yt,r,targ(2,*)
      complex *16 pot(*),grad(2,*),zk,ima,z,h0,h1

      data ima/(0.0d0,1.0d0)/


      do j = 1,nt
         xt = targ(1,j)
         yt = targ(2,j)
         
         pot(j) = 0.0d0
         grad(1,j) = 0.0d0
         grad(2,j) = 0.0d0
         ifexpon = 1
         do i=1,k
            r = sqrt((xt-zout(1,i))**2 + (yt-zout(2,i))**2)
            z = r*zk
            call hank103(z,h0,h1,ifexpon)

            pot(j) = pot(j) + cs(i)*(log(r) + h0) 
c     c         pot = pot + cs(i)*log(r) 
            grad(1,j) = grad(1,j) + 
     1           cs(i)*(xt-zout(1,i))*(1/r**2-h1*zk/r)
            grad(2,j) = grad(2,j) + 
     1           cs(i)*(yt-zout(2,i))*(1/r**2-h1*zk/r)
         enddo
      enddo
      return 
      end
c-------------------------------------------------------------------     
      subroutine wevaltargstoke(zk,nt,targ,wgeo,ncomp,nchs,cms,
     1     q1,q2,soln,wvals,wbdry,wintex)

      implicit real *8 (a-h,o-z)
      real *8 targ(2,*),wgeo(*),cms(2,*)
      integer ncomp, nchs(*),nt
      complex *16 soln(*),zk,wvals(*),q1,q2,wbdry(*),wintex
c     local
      complex *16, allocatable :: sneu(:,:),
     1     slay(:,:), b1(:), b2(:), b3(:), work(:),
     2     streammat(:,:)
      real *8, allocatable :: rnorms(:,:,:), whts(:)
      real *8 :: p1,p2,p3,p4,errs(1000)
      integer k, nch, ichunks, iadjs, iders, iders2, ihs
      integer npts, i, j
      complex *16 :: zero, one, oneint, ztemp, wint, zmatt(2,2), mu(2)
      complex *16 :: zval, zgrad(2), zhess(2,2)
      data zero, one / (0.0d0,0.0d0), (1.0d0,0.0d0) /

      external fgreenlap, zkernel_sprime, zkernel_slp, fgreensdummy
      external zhelmstokes_stream_kern

      ier = 0

c     get dimensions
      call chunkunpack1(wgeo,k,nch,ichunks,iadjs, 
     1     iders,iders2,ihs)

      npts = k*nch
      call prinf('k *',k,1)
      call prinf('nch *',nch,1)      

c     normals and smooth integration weights
      allocate(rnorms(2,k,nch),whts(k*nch))
      call chunknormals(wgeo,rnorms)
      call chunkwhts(k,nch,wgeo(ichunks),wgeo(iders), 
     1     wgeo(ihs),whts)

c     needed submatrices
      allocate(sneu(npts,npts),slay(npts,npts), 
     1     b1(npts),b2(npts))

      do i = 1,npts
         b1(i) = zero
         b2(i) = zero
      enddo

      do j = 1,npts
         do i = 1,npts
            sneu(i,j) = zero
            slay(i,j) = zero
         enddo
      enddo

c     build kernel matrices

c     neumann problem system matrix (0.5+S'+1*w^T)
      call zbuildmat(k,wgeo,zkernel_sprime,q1,q2,
     1     fgreenlap,zk,pars1,pars2,npts,sneu)
      do i = 1,npts
         sneu(i,i) = sneu(i,i) + 0.5d0
      enddo
      do j = 1,npts
         do i = 1,npts
            sneu(i,j) = sneu(i,j) + whts(j)
         enddo
      enddo

c     single layer evaluation matrix
      call zbuildmat(k,wgeo,zkernel_slp,q1,q2,
     1     fgreenlap,zk,pars1,pars2,npts,slay)

c     grab normal part of density
      ii = 1
      do i = 1,nch
         do j = 1,k
            b1(ii) = -q2*(soln(2*ii-1)*rnorms(1,j,i) +
     1           soln(2*ii)*rnorms(2,j,i))
            ii = ii+1
         enddo
      enddo

c     evaluate corresponding single layer potential
      call multa(slay,p1,p2,p3,p4,b1,b2,npts)
c     evaluate its tangential derivative
      ndim = 2
      call chunkderf(b2,b1,ndim,k,nch,wgeo(ichunks),
     1     wgeo(iders),wgeo(ihs))

      ztemp = zero
      do i = 1,npts
         ztemp = ztemp+whts(i)*b1(i)
      enddo
      call prin2('ztemp*',ztemp,2)

c     solve neumann problem with tau deriv of S as data

      ngmrec = 200
      numit = 200
      lw = (ngmrec*2 + 4)*npts
      allocate(work(lw))
      job = 0
      ier = 0
      eps = 1.0d-14

      call prinf('solving neumann problem *',ier,0)
      
      call cgmres(ier,npts,sneu,multa,p1,p2,p3,p4,b1,
     1     eps,numit,b2,niter,errs,ngmrec,work)

      call prin2('errs *',errs,niter)

      ztemp = zero
      do i = 1,npts
         ztemp = ztemp+whts(i)*b2(i)
      enddo

      call prin2('ztemp*',ztemp,2)

c     evaluate this potential back on the boundary 
c     this is the non-stream function part of a stokes
c     rep converted into a stream function 

      call multa(slay,p1,p2,p3,p4,b2,b1,npts)

c     evaluate stream function part of stokes on boundary
      npts2 = npts*2
      allocate(streammat(npts2,npts2),b3(npts2))
      ndim = 2
c     note the stream_kern function is a hack:
c     in order to use the vector type matrix builder,
c     every other row is zero ( the stream
c     function is scalar while the density is a vector)
      call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_stream_kern,
     1     q1,q2,fgreensdummy,zk,pars1,pars2,npts2,streammat)

      call multa(streammat,p1,p2,p3,p4,soln,b3,npts2)

      do i = 1,npts
         wbdry(i) = b1(i) + b3(2*i-1)
      enddo

c     get integral of w along boundary 

c      wint = zero
c      oneint = zero
c      do i = 1,npts
c         oneint = oneint + whts(i)*one
c         wint = wint + whts(i)*wbdry(i)
c      enddo

c      call prin2('wintex *',wintex,2)
c      call prin2('wint *',wint,2)
c      call prin2('oneint *',oneint,2)

c      ztemp = (wintex-wint)/oneint


c      call prin2('ztemp *',ztemp,2)
      
c      write(*,*) b1(1)+ztemp
c      write(*,*) b3(1)

      do i = 1,npts
         wbdry(i) = wbdry(i) + soln(npts2+1)
      enddo

      do i = 1,npts
         do iii = 2,ncomp
            isrc = ichunks+2*(i-1)
            call fgreenlap(zk,cms(1,iii),wgeo(isrc),
     1           pars1,pars2,zval,zgrad,zhess)
            wbdry(i) = wbdry(i) + zval*soln(npts2+iii)
         enddo
      enddo
      
c     get value at targets

      do i = 1,nt
         wvals(i) = soln(npts2+1)
      enddo

      do ii = 1,nt
         do iii = 2,ncomp
            call fgreenlap(zk,cms(1,iii),targ(1,ii),
     1           pars1,pars2,zval,zgrad,zhess)
            wvals(ii) = wvals(ii) + zval*soln(npts2+iii)
         enddo
         do i = 1,nch
            do j = 1,k
               isrc = ichunks + (i-1)*2*k+2*(j-1)
               iii = (i-1)*k+j
               call zhelmstokes_stream_kern(zk,wgeo(isrc),
     1              targ(1,ii),rnorms(1,j,i),rnorms(1,j,i),
     2              q1,q2,fgreensdummy,pars1,pars2,zmatt)
               call zkernel_slp(zk,wgeo(isrc),
     1              targ(1,ii),rnorms(1,j,i),rnorms(1,j,i),
     2              q1,q2,fgreenlap,pars1,pars2,zval)

               mu(1) = soln(2*iii-1)
               mu(2) = soln(2*iii)
               wvals(ii) = wvals(ii) + whts(iii)*(mu(1)*zmatt(1,1)+
     1              mu(2)*zmatt(1,2)+zval*b2(iii))
            enddo
         enddo
      enddo

      
      return
      end
c-----------------------------------------------------------------      

      subroutine multa(a,p1,p2,p3,p4,x,y,n)
      implicit real *8 (a-h,o-z)
c
cc       computes the product a*x = y
c
      complex *16 a(n,*),p1,p2,p3,p4,x(*),y(*)
      complex *16 alpha, beta

      beta = 0.0d0
      incx = 1
      incy = 1
      alpha = 1.0d0

      call zgemv ('N', n, n, alpha, a, n, x, incx, beta, y, incy)

      return
      end
c---------------------------------------------
      subroutine sqrtscalemat(ntot,n,xmat,qwts)
      implicit real *8 (a-h,o-z)
      real *8 qwts(*)
      complex *16 xmat(ntot,*)

      call prinf('ntot=*',ntot,1)
      call prinf('n=*',n,1)
      call prin2('qwts=*',qwts,12)

      do i=1,n
      do j=1,n

      xmat(i,j) = xmat(i,j)*sqrt(qwts(i)/qwts(j))
      xmat(i,j+n) = xmat(i,j+n)*sqrt(qwts(i)/qwts(j))
      xmat(i+n,j) = xmat(i+n,j)*sqrt(qwts(i)/qwts(j))
      xmat(i+n,j+n) = xmat(i+n,j+n)*sqrt(qwts(i)/qwts(j))

      enddo
      enddo

      return
      end
c-------------------------------------------------------------      
      subroutine genlatsrctarg(ncol,nsrc,ntarg,r,a,b,pert,
     1    targ,cms,ck,xy)
      implicit real *8 (a-h,o-z)
      real *8 targ(2,*),cms(2,*),xy(2,*),ck(*),xstart(10000)

      nrow = 2
      nsrc = nrow*ncol


      xsep = 5*r
      xtot = xsep*(ncol-1)+2*r
      xleft = a - xtot

      xstart(1) = xleft/2 + 0.02 + r
      xstart(2) = xleft/2 - 0.02 + r

      ysep = 5*r
      ytot = ysep*(nrow-1)+2*r
      yleft = b - ytot


      ystart = yleft/2 + r


      icomp = 2
      cms(1,1) = 3*a
      cms(2,1) = 2*b
      ck(1) = hkrand(0)
      
      ntarg = 0


      
      do 1100 i=1,nrow
      do 1000 j=1,ncol

      cms(1,icomp) = xstart(i) + (j-1)*xsep
      cms(2,icomp) = ystart + (i-1)*ysep

      if(i.eq.1.and.j.ne.ncol) then

      ntarg = ntarg+1
      targ(1,ntarg) = cms(1,icomp) + 2.5*r + pert*(hkrand(0)-0.5)*r
      targ(2,ntarg) = cms(2,icomp) + 2.5*r + pert*(hkrand(0)-0.5)*r

      ntarg = ntarg+1
      targ(1,ntarg) = cms(1,icomp) + 2.5*r + pert*(hkrand(0)-0.5)*r
      targ(2,ntarg) = cms(2,icomp) - 2.5*r + pert*(hkrand(0)-0.5)*r

      endif

      if(i.eq.2.and.j.ne.ncol) then

      ntarg = ntarg+1
      targ(1,ntarg) = cms(1,icomp) + 2.5*r + pert*(hkrand(0)-0.5)*r
      targ(2,ntarg) = cms(2,icomp) + 2.5*r + pert*(hkrand(0)-0.5)*r

      endif
      
      xy(1,icomp-1) = cms(1,icomp) + pert*(hkrand(0)-0.5)*r
      xy(2,icomp-1) = cms(2,icomp) + pert*(hkrand(0)-0.5)*r

      ck(icomp-1) = hkrand(0)
      icomp = icomp + 1

 1000 continue      
 1100 continue      

      return
      end
c----------------------------------------------------------------      
      subroutine getchunks(iort,k,nch,chunks,adjs,ders,ders2,hs)
      implicit real *8 (a-h,o-z)
      real *8 chunks(2,k,*),ders(2,k,*),ders2(2,k,*),hs(*)
      integer adjs(2,*)
      real *8 ts(1000),uk(10000),vk(100000),whts(1000)

      done = 1
      pi = atan(done)*4

      h = 2*pi/nch
      itype = 0
      call legeexps(itype,k,ts,uk,vk,whts)

      do 1100 i=1,nch

      tstart = (i-1)*h
      tend = i*h

      hs(i) = (tend-tstart)/2
      adjs(1,i) = i-1
      adjs(2,i) = i+1

      do 1000 j=1,k

      theta = tstart + (ts(j)+1)/2*(tend-tstart)
      if(iort.eq.-1) theta = 2*pi - theta
      chunks(1,j,i) = cos(theta)
      chunks(2,j,i) = sin(theta)

      ders(1,j,i) = -sin(theta)
      ders(2,j,i) = cos(theta)

      ders2(1,j,i) = -cos(theta)
      ders2(2,j,i) = -sin(theta)

 1000 continue      
 1100 continue      

      adjs(1,1) = nch
      adjs(2,nch) = 1

      return
      end
c-----------------------------------------------------------------------      

      subroutine wevaltargstoke2(zk,nt,targ,isint,wgeo,ncomp,nchs,cms,
     1     q1,q2,soln,wvals,wbdry,wintex)

      implicit real *8 (a-h,o-z)
      real *8 targ(2,*),wgeo(*),cms(2,*)
      integer ncomp, nchs(*),nt,isint(*)
      complex *16 soln(*),zk,wvals(*),q1,q2,wbdry(*),wintex
c     local
      complex *16, allocatable :: sneu(:,:),
     1     slay(:,:), b1(:), b2(:), b3(:), work(:),
     2     streammat(:,:)
      real *8, allocatable :: rnorms(:,:,:), whts(:)
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),ders2(:,:,:)
      real *8, allocatable :: hs(:),dsdt(:,:),rnx(:),rny(:)
      integer, allocatable :: adjs(:,:)

      real *8, allocatable :: xcoefs(:,:),ycoefs(:,:),rnxc(:,:)
      real *8, allocatable :: rnyc(:,:),dsdtc(:,:)
      complex *16, allocatable :: soln1c(:,:),soln2c(:,:),soln3c(:,:)
      complex *16, allocatable :: soln1(:),soln2(:)
      real *8, allocatable :: ts(:),uk(:,:),vk(:,:),whts0(:)

      real *8 :: p1,p2,p3,p4,errs(1000)
      integer k, nch, ichunks, iadjs, iders, iders2, ihs
      integer npts, i, j
      complex *16 :: zero, one, oneint, ztemp, wint, zmatt(2,2), mu(2)
      complex *16 :: zval, zgrad(2), zhess(2,2),pottmp,pottmp2
      data zero, one / (0.0d0,0.0d0), (1.0d0,0.0d0) /

      external fgreenlap, zkernel_sprime, zkernel_slp, fgreensdummy
      external zhelmstokes_stream_kern

      ier = 0

c     get dimensions
      call chunkunpack1(wgeo,k,nch,ichunks,iadjs, 
     1     iders,iders2,ihs)

      npts = k*nch
      call prinf('k *',k,1)
      call prinf('nch *',nch,1)      

c     normals and smooth integration weights
      allocate(rnorms(2,k,nch),whts(k*nch))
      call chunknormals(wgeo,rnorms)
      call chunkwhts(k,nch,wgeo(ichunks),wgeo(iders), 
     1     wgeo(ihs),whts)

c     needed submatrices
      allocate(sneu(npts,npts),slay(npts,npts), 
     1     b1(npts),b2(npts))

      do i = 1,npts
         b1(i) = zero
         b2(i) = zero
      enddo

      do j = 1,npts
         do i = 1,npts
            sneu(i,j) = zero
            slay(i,j) = zero
         enddo
      enddo

c     build kernel matrices

c     neumann problem system matrix (0.5+S'+1*w^T)
      call zbuildmat(k,wgeo,zkernel_sprime,q1,q2,
     1     fgreenlap,zk,pars1,pars2,npts,sneu)
      do i = 1,npts
         sneu(i,i) = sneu(i,i) + 0.5d0
      enddo
      do j = 1,npts
         do i = 1,npts
            sneu(i,j) = sneu(i,j) + whts(j)
         enddo
      enddo

c     single layer evaluation matrix
      call zbuildmat(k,wgeo,zkernel_slp,q1,q2,
     1     fgreenlap,zk,pars1,pars2,npts,slay)

c     grab normal part of density
      ii = 1
      do i = 1,nch
         do j = 1,k
            b1(ii) = -q2*(soln(2*ii-1)*rnorms(1,j,i) +
     1           soln(2*ii)*rnorms(2,j,i))
            ii = ii+1
         enddo
      enddo

c     evaluate corresponding single layer potential
      call multa(slay,p1,p2,p3,p4,b1,b2,npts)
c     evaluate its tangential derivative
      ndim = 2
      call chunkderf(b2,b1,ndim,k,nch,wgeo(ichunks),
     1     wgeo(iders),wgeo(ihs))

      ztemp = zero
      do i = 1,npts
         ztemp = ztemp+whts(i)*b1(i)
      enddo
      call prin2('ztemp*',ztemp,2)

c     solve neumann problem with tau deriv of S as data

      ngmrec = 200
      numit = 200
      lw = (ngmrec*2 + 4)*npts
      allocate(work(lw))
      job = 0
      ier = 0
      eps = 1.0d-14

      call prinf('solving neumann problem *',ier,0)
      
      call cgmres(ier,npts,sneu,multa,p1,p2,p3,p4,b1,
     1     eps,numit,b2,niter,errs,ngmrec,work)

      call prin2('errs *',errs,niter)

      ztemp = zero
      do i = 1,npts
         ztemp = ztemp+whts(i)*b2(i)
      enddo

      call prin2('ztemp*',ztemp,2)

c     evaluate this potential back on the boundary 
c     this is the non-stream function part of a stokes
c     rep converted into a stream function 

      call multa(slay,p1,p2,p3,p4,b2,b1,npts)

c     evaluate stream function part of stokes on boundary
      npts2 = npts*2
      allocate(streammat(npts2,npts2),b3(npts2))
      ndim = 2
c     note the stream_kern function is a hack:
c     in order to use the vector type matrix builder,
c     every other row is zero ( the stream
c     function is scalar while the density is a vector)
      call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_stream_kern,
     1     q1,q2,fgreensdummy,zk,pars1,pars2,npts2,streammat)

      call multa(streammat,p1,p2,p3,p4,soln,b3,npts2)

      do i = 1,npts
         wbdry(i) = b1(i) + b3(2*i-1)
      enddo

c     get integral of w along boundary 

c      wint = zero
c      oneint = zero
c      do i = 1,npts
c         oneint = oneint + whts(i)*one
c         wint = wint + whts(i)*wbdry(i)
c      enddo

c      call prin2('wintex *',wintex,2)
c      call prin2('wint *',wint,2)
c      call prin2('oneint *',oneint,2)

c      ztemp = (wintex-wint)/oneint


c      call prin2('ztemp *',ztemp,2)
      
c      write(*,*) b1(1)+ztemp
c      write(*,*) b3(1)

      do i = 1,npts
         wbdry(i) = wbdry(i) + soln(npts2+1)
      enddo

      do i = 1,npts
         do iii = 2,ncomp
            isrc = ichunks+2*(i-1)
            call fgreenlap(zk,cms(1,iii),wgeo(isrc),
     1           pars1,pars2,zval,zgrad,zhess)
            wbdry(i) = wbdry(i) + zval*soln(npts2+iii)
         enddo
      enddo
c
cc       setup everuything for adaptive integration
c
      nn = k*nch
      allocate(chunks(2,k,nch),ders(2,k,nch),ders2(2,k,nch),hs(nch))
      allocate(dsdt(k,nch),rnx(nn),rny(nn),adjs(2,nch))
      allocate(soln1(nn),soln2(nn))

      allocate(xcoefs(k,nch),ycoefs(k,nch),rnxc(k,nch),rnyc(k,nch))
      allocate(dsdtc(k,nch),soln1c(k,nch),soln2c(k,nch))
      allocate(soln3c(k,nch))

      call chunkunpack(wgeo,k,nch,chunks,adjs,ders,ders2,hs)

      allocate(ts(k),uk(k,k),vk(k,k),whts0(k))

      itype = 2
      call legeexps(itype,k,ts,uk,vk,whts0)

      ii = 1
      do ich=1,nch
      do j=1,k

      dsdt(j,ich) = sqrt(ders(1,j,ich)**2 + ders(2,j,ich)**2)
      rnx(ii) = rnorms(1,j,ich)
      rny(ii) = rnorms(2,j,ich)

      ii = ii+1

      enddo
      enddo


      call prinf('nn=*',nn,1)
      do i=1,nn

      soln1(i) = soln(2*i-1)
      soln2(i) = soln(2*i)

      enddo


      do ich=1,nch

      istart = (ich-1)*k+1
      call chunksexps_fast(k,uk,chunks(1,1,ich),xcoefs(1,ich),
     1       ycoefs(1,ich))

      call matvec(k,uk,rnx(istart),rnxc(1,ich))
      call matvec(k,uk,rny(istart),rnyc(1,ich))
      call matvec(k,uk,dsdt(1,ich),dsdtc(1,ich))

      call matvecc(k,uk,soln1(istart),soln1c(1,ich))
      call matvecc(k,uk,soln2(istart),soln2c(1,ich))
      call matvecc(k,uk,b2(istart),soln3c(1,ich))

      enddo
      

cc         call prin2('xcoefs=*',xcoefs,k)
cc         call prin2('ycoefs=*',ycoefs,k)
cc         call prin2_long('rnx=*',rnx,k)
cc         call prin2_long('rny=*',rny,k)
cc         call prin2('chunks=*',chunks,2*k)
cc         call prin2('rnxc=*',rnxc,k)
cc         call prin2('rnyc=*',rnyc,k)
cc         call prin2('dsdtc=*',dsdtc,k)
cc         call prin2('soln1c=*',soln1c,2*k)
cc         call prin2('soln2c=*',soln2c,2*k)
cc         call prin2('soln3c=*',soln3c,2*k)
cc         call prin2('soln1=*',soln1,2*k)
cc         call prin2('soln2=*',soln2,2*k)
cc         call prin2('b2=*',b2,2*k)

c     get value at targets


      do ii = 1,nt
         if(isint(ii).eq.1) then

         wvals(ii) = soln(npts2+1)
         do iii = 2,ncomp
            call fgreenlap(zk,cms(1,iii),targ(1,ii),
     1           pars1,pars2,zval,zgrad,zhess)
            wvals(ii) = wvals(ii) + zval*soln(npts2+iii)
         enddo
         pottmp = 0
         call potevaladap(targ(1,ii),targ(2,ii),k,nch,xcoefs,ycoefs,
     1           rnxc,rnyc,dsdtc,hs,q1,q2,zk,soln1c,soln2c,soln3c,
     2           pottmp)

         wvals(ii) = wvals(ii) + pottmp

         endif
      enddo

      
      return
      end
c---------------------------------------
       subroutine potevaladap(targx,targy,k,nch,xcoefs,ycoefs,
     1           rnxc,rnyc,dsdtc,hs,q1,q2,zk,soln1c,soln2c,soln3c,
     2           pot)

       implicit real *8 (a-h,o-z)
c
cc       this subroutine computes the modified biharmonic 
c        stream function using adaptive gaussian integration 
c       
       real *8 xcoefs(k,*),ycoefs(k,*),rnxc(k,*),rnyc(k,*),dsdtc(k,*)
       real *8 hs(*)
       complex *16 q1,q2,zk,soln1c(k,*),soln2c(k,*),soln3c(k,*),pot
       complex *16 cint
       real *8 pars(100000),pars2(100000)
       external fint

       pot = 0
       a = -1.0d0
       b = 1.0d0
       m = 16
       eps = 1.0d-8

       pars2(1) = k +0.1d0
       pars2(2) = targx
       pars2(3) = targy
       pars2(4) = real(zk)
       pars2(5) = imag(zk)
       pars2(6) = real(q1)
       pars2(7) = imag(q1)
       pars2(8) = real(q2)
       pars2(9) = imag(q2)

       do i=1,nch

       do j=1,k

       pars(j) = xcoefs(j,i)
       pars(j+k) = ycoefs(j,i)
       pars(j+2*k) = rnxc(j,i)
       pars(j+3*k) = rnyc(j,i)
       pars(j+4*k) = dsdtc(j,i)
       pars(j+5*k) = real(soln1c(j,i))
       pars(j+6*k) = imag(soln1c(j,i))
       pars(j+7*k) = real(soln2c(j,i))
       pars(j+8*k) = imag(soln2c(j,i))
       pars(j+9*k) = real(soln3c(j,i))
       pars(j+10*k) = imag(soln3c(j,i))

       enddo

       ier = 0
       cint = 0
       call cadapgau(ier,a,b,fint,pars,pars2,m,eps,cint,maxrec,numint)

       if(ier.ne.0) then
       call prinf('i=*',i,1)
       call prin2('pars2=*',pars2,3)

       stop

       endif

       pot = pot + cint*hs(i)

       enddo

       return
       end
c-----------------------------------------------------------------    

      function fint(t,pars,pars2)
      implicit real *8 (a-h,o-z)
      real *8 pars(*),pars2(*),targ(2),rnorm(2),src(2)
      complex *16 fint,soln1,soln2,soln3,zk,q1,q2,zmatt(2,2),zval

      external fgreenlap, fgreensdummy

      k = pars2(1)
      targ(1) = pars2(2)
      targ(2) = pars2(3)
      zkr = pars2(4)
      zki = pars2(5)
      q1r = pars2(6)
      q1i = pars2(7)
      q2r = pars2(8)
      q2i = pars2(9)

      zk = cmplx(zkr,zki)
      q1 = cmplx(q1r,q1i)
      q2 = cmplx(q2r,q2i)

      call legeexev(t,src(1),pars(1),k-1)
      call legeexev(t,src(2),pars(k+1),k-1)
      call legeexev(t,rnorm(1),pars(2*k+1),k-1)
      call legeexev(t,rnorm(2),pars(3*k+1),k-1)
      call legeexev(t,dsdt,pars(4*k+1),k-1)
      call legeexev(t,soln1r,pars(5*k+1),k-1)
      call legeexev(t,soln1i,pars(6*k+1),k-1)
      call legeexev(t,soln2r,pars(7*k+1),k-1)
      call legeexev(t,soln2i,pars(8*k+1),k-1)
      call legeexev(t,soln3r,pars(9*k+1),k-1)
      call legeexev(t,soln3i,pars(10*k+1),k-1)


      soln1 = cmplx(soln1r,soln1i)
      soln2 = cmplx(soln2r,soln2i)
      soln3 = cmplx(soln3r,soln3i)

      call zhelmstokes_stream_kern(zk,src,targ,rnorm,rnorm,q1,q2,
     1       fgreensdummy,p1,p2,zmatt)

      call zkernel_slp(zk,src,targ,rnorm,rnorm,q1,q2,fgreenlap,
     1         p1,p2,zval)


       fint = (soln1*zmatt(1,1) + soln2*zmatt(1,2) + soln3*zval)*dsdt


      return
      end

c-------------------------------------------------      

      subroutine matvec(n,a,x,y)
      implicit real *8 (a-h,o-z)
      real *8 a(n,n),x(*),y(*)

      do i=1,n

      y(i) = 0

      do j=1,n

      y(i)  = y(i) + a(i,j)*x(j)

      enddo
      enddo

      return
      end
c-------------------------------------------------      

      subroutine matvecc(n,a,x,y)
      implicit real *8 (a-h,o-z)
      real *8 a(n,n)
      complex *16 x(*),y(*)

      do i=1,n

      y(i) = 0

      do j=1,n

      y(i)  = y(i) + a(i,j)*x(j)

      enddo
      enddo

      return
      end
