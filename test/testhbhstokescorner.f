      program main
      implicit real *8 (a-h,o-z)
      real *8 w(1000000)

      call prini(6,13)
      call prin2('Enter icase*',n,0)
      read *, icase


      aa = 1
      bb = 0.5

 2100 format('data for ff=',e11.5)
 2110 format('npts=',i7)
 2120 format('sing=')
 2130 format(2x,e11.5)
 2150 format('c  ')

      ff = 2.5*2**(icase-1.0d0)
      
      ising = 1

      iw2 = 26+icase
      call testhbhstokescorner(iw2,aa,bb,ff,ntot,w(ising))
      

      iw = 16 + icase
      write(iw,2100) ff
      write(iw,2110) ntot
      write(iw,2150)
      write(iw,2120)
      do i=1,ntot

      write(iw,2130) w(ising+2*i-2)

      enddo

      stop
      end

c------------------------------------------------------     
c
      subroutine testhbhstokescorner(ifile,aa,bb,ff,ntot,sing)

      implicit real *8 (a-h,o-z)

      real *8, allocatable :: chunks(:,:,:)
      real *8, allocatable :: ders(:,:,:)
      real *8, allocatable :: rn(:,:,:)
      real *8 sc
      real *8, allocatable :: ders2(:,:,:),dsdta(:,:)
      real *8, allocatable :: hs(:)
      integer, allocatable :: adjs(:,:)
      real *8, allocatable :: ck(:)

      real *8, allocatable :: wgeos(:)
      real *8, allocatable :: wgeos3(:)

      complex *16, allocatable :: work(:),workspace(:)
      complex *16, allocatable :: rhs(:),soln(:),soln2(:)
      integer, allocatable :: nchs(:)

      complex *16 zk

c     matrix formation, eigenvalue calculation

      real *8, allocatable :: whts(:,:)
      complex *16, allocatable :: xmat(:,:),xmat2(:,:), onesmat(:,:)
      complex *16 p1(10),p2,p3,p4
      real *8, allocatable :: cms(:,:)
      
      real *8, allocatable :: xs(:),ys(:),cxs(:),cys(:),qwts(:)
      real *8, allocatable :: dxdt(:),dydt(:),dsdt(:),rnx(:),rny(:)
      real *8, allocatable :: rkappa(:)

      integer novers(1000), iparsindeces(1000)
      integer ifcloseds(1000)
      real *8 chsmalls(1000), tas(1000), tbs(1000)
      real *8 epss(1000), pars(10 000), pars1(100),pars2(100)

      real *8 src(2),targ(2),verts(2,10000),widths(10000)
      
      complex *16, allocatable :: pot(:,:)
      complex *16, allocatable :: pottau(:,:)
      complex *16, allocatable :: grad(:,:,:)
      complex *16, allocatable :: potn(:,:)
      complex *16, allocatable :: wbdry(:,:)
      complex *16 gradtmp(2), wintex
      
      complex *16 sing(*), q1, q2

      integer ncomp, nwiggles
      real *8 tmp, errp, rp
      real *8 errs(1000)

      complex *16 wex,wval,eye

      data eye/(0.0d0,1.0d0)/

      external multa, zhelmstokes_kern, fgreensdummy


      done = 1
      pi = atan(done)*4
      nchunkmax = 10000

      zk = 1.2d0

      eps = 1.0d-13
      ifclosed = 1
      chsmall = 1.0d-5
      nover = 2
      k = 16

      ncomp = 1
      allocate(chunks(2,k,nchunkmax))
      allocate(ders(2,k,nchunkmax),dsdta(k,nchunkmax))
      allocate(rn(2,k,nchunkmax))
      allocate(ders2(2,k,nchunkmax))
      allocate(hs(nchunkmax))
      allocate(adjs(2,nchunkmax))

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set up geometry
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      lwgeos = 1000000

      allocate(wgeos3(lwgeos),wgeos(lwgeos))

      nverts = 5

      verts(1,1) = 0
      verts(2,1) = 0
      widths(1) = bb/ff

      verts(1,2) = aa
      verts(2,2) = 0
      widths(2) = bb/ff

      verts(1,3) = aa
      verts(2,3) = bb
      widths(3) = bb/ff

      verts(1,4) = 0
      verts(2,4) = bb
      widths(4) = bb/ff 

      verts(1,5) = 0
      verts(2,5) = 0
      widths(5) = bb/ff

      p1 = 0
      p2 = 0

      i1 = 0
      i2 = 0

      chsmall = 10
      ta = 0
      tb = 2*pi
      nover = 2
      pars(1) = 1.0
      pars(2) = 0.1d0

      ifbell = 1

      eps = 1.0d-6
      ifclosed = 1

      nch = 0

      call chunkpolysmooth(ier,eps,widths,ifbell,p1,p2,i1,i2,nverts,
     1      verts,ifclosed,nover,k,nch,chunks,adjs,ders,ders2,hs) 

      pert = 1.0d-2

      allocate(cms(2,ncomp),ck(ncomp))
      cms(1,1) = 10*aa/2 + hkrand(0)*pert
      cms(2,1) = 3*bb/2 + hkrand(0)*pert

      ck(1) = hkrand(0)



c     Allocate potential, gradient and potn arrays

      allocate(nchs(ncomp))

      nchs(1) = nch

      irefinelev= 0
      do ii=1,irefinelev
            nch1 = nch
            do i=1,nch1
               call chunksplit1(i,k,nch,chunks,adjs,ders,ders2,hs)
            enddo
            do i=1,ncomp

            nchs(i) = 2*nchs(i)

            enddo
        enddo

      call chunkpack(k,nch,chunks,adjs,ders,ders2,hs,wgeos,lused)
      n = k*nch
      allocate(xs(n),ys(n))
      allocate(dsdt(n),rkappa(n),dxdt(n),dydt(n),rnx(n),rny(n))
      allocate(whts(k,nch),qwts(n))

      call chunkunpack1(wgeos,k,nch,ichunks,iadjs,iders,iders2,ihs)
      
      call chunkwhts(k,nch,wgeos(ichunks),wgeos(iders),wgeos(ihs),whts)

      
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

      call pyplot(ifile,xs,ys,n,2,'a*')


      ntot = 2*k*nch
      allocate(xmat(ntot,ntot),xmat2(ntot,ntot),
     1     onesmat(ntot,ntot))


      call prinf('k=*',k,1)
      call prinf('ntot=*',ntot,1)
      call prinf('n=*',n,1)
      q1 = (0.0d0,0.0d0)
      q2 = (1.0d0,0.0d0)
      ndim = 2
      call zbuildmat_vec(ndim,k,wgeos,zhelmstokes_kern,
     1     q1,q2,fgreensdummy,zk,pars1,pars2,ntot,
     2     xmat)

      call normalonesmat(onesmat,whts,rn,k,nch)
      do j = 1,ntot
         do i = 1,ntot
            xmat(i,j) = xmat(i,j) + onesmat(i,j)
            if (i .eq. j) xmat(i,j) = xmat(i,j)-0.5d0
            xmat2(i,j) = xmat(i,j)
         enddo
      enddo
      

c
c     Allocate potential, gradient and potn arrays
c
      allocate(pot(k,nch),potn(k,nch),grad(2,k,nch),pottau(k,nch))
      allocate(rhs(ntot),soln(ntot),soln2(ntot2))
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
             call gethboundarydata(zk,ncomp,cms,ck,chunks(1,j,i),
     1                            chunks(2,j,i),pot(j,i),gradtmp)
             potn(j,i) = gradtmp(1)*rn(1,j,i)+gradtmp(2)*rn(2,j,i)
             ii = ii + 1
             rhs(2*ii-1) = -gradtmp(2)
             rhs(2*ii) = gradtmp(1)
         enddo
      enddo

      call prin2('rhs=*',rhs,24)

      
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
      eps = 1.0d-12

      call cgmres(ier,ntot,xmat,multa,p1,p2,p3,p4,rhs,eps,numit,soln,
     1        niter,errs,ngmrec,work)

      call prinf('ier=*',ier,1)
      call prinf('niter=*',niter,1)
      call prin2('errs=*',errs,niter)


      call prin2('soln=*',soln,24)

c     get integral of exact solution on outer boundary

      wintex = 0
      do i = 1,nch
         do j = 1,k
            wintex = wintex + whts(j,i)*pot(j,i)
         enddo
      enddo

      call prin2('wintex *',wintex,2)
      
      wval = 0
      wex = 0
      targ(1) = 0.7d0
      targ(2) = 0.27d0
      nt = 1

      allocate(wbdry(k,nch))

      call wevaltargstoke(zk,nt,targ,wgeos,ncomp,nchs,cms,
     1     q1,q2,soln,wval,wbdry,wintex)

      call prin2('wbdry*',wbdry,12)
      call prin2('pot*',pot,12)

      call prin2('wval=*',wval,2)
      call gethboundarydata(zk,ncomp,cms,ck,targ(1),targ(2),wex,gradtmp)

      call prin2('wex=*',wex,2)

      erra = abs(wex-wval)/abs(wex)

      call prin2('erra=*',erra,1)

      call prinf('ntot=*',ntot,1)

      ier = 0
      eps = 1.0d-12
      ncols = 0
      ltot = 0

      lw = ntot*ntot*10
      allocate(workspace(lw))
      call prinf('computing svd ...*',ier,0)
      time1 = second()
c      call csvdpiv2(ier,xmat2,ntot,ntot,sing,ncols,eps,workspace,
c     1     lw,ltot)
      time2 = second()      
      call prin2('done. time=*',time2-time1,1)      

      call prinf('ier=*',ier,1)
      call prinf('ncols=*',ncols,1)


      return
      end
c----------------------------------------------------------------------

c-----------------------------------------------------------------------      


      subroutine gethboundarydata(zk,k,zout,cs,xt,yt,pot,grad)
c     This subroutine computes the biharmic potential given by
c     \phi = \sum_{j=0}^{k} cs(j) (log (r_j) + h0(r_{j}) where
c     r_j = sqrt((xt - zout(1,j)**2 + (yt-zout(j))**2). The
c     subroutine returns the potential and the corresponding
c     gradient

      implicit real *8 (a-h,o-z)
      integer k,i
      real *8 zout(2,1:k), cs(1:k), xt,yt,r
      complex *16 pot,grad(2),zk,ima,z,h0,h1

      data ima/(0.0d0,1.0d0)/


      
      pot = 0.0d0
      grad(1) = 0.0d0
      grad(2) = 0.0d0
      ifexpon = 1
      do i=1,k
         r = sqrt((xt-zout(1,i))**2 + (yt-zout(2,i))**2)
         z = r*zk
         call hank103(z,h0,h1,ifexpon)

         pot = pot + cs(i)*(log(r) + h0) 
cc         pot = pot + cs(i)*log(r) 
         grad(1) = grad(1) + 
     1              cs(i)*(xt-zout(1,i))*(1/r**2-h1*zk/r)
         grad(2) = grad(2) + 
     1              cs(i)*(yt-zout(2,i))*(1/r**2-h1*zk/r)
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
      complex *16 :: zval
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
            b1(ii) = -(soln(2*ii-1)*rnorms(1,j,i) +
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
      ntot = npts*2
      allocate(streammat(ntot,ntot),b3(ntot))
      ndim = 2
c     note the stream_kern function is a hack:
c     in order to use the vector type matrix builder,
c     every other entry is zero ( the stream
c     function is scalar while the density is a vector)
      call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_stream_kern,
     1     q1,q2,fgreensdummy,zk,pars1,pars2,ntot,streammat)

      call multa(streammat,p1,p2,p3,p4,soln,b3,ntot)

      do i = 1,npts
         wbdry(i) = b1(i) + b3(2*i-1)
      enddo

c     get integral of w along boundary 

      wint = zero
      oneint = zero
      do i = 1,npts
         oneint = oneint + whts(i)*one
         wint = wint + whts(i)*wbdry(i)
      enddo

      call prin2('wintex *',wintex,2)
      call prin2('wint *',wint,2)
      call prin2('oneint *',oneint,2)

      ztemp = (wintex-wint)/oneint
      write(*,*) b1(1)+ztemp
      write(*,*) b3(1)

      do i = 1,npts
         wbdry(i) = wbdry(i) + ztemp
      enddo

c     get value at targets

      do i = 1,nt
         wvals(i) = ztemp
      enddo
      
      do ii = 1,nt
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
