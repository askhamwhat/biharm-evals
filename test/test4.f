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
      call biharmdmatgen(iw2,aa,bb,ff,ntot,w(ising))
      

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
      subroutine biharmdmatgen(ifile,aa,bb,ff,ntot,sing)

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
      complex *16, allocatable :: xmat(:,:),xmat2(:,:)
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
      complex *16 gradtmp(2)
      
      complex *16 sing(*)

      integer ncomp, nwiggles
      real *8 tmp, errp, rp
      real *8 errs(1000)

      complex *16 wex,wval,eye

      data eye/(0.0d0,1.0d0)/

      external multa


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
c     Get geometry from file
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
      allocate(xmat(ntot,ntot),xmat2(ntot,ntot))


      call prinf('k=*',k,1)
      call prinf('ntot=*',ntot,1)
      call prinf('n=*',n,1)

      call formhbiharmdmatfark(k,ntot,n,zk,wgeos,xmat,rkappa,qwts)

      call prin2('xmat=*',xmat,24)
      call prin2('rkappa=*',rkappa,12)

c
cc
cc      copy xmat
c 
       do i=1,ntot
       do j=1,ntot

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
             rhs(ii) = pot(j,i)*sqrt(qwts(ii))
             rhs(ii+n) = potn(j,i)*sqrt(qwts(ii))
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

      do i=1,n

      soln(i) = soln(i)/sqrt(qwts(i))
      soln(i+n) = soln(i+n)/sqrt(qwts(i))

      enddo


      call prin2('soln=*',soln,24)
      

      wval = 0
      wex = 0
      targ(1) = 0.7d0
      targ(2) = 0.27d0
      call wevaltargfark(zk,targ,xs,ys,dxdt,dydt,dsdt,rnx,rny,rkappa,
     1        qwts,soln,n,ntot,wval)

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
      call csvdpiv2(ier,xmat2,ntot,ntot,sing,ncols,eps,workspace,
     1        lw,ltot)

      call prinf('ier=*',ier,1)
      call prinf('ncols=*',ncols,1)


      return
      end
c----------------------------------------------------------------------

      subroutine formhbiharmdmatfark(k,n,ns,zk,wgeo,xmat,rkappa,qwts)
      implicit real *8 (a-h,o-z)
      real *8 wgeo(*),rkappa(*),qwts(*)
      complex *16 xmat(n,*)

      complex *16 p1(10),p2,zk,par0,pars1,pars2

c     matrix formation, eigenvalue calculation

      complex *16, allocatable :: tempmat(:,:)
      complex *16 q1,q2 
      
      external zfark_kernel, zfark_kernelp
      external zhfark_ck1,zhfark_ck2 

      done = 1

      allocate(tempmat(ns,ns))

      call prinf('k=*',k,1)
      call prinf('n=*',n,1)
      call prinf('ns=*',ns,1)

      do i=1,n
      do j=1,n

      xmat(i,j) = 0

      enddo
      enddo

      do i=1,ns
      do j=1,ns

      tempmat(i,j) = 0

      enddo
      enddo


      q1 = (0,0)
      q2 = (0,0)

      par0 = zk

      call zbuildmat(k, wgeo, zfark_kernel, q1, q2,
     1     zhfark_ck1, par0, pars1, pars2, ntot, tempmat)
      
      xmat(1:ns,1:ns) = tempmat
      

      call zbuildmat(k, wgeo, zfark_kernel, q1, q2,
     1     zhfark_ck2, par0, pars1, pars2, ntot, tempmat)

      xmat(1:ns,ns+1:n) = tempmat
      
      call zbuildmat(k, wgeo, zfark_kernelp, q1, q2,
     1     zhfark_ck1, par0, pars1, pars2, ntot, tempmat)

      xmat(ns+1:n,1:ns) = tempmat

      call zbuildmat(k, wgeo, zfark_kernelp, q1, q2,
     1     zhfark_ck2, par0, pars1, pars2, ntot, tempmat)

      call prin2('par0=*',par0,2)

      xmat(ns+1:n,ns+1:n) = tempmat

      do i=1,n

      xmat(i,i) = xmat(i,i) + 0.5d0

      enddo

      do i=1,ns

      xmat(i+ns,i) = xmat(i+ns,i)-rkappa(i)

      enddo

      call sqrtscalemat(n,ns,xmat,qwts)

      return
      end

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
      subroutine wevaltargfark(zk,targ,xs,ys,dxdt,dydt,dsdt,rnx,rny,
     1        rkappa,qwts,soln,n,ntot,wval)

      implicit real *8 (a-h,o-z)
      real *8 targ(2),xs(*),ys(*),dxdt(*),dydt(*),dsdt(*),rnx(*),
     1          rny(*),rkappa(*),qwts(*),src(2),rnsrc(2)
      complex *16 soln(*),wval,zk,z,val,val2,grad(2),hess(4)

      wval = 0

      done = 1
      pi = atan(done)*4


      do i=1,n
      
      src(1) = xs(i)
      src(2) = ys(i)
      rnsrc(1) = rnx(i)
      rnsrc(2) = rny(i)


      call zhfark_ck1(zk,src,targ,rnsrc,pars2,val,grad,hess)
      call zhfark_ck2(zk,src,targ,rnsrc,pars2,val2,grad,hess)

      wval = wval + val*qwts(i)*soln(i)
      wval = wval + val2*qwts(i)*soln(i+n)

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


      do i=1,n

      y(i) = 0

      do j=1,n

      y(i) = y(i) + a(i,j)*x(j)

      enddo
      enddo

      
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
