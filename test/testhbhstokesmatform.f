
      program testhbhstokesmatform
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This program tests the matrix forming routines
c     in cmat_build_vec.f90 as applied to the Helmholtz
c     Stokes double layer potential
c
      
      implicit real *8 (a-h,o-z)
      parameter (nmax = 100000)
      integer adjs(2,nmax)
      real *8 hs(nmax),xyin(10),xyout(10)
c     
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),
     1     ders2(:,:,:),dsdt(:,:),rnorms(:,:,:),
     2     whts(:,:),taus(:,:,:),wgeo(:)
      complex *16, allocatable :: smu(:,:,:), dmu(:,:,:)
      complex *16, allocatable :: rhs(:,:,:)      

      complex *16, allocatable :: sysmat(:,:), dmat(:,:),
     1     onesmat(:,:)
      
      complex *16 :: veltemp(2), stresstemp(2,2), muexact(2)
      complex *16 :: veldertemp(2,2)
      complex *16 :: velsum(2), velexact(2), velsums(2), velsumd(2)
      complex *16 :: mu(2), nu(2)

      integer, allocatable :: ipiv(:)
      
      complex *16 :: zk, q1, q2
      
      real *8 pars(100), src(2), targ(2)
c     
      external fcurve,fcurve2, fcurve3
      external zhelmstokes_kern, fgreensdummy
c
      character trans

      trans = 'T'
      call prini(6,13)
      done=1
      pi=4*atan(done)
      ima=(0,1)

c     
c     
c     define points in exterior and interior of R.
c     
      xout=-6.1d0
      yout=-7.2d0

      xin = 0.1d0
      yin = -0.21d0
c     
      xyout(1)=xout
      xyout(2)=yout

      call prin2('xyout=*',xyout,2)

      xyin(1) = xin
      xyin(2) = yin

c     
c     exact solution is a helmholtz stokes velocity 
c     field induced by the vector density muexact
c     

      zk = (1.0d0,0.0d0)

      muexact(1) = -0.5d0+hkrand(0)
      muexact(2) = -0.5d0+hkrand(0)

c     
c     define then chunk up domain
c     
      eps10=1.0d-12
      nover1=4
      k=16
c     
      maxchunks=nmax
      allocate(chunks(2,k,maxchunks))
      allocate(ders(2,k,maxchunks))
      allocate(ders2(2,k,maxchunks))
c     
      ifclosed=1
      chsmall=1000
      ta=0
      tb=2*pi

      t = 0.3d0

      pert = 0.2
      do i=1,100
         pars(i) = pert*hkrand(0)
      enddo
c     
      call chunkfunc(eps10,ifclosed,chsmall,ta,tb,fcurve3,
     1     pars,nover1,k,nch,chunks,adjs,ders,ders2,hs)
      hsmax = 0
      do i=1,nch
         if(hsmax.le.hs(i)) hsmax=hs(i)
         do j=1,k
c     c              write(72,*) chunks(1,j,i), chunks(2,j,i)
         enddo
      enddo
      call prin2('hsmax=*',hsmax,1)

      nt = k*nch

      lwgeo = 6*k*nch + 3*nch + 1000
      call prinf('lwgeo = *',lwgeo,1)
      allocate(wgeo(lwgeo))
      call chunkpack(k,nch,chunks,adjs,ders,ders2,
     1     hs,wgeo,lused)
      if (lused .gt. lwgeo) then
         write(*,*) 'INSUFFICIENT MEMORY FOR WGEO: ABORT'
         stop
      endif

      
c
c     allocate variables defined on boundary
c

      allocate(dsdt(k,maxchunks))
      allocate(smu(2,k,maxchunks),dmu(2,k,maxchunks))
      allocate(rhs(2,k,maxchunks))
      allocate(rnorms(2,k,maxchunks))
      allocate(taus(2,k,maxchunks))

      allocate(whts(k,nch))
      call chunkwhts(k,nch,chunks,ders,hs,whts)
      
      do i=1,nch
         do j=1,k
            rnorm = sqrt(ders(1,j,i)**2 + ders(2,j,i)**2)
            rnorms(1,j,i) = ders(2,j,i)/rnorm
            rnorms(2,j,i) = -ders(1,j,i)/rnorm
            taus(1,j,i) = ders(1,j,i)/rnorm
            taus(2,j,i) = ders(2,j,i)/rnorm
c     
            dsdt(j,i) = rnorm*hs(i)
c     
            ifstress=1
            call zhelmstokeslet(zk,xyout,chunks(1,j,i),muexact,
     1           veltemp,ifstress,stresstemp)

c     right hand side is the value of the velocity field

            rhs(1,j,i) = veltemp(1)
            rhs(2,j,i) = veltemp(2)

         enddo
      enddo

c     allocate matrices

      nmat = 2*k*nch
      call prinf('nmat = *',nmat,1)
      allocate(sysmat(nmat,nmat),stat=mystat)
      if (mystat .ne. 0) write(*,*) 'fail 1'
      allocate(dmat(nmat,nmat),stat=mystat)
      if (mystat .ne. 0) write(*,*) 'fail 2'      
      allocate(onesmat(nmat,nmat),stat=mystat)
      if (mystat .ne. 0) write(*,*) 'fail 3'      

c     form matrices

c     double layer mat
      
      q1 = (0.0d0,0.0d0)
      q2 = (1.0d0,0.0d0)
      ndim = 2
      write(*,*) 'building double layer mat ...'
      time1 = second()
      call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_kern,
     1     q1,q2,fgreensdummy,zk,pars1,pars2,ntot,
     2     dmat)
      time2 = second()
      write(*,*) 'done in ', time2-time1, ' seconds.'      
c     "ones" mat ( n(x) \outer \int  n(y) \cdot )
      write(*,*) 'building normal ones mat ...'
      call normalonesmat(onesmat,whts,rnorms,k,nch)
      write(*,*) 'done.'
      write(*,*) 'forming sysmat ...'
      do j = 1,ntot
         do i = 1,ntot
            sysmat(i,j) = dmat(i,j) + onesmat(i,j)
            if (i .eq. j) sysmat(i,j) = sysmat(i,j)-0.5d0
         enddo
      enddo
      write(*,*) 'done.'

      write(*,*) 'computing LU factors of sysmat ...'
      time1 = second()
      allocate(ipiv(ntot))
      call zgetrf( ntot, ntot, sysmat, ntot, ipiv, info )
      time2 = second()
      write(*,*) 'done in ', time2-time1, ' seconds.'

      write(*,*) 'triangular solve for density ...'
      do i = 1,nch
         do j = 1,k
            dmu(1,j,i) = rhs(1,j,i)
            dmu(2,j,i) = rhs(2,j,i)
            smu(1,j,i) = 0.0d0
            smu(2,j,i) = 0.0d0
         enddo
      enddo
      time1 = second()
      nrhs = 1
      call zgetrs('N', ntot, nrhs, sysmat, ntot, ipiv, dmu, ntot, info )
      time2 = second()
      write(*,*) 'done in ', time2-time1, ' seconds.'

c     
c     Compute velocity due to double layer
c     

      velsum(1) = 0.0d0
      velsum(2) = 0.0d0
      velsums(1) = 0.0d0
      velsums(2) = 0.0d0
      velsumd(1) = 0.0d0
      velsumd(2) = 0.0d0

      dmaxwht = 0.0d0
      sumwht = 0.0d0
      
      do i = 1,nch
         do j = 1,k
            ifder = 0
            ifstress=0
            src(1) = chunks(1,j,i)
            src(2) = chunks(2,j,i)
            mu(1) = smu(1,j,i)
            mu(2) = smu(2,j,i)
            call zhelmstokeslet(zk,src,xyin,mu,
     1           veltemp,ifstress,stresstemp)
            velsums(1) = velsums(1)+veltemp(1)*whts(j,i)
            velsums(2) = velsums(2)+veltemp(2)*whts(j,i)
            
            mu(1) = dmu(1,j,i)
            mu(2) = dmu(2,j,i)
            nu(1) = rnorms(1,j,i)
            nu(2) = rnorms(2,j,i)
            trans='T'
            call zhelmstresslet(zk,src,xyin,mu,
     1           nu,veltemp,trans)
            velsumd(1) = velsumd(1)+veltemp(1)*whts(j,i)
            velsumd(2) = velsumd(2)+veltemp(2)*whts(j,i)
         enddo
      enddo
      
      velsum(1) = velsums(1) + velsumd(1)
      velsum(2) = velsums(2) + velsumd(2)

      call zhelmstokeslet(zk,xyout,xyin,muexact,
     1     velexact,ifstress,stresstemp)

c      write(*,*) velexact(1), velexact(2)
c      write(*,*) velsum(1), velsum(2)
      write(*,*) velexact(1)/velsum(1), velexact(2)/velsum(2)
c      write(*,*) velsums(1), velsums(2)
c      write(*,*) velsumd(1), velsumd(2)

      stop
      end
c     
c     
c
      subroutine normalonesmat(onesmat,whts,rnorms,k,nch)
c
c     This subroutine returns the matrix which
c     gives n(x) \int n(y) \cdot f(y) \, dy
c
c     Input:
c
c     whts - the smooth integration weights on the
c     boundary chunks
c     rnorms - the outward normal along the boundary
c     k - the number of points per chunk
c     nch - the number of chunks on the boundary
c     
      implicit none
c     global
      complex *16 onesmat(2*k*nch,2*k*nch)
      real *8 rnorms(2,k,nch), whts(k,nch)
      integer  k, nch
c     local
      integer i1,i2,i3,j1,j2,j3,imat,jmat
      real *8 rnormtemp, wht


      jmat = 1

      do j1 = 1,nch
      do j2 = 1,k
      do j3 = 1,2
         imat = 1
         rnormtemp = rnorms(j3,j2,j1)
         wht = whts(j2,j1)
         do i1 = 1,nch
         do i2 = 1,k
         do i3 = 1,2
            onesmat(imat,jmat) =
     1           wht*rnormtemp*rnorms(i3,i2,i1)
            imat = imat+1
         enddo
         enddo
         enddo
         jmat = jmat+1
      enddo
      enddo
      enddo

      return
      end
c     
c     
c     
      subroutine fcurve3(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
      implicit real *8 (a-h,o-z)
      integer n
      real *8 pars(100), t, pert, x0,y0, r0,r,rp,rpp
c     
      x0=0.0d0
      y0=0.0d0
      n = 6
      r0 = 5
c     
      r = r0
      rp = 0
      rpp = 0
      do i=1,n
         r = r + pars(i)*dsin(t*i)
         rp = rp + pars(i)*dcos(t*i)*i
         rpp = rpp - pars(i)*dsin(t*i)*i*i
      enddo
      x=x0+r*cos(t)
      y=y0+r*sin(t)

c     call prin2('t=*',t,1)


      dxdt=-r*sin(t) + rp*cos(t)
      dydt=r*cos(t) + rp*sin(t)

      dxdt2= -r*cos(t) + rpp*cos(t) - 2*rp*sin(t)
      dydt2= -r*sin(t) + rpp*sin(t) + 2*rp*cos(t)
c     
      return
      end
c     
c
