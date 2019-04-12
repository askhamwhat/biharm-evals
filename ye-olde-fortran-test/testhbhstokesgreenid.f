
      program testhbhstokesgreenid
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This program tests the Green's identity for the 
c     Helmholtz Stokes single and double layer kernels
c
      
      implicit real *8 (a-h,o-z)
      parameter (nmax = 100000)
      integer adjs(2,nmax)
      real *8 hs(nmax),xyin(10),xyout(10)
c     
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),
     1     ders2(:,:,:),dsdt(:,:),rnorms(:,:,:),
     2     whts(:,:),taus(:,:,:)
      complex *16, allocatable :: smu(:,:,:), dmu(:,:,:)
      complex *16, allocatable :: smutau(:,:), dmutau(:,:)
      complex *16, allocatable :: smunu(:,:), dmunu(:,:)

      complex *16 :: veltemp(2), stresstemp(2,2), muexact(2)
      complex *16 :: veldertemp(2,2)
      complex *16 :: velsum(2), velexact(2), velsums(2), velsumd(2)
      complex *16 :: mu(2), nu(2)

      complex *16 :: zk
      
      real *8 pars(100), src(2), targ(2)
c     
      external fcurve,fcurve2, fcurve3
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
      nover1=6
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

c
c     set up densities for Green's identity
c

      allocate(dsdt(k,maxchunks))
      allocate(smu(2,k,maxchunks),dmu(2,k,maxchunks))
      allocate(smutau(k,maxchunks),dmutau(k,maxchunks))
      allocate(smunu(k,maxchunks),dmunu(k,maxchunks))
      allocate(rnorms(2,k,maxchunks))
      allocate(taus(2,k,maxchunks))

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

c     density for single layer (normal component of stress)
            
            smu(1,j,i) = stresstemp(1,1)*rnorms(1,j,i) 
     1           + stresstemp(1,2)*rnorms(2,j,i)
            smu(2,j,i) = stresstemp(2,1)*rnorms(1,j,i) 
     1           + stresstemp(2,2)*rnorms(2,j,i)

            smutau(j,i) = smu(1,j,i)*taus(1,j,i)
     1           +smu(2,j,i)*taus(2,j,i)
            smunu(j,i) = smu(1,j,i)*rnorms(1,j,i)
     1           +smu(2,j,i)*rnorms(2,j,i)

c     density for double layer (velocity)

            dmu(1,j,i) = veltemp(1)
            dmu(2,j,i) = veltemp(2)

            dmutau(j,i) = dmu(1,j,i)*taus(1,j,i)
     1           +dmu(2,j,i)*taus(2,j,i)
            dmunu(j,i) = dmu(1,j,i)*rnorms(1,j,i)
     1           +dmu(2,j,i)*rnorms(2,j,i)


         enddo
      enddo


c
c     Compute velocity from Green's I.D.
c     

      allocate(whts(k,nch))
      call chunkwhts(k,nch,chunks,ders,hs,whts)

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
      
      velsum(1) = -(velsums(1) + velsumd(1))
      velsum(2) = -(velsums(2) + velsumd(2))

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
