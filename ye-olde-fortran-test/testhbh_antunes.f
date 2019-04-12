c-----------------------2--------------------------
c
c This file tests the routines for the eigenvalues
c of a clamped plate on a domain from Antunes 2011, 
c for which the first eigenvalue should be 9.466464562432
c (k = 3.07676202564)
c this domain is interesting b/c it is a relatively
c simple shape for which the first eigenfunction has 
c a sign change (it is a convex domain with analytic
c boundary)
c
c In particular, this routine takes in points A and B
c and an integer N (a power of 2 is a good choice)
c and computes the determinant of the appropriately
c scaled system matrix for the value k at each of the 
c N+1 second kind chebyshev points on [a,b]
c
c


      program main
      implicit real *8 (a-h,o-z)
      real *8 w(1000000)
      real *8 a, b, pi
      complex *16, allocatable :: zks(:), zdets(:)
      complex *16 q1, q2
      integer n, nchs(100)
      real *8 cms(2,100)

      integer, allocatable :: adjs(:,:)
      real *8, allocatable :: hs(:)
c     
      real *8, allocatable :: chunks(:,:,:),ders(:,:,:),
     1     ders2(:,:,:), wgeo(:), xs(:), ys(:)

      real *8 pars(100)
c     
      external fcurve_antunes

      character *100 filename, format_string

      call prini(6,13)

      pi = 4.0d0*datan(1.0d0)

c
c     define then chunk up domain
c     
      eps10=1.0d-14
      nover1=1
      k=16
c     
      maxchunks=1000
      allocate(chunks(2,k,maxchunks),hs(maxchunks),adjs(2,maxchunks))
      allocate(ders(2,k,maxchunks))
      allocate(ders2(2,k,maxchunks))
c     
      ifclosed=1
      chsmall=1000
      ta=0
      tb=2*pi

      call chunkfunc(eps10,ifclosed,chsmall,ta,tb,fcurve_antunes,
     1     pars,nover1,k,nch,chunks,adjs,ders,ders2,hs)

      if (nch .gt. maxchunks) then
         write(*,*) 'MAX CHUNKS EXCEEDED. ABORT'
         stop
      endif

      call prinf('k*',k,1)
      call prinf('nch*',nch,1)
      call prinf('npts*',nch*k,1)

      allocate(xs(nch*k),ys(nch*k))
      do i = 1,nch
         do j = 1,k
            xs((i-1)*k+j) = chunks(1,j,i)
            ys((i-1)*k+j) = chunks(2,j,i)
         enddo
      enddo

      ifile = 37
      call pyplot(ifile,xs,ys,nch*k,2,'a*')

      lwgeo = 1000000
      allocate(wgeo(lwgeo))

      call chunkpack(k,nch,chunks,adjs,ders,ders2,hs,wgeo,lused)
      if (lused .gt. lwgeo) then
         write(*,*) 'LWGEO EXCEEDED. ABORT'
         stop
      endif

c     this is a simply connected domain (cms(1:2,1) is a point outside)

      ncomp = 1
      nchs(1) = nch
      cms(1,1) = 3.0d0
      cms(2,1) = 3.0d0
      


      call prin2('Enter a (left end-point)*',n,0)
      read *, a
      call prin2('Enter b* (right end-point)*',n,0)
      read *, b
      call prin2('Enter n (integer bigger than 0)*',n,0)
      read *, n

      if (b .lt. a .or. a .lt. 1d-15 .or. n .lt. 1) then
         write(*,*) 'bad values for a,b,or n'
         stop
      endif

      allocate(zks(0:n),zdets(0:n))

      do i = 0,n
         zks(i) = 0.5*(a+b) + 0.5*(b-a)*cos((n-i)*pi/n)
      enddo

c     set q2 so that the operator is of the form I + K

      q1 = (0.0d0,0.0d0)
      q2 = (-2.0d0,0.0d0)

c     flag for square root scaling
      ifss = 1

      do i = 0,n
         write(*,*) '******************************'
         write(*,*) '******************************'
         write(*,*) 'i = ', i
         write(*,*) 'zks(i) ', zks(i)
         
         call prin2('evaluating determinant ...*',zdets,0)
         call zhbhdir_det(zks(i),wgeo,ncomp,nchs,cms,q1,q2,
     1     ifss,zdets(i))

         call prin2('found determinant*',zdets(i),2)
      enddo
      

      iw = 16

      format_string = "(A7,I0.6)"
      write(filename,format_string) "antunes", n
      open(unit=iw,FILE=filename)


      do i=0,n

      write(iw,*) real(zks(i)), imag(zks(i)), real(zdets(i)), 
     1        imag(zdets(i))

      enddo

      stop
      end

c------------------------------------------------------     
c
      subroutine zhbhdir_det(zk,wgeo,ncomp,nchs,cms,q1,q2,
     1     ifss,zdetval)

      implicit none
c     inputs 
      real *8 :: wgeo(*), cms(2,*)
      complex *16 zk, q1, q2, zdetval
      integer :: ncomp, nchs(*), ifss

c     local
      complex *16, allocatable :: sysmat(:,:)
      real *8, allocatable :: whts(:,:)
      integer ntot, k, nch, ichunks,iadjs,iders,iders2,ihs
      integer :: ier, ich, jch, inode, jnode, info, ipt, jpt
      integer, allocatable :: ipiv(:)
      real *8 :: time1, time2, rr
      

      call chunkunpack1(wgeo,k,nch,ichunks,iadjs,iders,iders2,ihs)

      allocate(whts(k,nch))

      call chunkwhts(k,nch,wgeo(ichunks),wgeo(iders),wgeo(ihs),whts)

      ntot = 2*k*nch + ncomp
      allocate(sysmat(ntot,ntot))

      call prinf('k=*',k,1)
      call prinf('ntot=*',ntot,1)

      call prinf('building system matrix ...*',ier,0)
      time1 = second()
      call zhbhstokesmatbuild(zk,wgeo,ncomp,nchs,cms,
     1     q1,q2,ntot,sysmat,ier)
      time2 = second()
      call prin2('done*',time2-time1,1)

c
c     square root scale the matrix (if requested)
c
      if (ifss .eq. 1) then
         
         call prinf('square root scaling matrix ...*',ifss,0)
         time1 = second()

         do ich = 1,nch
            do inode =1,k

               ipt = (ich-1)*k + inode

               do jch = 1,nch
                  do jnode = 1,k

                     jpt = (jch-1)*k + jnode

                     rr = sqrt(whts(inode,ich)/whts(jnode,jch))
                     sysmat(2*ipt-1,2*jpt-1) = 
     1                    sysmat(2*ipt-1,2*jpt-1)*rr
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
         time2 = second()
         call prin2('done*',time2-time1,1)

      endif

      call prinf('calling zgetrf ...*',ier,0)
      allocate(ipiv(ntot))
      time1 = second()
      call zgetrf( ntot, ntot, sysmat, ntot, ipiv, info )
      time2 = second()
      call prin2('done*',time2-time1,1)
      call prinf('info*',info,1)
      call prinf('extracting determinant ...*',info,0)
      call zgetr_det(ntot,sysmat,ipiv,info,zdetval)
      

      return
      end
c----------------------------------------------------------------------


      subroutine zgetr_det(n,a,ipiv,info,zdet)
      implicit none
      integer :: n, info
      complex *16 :: a(n,n), zdet
      integer :: ipiv(*)

      integer i

      zdet = 0d0
      if (info.ne.0) then
         return
      endif
      zdet = 1d0
      do i=1,n
         if (ipiv(i).ne.i) then
            zdet = -zdet * a(i,i)
         else
            zdet = zdet * a(i,i)
         endif
      enddo
      end


      subroutine fcurve_antunes(t,pars,x,y,dxdt,dydt,dxdt2,dydt2)
      implicit real *8 (a-h,o-z)
      integer n
      real *8 pars(100), t, pert, x0,y0, r0,r,rp,rpp
c     
      x = 2.0d0*cos(t)
      dxdt = -2.0d0*sin(t)
      dxdt2 = -2.0d0*cos(t)

      y = sin(t) + 0.25d0*sin(2.0d0*t)
      dydt = cos(t) + 0.5d0*cos(2.0d0*t)
      dydt2 = -sin(t) - sin(2.0d0*t)
c     
      return
      end
