
      program testhelmbhders
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This program tests the difference Green's function
c     evaluation routines in helmbhrouts.f
c     against high-precision computations 
c     produced in Mathematica.
c
c     The test data files are formatted as
c     follows
c
c     test2pars.txt contains the parameters
c
c     ntheta, nrscale, nperr
c
c     The data loops over the zk values nrsclae*nperr
c     times
c
c     zk = e^(i*k*pi/ntheta), k = 0,ntheta-1 (this is different
c                       from the testhelmbhrouts.f file!!!)
c
c     each row of test2vals.txt is 
c
c     Re[zk] Im[zk] x1 x2 y1 y2 Re[g] Im[g] 
c
c     Where 
c
c     g(zk,x1,x2,y1,y2) = i/4*(H_0(zk*r)-2*i*log(r)/pi)
c
c     r = sqrt((x1-y1)^2 + (x2-y2)^2)
c
c     the nth order derivative data is stored in
c     test2valsdn.txt up to n = 5
c

      implicit real *8 (a-h,o-z)
      real *8 dtoobig, dtoosmall
      complex *16 eye
      data eye / (0.0d0,1.0d0) /
      complex *16, allocatable :: zks(:,:)
      real *8, allocatable :: zx(:,:,:), zy(:,:,:)
      real *8 errmax(0:5)
      complex *16 pot1, grad1(2), hess1(3), der31(4)
      complex *16 der41(5), der51(6)      
      complex *16, allocatable :: pot(:,:)
      complex *16, allocatable :: grad(:,:,:)
      complex *16, allocatable :: hess(:,:,:)
      complex *16, allocatable :: der3(:,:,:)
      complex *16, allocatable :: der4(:,:,:)
      complex *16, allocatable :: der5(:,:,:)            
      
c     parameters
      dtoobig = 1.0d300
      dtoosmall = 1.0d-300

      open(unit=20,file="../test/test-data/test2pars.txt")
      open(unit=21,file="../test/test-data/test2vals.txt")
      open(unit=31,file="../test/test-data/test2valsd1.txt")
      open(unit=32,file="../test/test-data/test2valsd2.txt")
      open(unit=33,file="../test/test-data/test2valsd3.txt")
      open(unit=34,file="../test/test-data/test2valsd4.txt")
      open(unit=35,file="../test/test-data/test2valsd5.txt")

      read(20,*) ntheta
      read(20,*) nrscale
      read(20,*) nperr

      write(*,*) "PARAMETERS (ntheta, nrscale, nperr)"
      write(*,*) ntheta, nrscale, nperr

      allocate(zks(nrscale*nperr,ntheta))
      allocate(zx(2,nrscale*nperr,ntheta))
      allocate(zy(2,nrscale*nperr,ntheta))
      allocate(pot(nrscale*nperr,ntheta))
      allocate(grad(2,nrscale*nperr,ntheta))
      allocate(hess(3,nrscale*nperr,ntheta))
      allocate(der3(4,nrscale*nperr,ntheta))
      allocate(der4(5,nrscale*nperr,ntheta))
      allocate(der5(6,nrscale*nperr,ntheta))

      ifpot = 1
      ifgrad = 1
      ifhess = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 1

      write(*,*) "READING IN DATA ..."
      do i = 1,ntheta
      do ii = 1,nrscale*nperr
         read(21,*,iostat=ios) dkr, dki, x1, x2, y1, y2,
     1        dgr, dgi 
         if (ios .ne. 0) then
            write(*,*) "bad read"
            stop
         endif
         read(31,*,iostat=ios) dgxr, dgxi, dgyr, dgyi
         if (ios .ne. 0) then
            write(*,*) "bad read d1"
            stop
         endif
         read(32,*,iostat=ios) dgxxr, dgxxi, dgxyr, dgxyi,
     1        dgyyr, dgyyi
         if (ios .ne. 0) then
            write(*,*) "bad read d2"
            stop
         endif
         read(33,*,iostat=ios) dgxxxr, dgxxxi, dgxxyr, dgxxyi,
     1        dgxyyr, dgxyyi, dgyyyr, dgyyyi
         if (ios .ne. 0) then
            write(*,*) "bad read d3"
            stop
         endif
         read(34,*,iostat=ios) dgxxxxr, dgxxxxi, dgxxxyr, dgxxxyi,
     1        dgxxyyr, dgxxyyi, dgxyyyr, dgxyyyi, dgyyyyr, dgyyyyi
         if (ios .ne. 0) then
            write(*,*) "bad read d4"
            stop
         endif
         read(35,*,iostat=ios) dgxxxxxr, dgxxxxxi, dgxxxxyr,dgxxxxyi,
     1        dgxxxyyr, dgxxxyyi, dgxxyyyr, dgxxyyyi, dgxyyyyr,
     2        dgxyyyyi, dgyyyyyr, dgyyyyyi
         if (ios .ne. 0) then
            write(*,*) "bad read d5"
            stop
         endif
         zks(ii,i) = dkr+eye*dki
         zx(1,ii,i) = x1
         zx(2,ii,i) = x2
         zy(1,ii,i) = y1
         zy(2,ii,i) = y2
         pot(ii,i) = (dgr + eye*dgi)
         grad(1,ii,i) = (dgxr + eye*dgxi)
         grad(2,ii,i) = (dgyr + eye*dgyi)
         hess(1,ii,i) = (dgxxr + eye*dgxxi)
         hess(2,ii,i) = (dgxyr + eye*dgxyi)
         hess(3,ii,i) = (dgyyr + eye*dgyyi)
         der3(1,ii,i) = (dgxxxr + eye*dgxxxi)
         der3(2,ii,i) = (dgxxyr + eye*dgxxyi)
         der3(3,ii,i) = (dgxyyr + eye*dgxyyi)
         der3(4,ii,i) = (dgyyyr + eye*dgyyyi)
         der4(1,ii,i) = (dgxxxxr + eye*dgxxxxi)
         der4(2,ii,i) = (dgxxxyr + eye*dgxxxyi)
         der4(3,ii,i) = (dgxxyyr + eye*dgxxyyi)
         der4(4,ii,i) = (dgxyyyr + eye*dgxyyyi)
         der4(5,ii,i) = (dgyyyyr + eye*dgyyyyi)
         der5(1,ii,i) = (dgxxxxxr + eye*dgxxxxxi)
         der5(2,ii,i) = (dgxxxxyr + eye*dgxxxxyi)
         der5(3,ii,i) = (dgxxxyyr + eye*dgxxxyyi)
         der5(4,ii,i) = (dgxxyyyr + eye*dgxxyyyi)
         der5(5,ii,i) = (dgxyyyyr + eye*dgxyyyyi)
         der5(6,ii,i) = (dgyyyyyr + eye*dgyyyyyi)
      enddo
      enddo

      do i = 0,5
         errmax(i) = 0.0d0
      enddo
               
      write(*,*) "RUNNING TESTS"

      write(*,*) "only for theta=0,pi/2"
      write(*,*) "other angles can have branch cut problems"

      do i = 1,ntheta/2+1,ntheta/2
      do ii = 1,nrscale*nperr
         call helmbhgreen_all(zks(ii,i),zx(1,ii,i),zy(1,ii,i),
     1        ifpot,pot1,ifgrad,grad1,ifhess,hess1,ifder3,der31,
     2        ifder4,der41,ifder5,der51)
         
         
         err = cdabs(pot1-pot(ii,i))/cdabs(pot(ii,i))
         if (err .gt. errmax(0)) errmax(0) = err
         err = cdabs(grad1(1)-grad(1,ii,i))/(cdabs(grad(1,ii,i)))
         if (err .gt. errmax(1)) errmax(1) = err
         err = cdabs(grad1(2)-grad(2,ii,i))/(cdabs(grad(2,ii,i)))
         if (err .gt. errmax(1)) errmax(1) = err
         err = cdabs(hess1(1)-hess(1,ii,i))/(cdabs(hess(1,ii,i)))
         if (err .gt. errmax(2)) errmax(2) = err
         err = cdabs(hess1(2)-hess(2,ii,i))/(cdabs(hess(2,ii,i)))
         if (err .gt. errmax(2)) errmax(2) = err
         err = cdabs(hess1(3)-hess(3,ii,i))/(cdabs(hess(3,ii,i)))
         if (err .gt. errmax(2)) errmax(2) = err
         err = cdabs(der31(1)-der3(1,ii,i))/(cdabs(der3(1,ii,i)))
         if (err .gt. errmax(3)) errmax(3) = err
         err = cdabs(der31(2)-der3(2,ii,i))/(cdabs(der3(2,ii,i)))
         if (err .gt. errmax(3)) errmax(3) = err
         err = cdabs(der31(3)-der3(3,ii,i))/(cdabs(der3(3,ii,i)))
         if (err .gt. errmax(3)) errmax(3) = err
         err = cdabs(der31(4)-der3(4,ii,i))/(cdabs(der3(4,ii,i)))
         if (err .gt. errmax(3)) errmax(3) = err
         err = cdabs(der41(1)-der4(1,ii,i))/(cdabs(der4(1,ii,i)))
         if (err .gt. errmax(4)) errmax(4) = err
         err = cdabs(der41(2)-der4(2,ii,i))/(cdabs(der4(2,ii,i)))
         if (err .gt. errmax(4)) errmax(4) = err
         err = cdabs(der41(3)-der4(3,ii,i))/(cdabs(der4(3,ii,i)))
         if (err .gt. errmax(4)) errmax(4) = err
         err = cdabs(der41(4)-der4(4,ii,i))/(cdabs(der4(4,ii,i)))
         if (err .gt. errmax(4)) errmax(4) = err
         err = cdabs(der41(5)-der4(5,ii,i))/(cdabs(der4(5,ii,i)))
         if (err .gt. errmax(4)) errmax(4) = err
         err = cdabs(der51(1)-der5(1,ii,i))/(cdabs(der5(1,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
         err = cdabs(der51(2)-der5(2,ii,i))/(cdabs(der5(2,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
         err = cdabs(der51(3)-der5(3,ii,i))/(cdabs(der5(3,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
         err = cdabs(der51(4)-der5(4,ii,i))/(cdabs(der5(4,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
         err = cdabs(der51(5)-der5(5,ii,i))/(cdabs(der5(5,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
         err = cdabs(der51(6)-der5(6,ii,i))/(cdabs(der5(6,ii,i)))
         if (err .gt. errmax(5)) errmax(5) = err
      enddo
      enddo

      write(*,*) "WORST CASE RELATIVE ERROR BY DERIVATIVE ORDER"
      
      do i = 0,5
         write(*,*) i, errmax(i)
      enddo

      stop
      end
