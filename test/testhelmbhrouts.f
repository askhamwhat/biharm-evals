
      program testhelmbhrouts
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This program tests the difference function
c     evaluation routines in helmbhrouts.f
c     against high-precision computations 
c     produced in Mathematica.
c
c     The test data files are formatted as
c     follows
c
c     test1pars.txt contains the parameters
c
c     ntheta, nrscale, nterms
c
c     The data loops over the zk, r, j values
c     zk is outermost, then r, then j is innermost
c
c     zk = e^(i*k*2*pi/ntheta), k = 1,ntheta
c     r = 10^(offset+k), k = 1,nrscale (offset = -15)
c     j = 0,nterms
c
c     each row of test1vals.txt is 
c
c     j Re[zk] Im[zk] r rscale Re[diffout_j[zk,r]] Im[diffout_j[zk,r] 
c     Re[diffin_j[zk,r]] Im[diffin_j[zk,r] 
c
c     Where 
c
c     diffout_j[zk,r] = H_j[zk*r]+i*2^(j-1)*(j-1)!*2/(pi*(zk*r)^j)
c
c     diffin_j[zk,r] = J_j[zk*r]-(zk*r)^j/(2^j*j!)
c

      implicit real *8 (a-h,o-z)
      real *8 dtoobig, dtoosmall
      complex *16 eye
      data eye / (0.0d0,1.0d0) /
      complex *16, allocatable :: zks(:,:)
      complex *16, allocatable :: diffout(:,:,:), diffin(:,:,:)
      real *8, allocatable :: rs(:,:), rscales(:,:)
      integer, allocatable :: iftest(:,:,:)
      complex *16, allocatable :: hvec(:)
      complex *16, allocatable :: ders(:)
      complex *16, allocatable :: diffs(:)
      real *8, allocatable :: errmax(:), rsmax(:)
      complex *16, allocatable :: zksmax(:)
      real *8, allocatable :: errmaxin(:), rsmaxin(:)
      complex *16, allocatable :: zksmaxin(:)
      
c     parameters
      dtoobig = 1.0d300
      dtoosmall = 1.0d-300

      open(unit=20,file="../test/test-data/test1pars.txt")
      open(unit=21,file="../test/test-data/test1vals.txt")

      read(20,*) ntheta
      read(20,*) nrscale
      read(20,*) nterms

      write(*,*) "PARAMETERS (ntheta, nrscale, nterms)"
      write(*,*) ntheta, nrscale, nterms

      allocate(zks(nrscale,ntheta),rs(nrscale,ntheta),
     1     rscales(nrscale,ntheta))
      allocate(diffout(0:nterms,nrscale,ntheta),
     1     diffin(0:nterms,nrscale,ntheta))
      allocate(hvec(0:nterms),ders(0:nterms),diffs(0:nterms))
      allocate(zksmax(0:nterms),errmax(0:nterms),rsmax(0:nterms))
      allocate(zksmaxin(0:nterms),errmaxin(0:nterms),
     1     rsmaxin(0:nterms))      

      write(*,*) "READING IN DATA ..."
      do i = 1,ntheta
         do j = 1,nrscale
            do k = 0,nterms
               read(21,*,iostat=ios) ktemp, dkr, dki, dr, drsc,
     1              dor, doi, dir, dii
               if (ios .ne. 0) then
                  write(*,*) "bad read"
                  stop
               endif
               if (ktemp .ne. k) write(*,*) "something went wrong"
               zks(j,i) = dkr+eye*dki
               rs(j,i) = dr
               rscales(j,i) = drsc
               diffout(k,j,i) = dor+eye*doi 
               diffin(k,j,i) = dir+eye*dii
            enddo
         enddo
      enddo

      do i = 0,nterms
         errmax(i) = 0.0d0
         zksmax(i) = 0.0d0
         rsmax(i) = 0.0d0
         errmaxin(i) = 0.0d0
         zksmaxin(i) = 0.0d0
         rsmaxin(i) = 0.0d0
      enddo
               
      write(*,*) "RUNNING TEST (ON FIRST QUADRANT)..."

      do i = 1,ntheta/4
         do j = 1,nrscale-2
            
            call diffsloghank(rs(j,i),zks(j,i),rscales(j,i),
     1           diffs,ifders,ders,hvec,nterms)
            do k = 0,nterms
               err = abs(diffs(k)-diffout(k,j,i))/abs(diffout(k,j,i))
               if (err .gt. errmax(k)) then
                  errmax(k) = err
                  zksmax(k) = zks(j,i)
                  rsmax(k) = rs(j,i)
               endif
            enddo

            call diffszkjfun(rs(j,i),zks(j,i),rscales(j,i),
     1           diffs,ifders,ders,hvec,nterms)
            do k = 0,nterms
               err = abs(diffs(k)-diffin(k,j,i))/abs(diffin(k,j,i))
               if (err .gt. errmaxin(k)) then
                  errmaxin(k) = err
                  zksmaxin(k) = zks(j,i)
                  rsmaxin(k) = rs(j,i)
               endif
            enddo

            
            
         enddo
      enddo

      write(*,*) "WORST CASE ERROR BY MODE"

      do i = 0,nterms
         write(*,*) i, errmax(i), errmaxin(i)
      enddo

      stop
      end
