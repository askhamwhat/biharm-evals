      implicit real *8 (a-h,o-z)
      real *8 pars(1000),roots(1000)

      external fdet1,fdet2

      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n

      r1 = 1.0d0
      r2 = 1.7d0

      pars(1) = r1
      pars(2) = r2


c
cc       find neumann eigenfunctions using bisection
c
c

      x = 5.0d0
      call fdet1(x,pars,f)

      a = 1.0d0
      b = 3.0d1

      h = 1.0d-2

      nmax = 10
      call findroots(fdet1,pars,a,b,h,nmax,roots,nroots)

      call prinf('nroots=*',nroots,1)
      call prin2('roots=*',roots,nroots)

      f = 0
      call fdet1(roots(1),pars,f)
      call prin2('f=*',f,1)


      a = 1.0d0
      b = 3.0d1
      h = 1.0d-2

      nmax = 100
      call findroots(fdet2,pars,a,b,h,nmax,roots,nroots)

      call prinf('nroots=*',nroots,1)
      call prin2_long('roots=*',roots,nroots)

      call fdet2(roots,pars,f)
      call prin2('f=*',f,1)

      return
      end
c----------------------------------      

      subroutine fdet1(x,pars,f)
      implicit real *8 (a-h,o-z)
      integer ifexpon
      real *8 pars(*),xtmp(6)
      complex *16 h1,h0,z1,z2

      ifexpon = 1
      r1 = pars(1)
      r2 = pars(2)

      z1 = x*r1

      h0 = 0
      h1 = 0
      ifexpon = 1 
      call hank103(z1,h0,h1,ifexpon)
      
      a11 = real(h1)
      a12 = imag(h1)

      z2 = cmplx(x*r2,0)

      h1 = 0
      call hank103(z2,h0,h1,ifexpon)
      a21 = -real(h1)
      a22 = -imag(h1)

      f = a11*a22 - a12*a21

      xtmp(1) = r1
      xtmp(2) = r2
      xtmp(3) = a11
      xtmp(4) = a12
      xtmp(5) = a21
      xtmp(6) = a22
cc      call prin2('xtmp=*',xtmp,6)
      
      

      return
      end
c------------------------------------------
      subroutine fdet2(x,pars,f)
      implicit real *8 (a-h,o-z)
      integer ifexpon,ipvt(1000)
      real *8 pars(*),xtmp(6)
      real *8 amat(4,4),work(10000),det(2)
      complex *16 h1r1,h0r1,h0r2,h1r2,z1,z2

      ifexpon = 1

      r1 = pars(1)
      r2 = pars(2)

      z1 = x*r1

      h0 = 0
      h1 = 0
      ifexpon = 1 
      call hank103(z1,h0r1,h1r1,ifexpon)

      z2 = cmplx(x*r2,0)

      h1 = 0
      call hank103(z2,h0r2,h1r2,ifexpon)

      amat(1,1) = 1
      amat(2,1) = 1
      amat(3,1) = 0
      amat(4,1) = 0

      amat(1,2) = log(r1)
      amat(2,2) = log(r2)
      amat(3,2) = 1.0d0/r1
      amat(4,2) = -1.0d0/r2

      amat(1,3) = real(h0r1)
      amat(2,3) = real(h0r2)
      amat(3,3) = -x*real(h1r1)
      amat(4,3) = x*real(h1r2)
      
      amat(1,4) = imag(h0r1)
      amat(2,4) = imag(h0r2)
      amat(3,4) = -x*imag(h1r1)
      amat(4,4) = x*imag(h1r2)

cc      call prin2('amat=*',amat,16)

      job = 10
      det = 0
      call dgeco(amat,4,4,ipvt,rcond,work)
cc      call prin2('rcond=*',rcond,1)
      call dgedi(amat,4,4,ipvt,det,work,job)

cc      call prin2_long('x=*',x,1)
cc      call prin2_long('det=*',det,2)
      f = det(1)*10.0d0**(det(2))

cc      call prin2('f=*',f,1)

      return
      end
c-------------------------------------------------------------
     
      subroutine findroots(fker,pars,a,b,h,nmax,roots,nroots)
      implicit real *8 (a-h,o-z)
c
cc       this subroutine computes the roots of a chebyshev expansion
c        on the interval [a,b] by testing at equispaced points with
c        spacing h

      real *8 roots(*),pars(*)
      external fker


        x0 = a

        nh = (b-a)/h

        call fker(x0,pars,val)
        ione = 1
        if(val.ge.0) ipar = 1
        if(val.lt.0) ipar = -1

        nroots = 0

        do 2000 i=1,nh
        x0 = x0 + h
        call fker(x0,pars,val)

        if(val.ge.0) ipar2 = 1
        if(val.lt.0) ipar2 = -1
cc      
c           .   .   . bissect
cc
        if (ipar .ne. ipar2) then
        
        xa = x0-h
        xb = x0
        do 1500 jj=1,64
        xm = (xa+xb)/2
        call fker(xm,pars,val)
        if(val.ge.0) iparm = 1
        if(val.lt.0) iparm = -1

        if (iparm .eq. ipar) xa = xm
        if (iparm .eq. ipar2) xb = xm
 1500 continue
        
        nroots = nroots+1

        if(nroots.gt.nmax) then

        call prin2('not enough memory*',i,0)
        stop

        endif
        roots(nroots) = xm

        ipar = ipar2

        endif

 2000 continue

      return
      end


c----------------------------------------------------------
