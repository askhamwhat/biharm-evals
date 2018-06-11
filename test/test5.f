      implicit real *8 (a-h,o-z)

      real *8 roots(10000),pars(10)
      external fval1

      call prini(6,13)
      call prinf('Enter n*',n,0)
      read *, n

      a = 1.1d0
      b = 2.3d0

c      n = 20
      h = 1.0d-4

      nmax = 100

      call chebroots(fval1,pars,a,b,n,h,rerr,nmax,nroots,roots)
      call prin2('roots=*',roots,nroots)

      do i=1,nroots

      call fval1(roots(i),pars,fx)
      call prin2('fx=*',fx,1)

      enddo

      stop
      end
c-------------------------------------------------------------
      
      subroutine chebroots(fker,pars,a,b,n,h,rerr,nmax,nroots,roots)
      implicit real *8 (a-h,o-z)

c
cc       this subrotine computes the roots of the function
c        fker on the interval [a,b] using boyd's method.
c
c        input 
c        fker - function whose roots we wish to find
c        the calling sequence for fker is 
c         subroutine fker(x,pars,fx)
c        pars - parameters for the function fker
c        a - left end of interval
c        b - right end of interval
c        n - number of terms in the chebyshev expansion of the function
c        h - spacing between points for finding the zeros
c
c       
c        output
c        rerr - estimate for error in function interpolation
c        nroots - number of roots computed
c        roots - location of the roots
c        
c
      real *8 tvals(n),fvals(n),fcoefs(n),roots(*)
      real *8 xnodes(n),whts(n),u(n,n),v(n,n),pars(*)
      real *8 pars2(n+10)

      external fker,fkerchebexp

      done = 1
      pi = atan(done)*4

      call chebexps(2,n,xnodes,u,v,whts)

      do i=1,n

      tvals(i) = a + (b-a)/2*(xnodes(i)+1)

      call fker(tvals(i),pars,fvals(i))

      enddo

      call pyplot(11,tvals,fvals,n,3,'a*')
      call chematvec(u,fvals,fcoefs,n)

      call prin2('fcoefs=*',fcoefs,n)

      rerr = abs(fcoefs(n)/fcoefs(1))

      pars2(1) = n+0.1
      pars2(2) = a
      pars2(3) = b
      do i=1,n

      pars2(i+3) = fcoefs(i)

      enddo


c
cc      find roots with spacing h
c

      call findroots(fkerchebexp,pars2,a,b,h,nmax,roots,nroots)


      
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

      subroutine fval1(x,pars,fx)
      implicit real *8 (a-h,o-z)
      real *8 pars(*)

      fx = sin(12.1d0*x) 


      return
      end
c------------------------------------------------

      subroutine fkerchebexp(x,pars,fx)
      implicit real *8 (a-h,o-z)
      real *8 pars(*)
c
cc      this subroutine evaluates a chebyshev expansion
c       at a point x on the interval [a,b]
c 
c     input
c     x - point in [a,b] where the chebyshev expansion is to be
c         evaluated
c     pars(1) - order of chebyshev expansion
c     pars(2) - a
c     pars(3) - b
c     pars(4:n+3) - chebyshev expansion coefficients
      

      n = pars(1)
      a = pars(2)
      b = pars(3)

      xeval = -1 + (x-a)/(b-a)*2
      fx = 0
      call chebexev(xeval,fx,pars(4),n-1)


      return
      end
c------------------------------------------------

      
      
      
      
