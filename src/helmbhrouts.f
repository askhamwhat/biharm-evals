cc Copyright (C) 2018: Travis Askham
cc Contact: askhamwhat@gmail.com
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the BSD-3-Clause License
      
      subroutine helmbhgreen_all(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green function 
c
c     G(x,y) = (i*H_0^(1)(zk*r)/4 + log(r)/(2*pi))/zk^2 ,
c
c     where H_0^(1) is the Hankel function of the first kind and
c     r = sqrt( (x_1-y_1)^2 + (x_2-y_2)^2),
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     zk - complex *16 helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c     ifder3 - if equal to 1, compute 3rd order derivatives
c     ifder4 - if equal to 1, compute 4th order derivatives
c     ifder5 - if equal to 1, compute 5th order derivatives
c
c     OUTPUTS
c
c     note: for derivatives, the ordering is such that,
c     for a given order, the first entry has all of the 
c     derivatives on the x1 variables, the second entry
c     has one derivative on the x2 variable and the rest 
c     on x1, and so on. e.g. 
c
c     hess(1) = d_x1 d_x1 G(x,y)
c     hess(2) = d_x1 d_x2 G(x,y)
c     hess(3) = d_x2 d_x2 G(x,y)
c
c     pot - complex *16, potential, if requested
c     grad(2) - complex *16, gradient, if requested
c     hess(3) - complex *16, hessian, if requested
c     der3(4) - complex *16, 3rd order derivatives, if requested.
c     der4(5) - complex *16, 4th order derivatives, if requested.
c     der5(6) - complex *16, 5th order derivatives, if requested.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk
      real *8 zx(2), zy(2)
      complex *16 pot, grad(2), hess(3)
      complex *16 der3(4), der4(5), der5(6)
      integer ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
c     local variables
      real *8 r, r2, r3, r4, r5, r6
      real *8 dx, dy, dx2, dx3, dx4, dx5, dy2, dy3, dy4, dy5
      integer ifders
      real *8 rscale
      integer nmax, ifder, j
      complex *16 z, eye, eyeo4, hvec(0:10), diffs(0:10), ders(0:10)
      complex *16 zkh, ztemp, zkr
      data eye /(0.0d0,1.0d0)/
      data eyeo4 /(0.0d0,0.25d0)/

      dx = zx(1) - zy(1)
      dy = zx(2) - zy(2)

      dx2 = dx*dx
      dx3 = dx2*dx
      dx4 = dx3*dx
      dx5 = dx4*dx
      dy2 = dy*dy
      dy3 = dy2*dy
      dy4 = dy3*dy
      dy5 = dy4*dy

      r2 = dx2 + dy2
      r = dsqrt(r2)
      r3 = r2*r
      r4 = r3*r
      r5 = r4*r

      nmax = 5

c     get values of difference functions up to nmax-th order
      
      rscale = 1.0d0
      ifders = 0
      call diffsloghank(r,zk,rscale,diffs,ifders,ders,hvec,nmax)

c     rescale

      zkr = zk*r

c     compute r derivatives using only like order terms to minimize
c     cancellation

      ztemp = eyeo4/(zk**2)
      ders(0) = diffs(0)*ztemp
      ztemp = -zk*ztemp
      ders(1) = diffs(1)*ztemp
      ztemp = -zk*ztemp
      ders(2) = ztemp*(diffs(2)-diffs(1)/zkr)
      ztemp = -zk*ztemp
      ders(3) = ztemp*(diffs(3)-3.0d0*diffs(2)/zkr)
      ztemp = -zk*ztemp
      ders(4) = ztemp*(diffs(4)-(5.25d0*diffs(3)
     1     -0.75d0*hvec(1))/zkr)
      ztemp = -zk*ztemp
      ders(5) = ztemp*(diffs(5)-(7.5d0*diffs(4)
     1     -2.5d0*hvec(2))/zkr)
      
c     evaluate potential and derivatives

      if (ifpot .eq. 1) then
         pot = ders(0)
      endif
      
      if (ifgrad .eq. 1) then
         grad(1) = dx*ders(1)/r
         grad(2) = dy*ders(1)/r
      endif

      if (ifhess .eq. 1) then
         hess(1) = dx2*ders(2)/r2+ders(1)*(1.0d0/r-dx2/r3)
         hess(2) = dx*dy*(ders(2)/r2-ders(1)/r3)
         hess(3) = dy2*ders(2)/r2+ders(1)*(1.0d0/r-dy2/r3)
      endif

      if (ifder3 .eq. 1) then
         der3(1) = (dx3*ders(3)+3.0d0*dy2*dx*(ders(2)/r-ders(1)/r2))/r3
         der3(2) = dx2*dy*(ders(3)/r3-3.0d0*(ders(2)/r4-ders(1)/r5)) + 
     1        dy*(ders(2)/r2-ders(1)/r3)
         der3(3) = dx*dy2*(ders(3)/r3-3.0d0*(ders(2)/r4-ders(1)/r5)) + 
     1        dx*(ders(2)/r2-ders(1)/r3)
         der3(4) = (dy3*ders(3)+3.0d0*dx2*dy*(ders(2)/r-ders(1)/r2))/r3
      endif

      if (ifder4 .eq. 1) then
         der4(1) = (dx4*(ders(4)-6.0d0*ders(3)/r+15.0d0*(ders(2)/r2
     1        -ders(1)/r3)))/r4 +
     1        (dx2*6.0d0*(ders(3)-3.0d0*(ders(2)/r-ders(1)/r2)))/r3 +
     2        (3.0d0*(ders(2)-ders(1)/r))/r2
         der4(2) = (dx3*dy*(ders(4)-6.0d0*ders(3)/r+15.0d0*(ders(2)/r2
     1        -ders(1)/r3)))/r4 + 
     1        (3.0d0*dx*dy*(ders(3)-3.0d0*(ders(2)/r-ders(1)/r2)))/r3
         der4(3) = dx2*dy2*(ders(4)-6.0d0*ders(3)/r+15.0d0*ders(2)/r2
     1        -15.0d0*ders(1)/r3)/r4
     1        + ders(3)/r - 2.0d0*ders(2)/r2 + 2.0d0*ders(1)/r3
         der4(4) = dx*dy3*(ders(4)-6.0d0*ders(3)/r+15.0d0*(ders(2)/r2
     1        -ders(1)/r3))/r4 + 
     1        3.0d0*dx*dy*(ders(3)-3.0d0*(ders(2)/r-ders(1)/r2))/r3
         der4(5) = dy4*(ders(4)-6.0d0*ders(3)/r+15.0d0*(ders(2)/r2
     1        -ders(1)/r3))/r4 +
     1        dy2*6.0d0*(ders(3)-3.0d0*(ders(2)/r-ders(1)/r2))/r3 +
     2        3.0d0*(ders(2)-ders(1)/r)/r2        
      endif

      if (ifder5 .eq. 1) then
         der5(1) = (dx5*ders(5)+10.0d0*dy2*dx3*ders(4)/r +
     1        (15.0d0*dy4*dx-30.0d0*dy2*dx3)*ders(3)/r2 +
     2        (60.0d0*dy2*dx3-45.0d0*dy4*dx)*ders(2)/r3 +
     3        (45.0d0*dy4*dx-60.0d0*dy2*dx3)*ders(1)/r4)/r5
         der5(2) = (dy*dx4*ders(5)
     1        +(6.0d0*dy3*dx2-4.0d0*dy*dx4)*ders(4)/r +
     1        (3.0d0*dy5+12.0d0*dy*dx4-30.0d0*dy3*dx2)*ders(3)/r2 +
     2        (72.0d0*dy3*dx2-9.0d0*dy5-24.0d0*dy*dx4)*ders(2)/r3 +
     3        (9.0d0*dy5-72.0d0*dy3*dx2+24.0d0*dy*dx4)*ders(1)/r4)/r5
         der5(3) = (dy2*dx3*ders(5)
     1        +(3.0d0*dy4*dx-6.0d0*dy2*dx3+dx5)*ders(4)/r +
     1        (27.0d0*dy2*dx3-15.0d0*dy4*dx-3.0d0*dx5)*ders(3)/r2 +
     2        (36.0d0*dy4*dx-63.0d0*dy2*dx3+6.0d0*dx5)*ders(2)/r3 +
     3        (63.0d0*dy2*dx3-36.0d0*dy4*dx-6.0d0*dx5)*ders(1)/r4)/r5
         der5(4) = (dx2*dy3*ders(5)
     1        +(3.0d0*dx4*dy-6.0d0*dx2*dy3+dy5)*ders(4)/r +
     1        (27.0d0*dx2*dy3-15.0d0*dx4*dy-3.0d0*dy5)*ders(3)/r2 +
     2        (36.0d0*dx4*dy-63.0d0*dx2*dy3+6.0d0*dy5)*ders(2)/r3 +
     3        (63.0d0*dx2*dy3-36.0d0*dx4*dy-6.0d0*dy5)*ders(1)/r4)/r5
         der5(5) = (dx*dy4*ders(5)
     1        +(6.0d0*dx3*dy2-4.0d0*dx*dy4)*ders(4)/r +
     1        (3.0d0*dx5+12.0d0*dx*dy4-30.0d0*dx3*dy2)*ders(3)/r2 +
     2        (72.0d0*dx3*dy2-9.0d0*dx5-24.0d0*dx*dy4)*ders(2)/r3 +
     3        (9.0d0*dx5-72.0d0*dx3*dy2+24.0d0*dx*dy4)*ders(1)/r4)/r5
         der5(6) = (dy5*ders(5)+10.0d0*dx2*dy3*ders(4)/r +
     1        (15.0d0*dx4*dy-30.0d0*dx2*dy3)*ders(3)/r2 +
     2        (60.0d0*dx2*dy3-45.0d0*dx4*dy)*ders(2)/r3 +
     3        (45.0d0*dx4*dy-60.0d0*dx2*dy3)*ders(1)/r4)/r5
      endif


      return
      end

      subroutine helmbhgreen(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green's function and its derivatives stably
c     (with the derivatives on the "x" or "target" variable).
c
c     INPUTS
c
c     zk - complex *16, helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - complex *16, potential, if requested
c     grad(2) - complex *16, gradient, if requested
c     hess(3) - complex *16, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 zx(2), zy(2)
      complex *16 zk, pot, grad(2), hess(3)
      integer ifpot, ifgrad, ifhess
c     local variables
      complex *16 der3(4), der4(5), der5(6)
      integer ifder3, ifder4, ifder5

      ifder3 = 0
      ifder4 = 0
      ifder5 = 0
      
      call helmbhgreen_all(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)

      return
      end

      subroutine helmbhgreend1(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates a directional derivative of the Green's function 
c     at the source ("y" variable) and its derivatives stably (with
c     the derivatives on the "x" or "target" variable).
c
c     The potential is (nabla_y * dir1) G
c     The gradient is grad_x (nabla_y * dir1) G
c     The hessian is hess_x (nabla_y * dir1) G
c
c     INPUTS
c
c     zk - complex *16, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     dir1(2) - real *8, vector determining direction of derivative
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - complex *16, potential, if requested
c     grad(2) - complex *16, gradient, if requested
c     hess(3) - complex *16, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, pot, grad(2), hess(3)
      real *8 zx(2), zy(2), dir1(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      complex *16 potloc, gradloc(2), hessloc(3), der3(4), der4(5)
      complex *16 der5(6)
      integer ifder3, ifder4, ifder5, ifpotloc, ifgradloc, ifhessloc

      ifpotloc = 0
      ifgradloc = 0
      ifhessloc = 0
      ifder3 = 0
      ifder4 = 0
      ifder5 = 0

      if (ifpot .eq. 1) ifgradloc = 1
      if (ifgrad .eq. 1) ifhessloc = 1
      if (ifhess .eq. 1) ifder3 = 1

      call helmbhgreen_all(zk,zx,zy,ifpotloc,potloc,ifgradloc,gradloc,
     1     ifhessloc,hessloc,ifder3,der3,ifder4,der4,ifder5,der5)

      if (ifpot .eq. 1) then
         pot = -dir1(1)*gradloc(1)-dir1(2)*gradloc(2)
      endif

      if (ifgrad .eq. 1) then
         grad(1) = -dir1(1)*hessloc(1) - dir1(2)*hessloc(2)
         grad(2) = -dir1(1)*hessloc(2) - dir1(2)*hessloc(3)
      endif

      if (ifhess .eq. 1) then
         hess(1) = -dir1(1)*der3(1) - dir1(2)*der3(2)
         hess(2) = -dir1(1)*der3(2) - dir1(2)*der3(3)
         hess(3) = -dir1(1)*der3(3) - dir1(2)*der3(4)
      endif
      
      return
      end

      subroutine helmbhgreend2(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1,dir2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates double directional derivatives (on y) of the Green's
c     function and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     zk - complex *16, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     dir1(2), dir2(2) - real *8, vectors determining directions of
c                        derivatives
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - complex *16, potential, if requested
c     grad(2) - complex *16, gradient, if requested
c     hess(3) - complex *16, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, pot, grad(2), hess(3)
      real *8 zx(2), zy(2), dir1(2), dir2(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 quadvec(3)
      integer if0, if1, if2, if3
      integer ifpot1,ifgrad1,ifhess1,ifder3,ifder4,ifder5
      complex *16 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)

c     convert two directional derivatives to xx, xy, yy

      quadvec(1) = dir1(1)*dir2(1)
      quadvec(2) = dir1(1)*dir2(2) + dir1(2)*dir2(1)
      quadvec(3) = dir1(2)*dir2(2)

      ifpot1 = 0
      ifgrad1 = 0
      ifhess1 = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

      call helmbhgreen_all(zk,zx,zy,ifpot1,potloc,ifgrad1,gradloc,
     1     ifhess1,hessloc,ifder3,der3,ifder4,der4,ifder5,der5)

      if (ifpot .eq. 1) then
         pot = quadvec(1)*hessloc(1)+quadvec(2)*hessloc(2)
     1        +quadvec(3)*hessloc(3)
      endif

      if (ifgrad .eq. 1) then
         grad(1) = quadvec(1)*der3(1)+quadvec(2)*der3(2)
     1        +quadvec(3)*der3(3)
         grad(2) = quadvec(1)*der3(2)+quadvec(2)*der3(3)
     1        +quadvec(3)*der3(4)
      endif

      if (ifhess .eq. 1) then
         hess(1) = quadvec(1)*der4(1)+quadvec(2)*der4(2)
     1        +quadvec(3)*der4(3)
         hess(2) = quadvec(1)*der4(2)+quadvec(2)*der4(3)
     1        +quadvec(3)*der4(4)
         hess(3) = quadvec(1)*der4(3)+quadvec(2)*der4(4)
     1        +quadvec(3)*der4(5)
      endif

      return
      end

      subroutine diffsloghank(r,zk,rscale,diffs,ifders,ders,
     1     hvec,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     WARNING: DERS NOT IMPLEMENTED YET
c
c     This subroutine evaluates the difference kernel
c
c     diffs(0) = H_0(zk * r) - 2*i*log(r)/pi 
c
c     diffs(1) = (H_n(zk*r)+2*i*2^(n-1)*(n-1)!/(pi*(zk*r)^n))*rscale^n
c
c     Where H_n is the type 1 Hankel function
c
c     For abs(zk*r/2) < 1.0, this routine performs the 
c     cancellation between the two parts explicitly.
c     Otherwise, the formula is evaluated directly
c
c     input:
c     r - real *8, as above
c     zk - complex *16, as above
c     rscale - real *8, scaling factor, as above
c     ifders - integer, flag: determines whether or not to compute 
c              derivatives
c     n - integer, order: compute the 0 through n values of diffs
c         and ders (if requested)
c     
c     output:
c     diffs - complex *16 array, as above
c     ders - complex *16 array, derivative with respect to r
c     hvec - complex *16 array, can be viewed as workspace, length
c            n+1, returns Hankel Function values (without differencing)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, diffs(0:n), ders(0:n), hvec(0:n)
      real *8 r, rscale
      integer n, ifders
c     local variables
      complex *16 eye
      integer ifder, p, j
      data eye / (0.0d0,1.0d0) /
      integer pmax
      parameter (pmax = 12)
      real *8 logcs0(0:pmax), logcs1(0:pmax), rscale2
      complex *16 cs0(0:pmax), cs1(0:pmax)
      complex *16 zkr, zkrh, zkrh2p, hder(0:1), z2iopi, ztemp, zlogtemp
      complex *16 ztemp0, ztemp1, ztempl0, ztempl1, zkrhsc, zkrh2, zkh
      complex *16 diffstop, diffs2
      data z2iopi / (0.0d0,0.636619772367581343075535053d0) /

c     H0-2*i*log(zk*r)/pi is sum of two power series (for small arg)

c     Note that we have to correct and add back 2*i*log(zk)/pi to get
c     true function...

c     coefficients of log(zk*r/2)*(zk*r/2)^(2*j) 
c     (multiplied later by 2*i/pi)

      data logcs0 /   0.000000000000000000D+000,
     1     -0.100000000000000000D+001,
     2     0.250000000000000000D+000,
     3     -0.277777777777777778D-001,
     4     0.173611111111111111D-002,
     5     -0.694444444444444444D-004,
     6     0.192901234567901235D-005,
     7     -0.393675988914084152D-007,
     8     0.615118732678256488D-009,
     9     -0.759405842812662331D-011,
     $     0.759405842812662331D-013,
     1     -0.627608134555919281D-015,
     2     0.435838982330499501D-017 /

c     coefficients of (zk*r/2)^(2*j)

      data cs0 / (0.100000000000000000D+001,0.367466905196615962D+000),
     1     ( -0.100000000000000000D+001 ,   0.269152867170965382D+000),
     2     (  0.250000000000000000D+000 ,  -0.146865688338689013D+000),
     3     ( -0.277777777777777778D-001 ,   0.222130373373319398D-001),
     4     (  0.173611111111111111D-002 ,  -0.166462549867334231D-002),
     5     ( -0.694444444444444444D-004 ,   0.754269612298167666D-004),
     6     (  0.192901234567901235D-005 ,  -0.229986793422831468D-005),
     7     ( -0.393675988914084152D-007 ,   0.505163934110747222D-007),
     8     (  0.615118732678256488D-009 ,  -0.838268240495125880D-009),
     9     ( -0.759405842812662331D-011 ,   0.108861603731588473D-010),
     $     (  0.759405842812662331D-013 ,  -0.113696131479448557D-012),
     1     ( -0.627608134555919281D-015 ,   0.975959972766743018D-015),
     2     (  0.435838982330499501D-017 ,  -0.700871957231362727D-017) /

c     Likewise, H1+i*2/(zk*r*pi) is sum of two power series
c     for small arg

c     coefficients of log(zk*r/2)*(zk*r/2)^(2*j+1)
c     (multiplied later by 2*i/pi)

      data logcs1 /   0.100000000000000000D+001,
     1     -0.500000000000000000D+000,
     2     0.833333333333333333D-001,
     3     -0.694444444444444444D-002,
     4     0.347222222222222222D-003,
     5     -0.115740740740740741D-004,
     6     0.275573192239858907D-006,
     7     -0.492094986142605190D-008,
     8     0.683465258531396098D-010,
     9     -0.759405842812662331D-012,
     $     0.690368948011511210D-014,
     1     -0.523006778796599401D-016,
     2     0.335260755638845770D-018 /

c     coefficients of (zk*r/2)^(2*j+1)

      data cs1 /(0.100000000000000000D+001,0.491570190128252900D-001),
     1     ( -0.500000000000000000D+000 ,   0.214153905131430359D+000),
     2     (  0.833333333333333333D-001 ,  -0.577971707291127453D-001),
     3     ( -0.694444444444444444D-002 ,   0.610588066451317710D-002),
     4     (  0.347222222222222222D-003 ,  -0.355029952941876147D-003),
     5     ( -0.115740740740740741D-004 ,   0.131851839051696746D-004),
     6     (  0.275573192239858907D-006 ,  -0.341083657955069719D-006),
     7     ( -0.492094986142605190D-008 ,   0.651034755017267366D-008),
     8     (  0.683465258531396098D-010 ,  -0.955581794844995840D-010),
     9     ( -0.759405842812662331D-012 ,   0.111278867605518515D-011),
     $     (  0.690368948011511210D-014 ,  -0.105357858265556574D-013),
     1     ( -0.523006778796599401D-016 ,   0.827173162991627227D-016),
     2     (  0.335260755638845770D-018 ,  -0.547341260406378455D-018) /


      rscale2 = rscale**2
      zkr = zk*r
      zkrh = 0.5d0*zkr
      zkrh2 = zkrh**2
      zkrhsc = rscale/zkrh

c     get Hankel function values 

      ifder = 0
      call h2dall(n,zkr,rscale,hvec,ifder,hder)

      if (abs(zkr) .le. 1.0d0) then

c     small arg, use power series to evaluate stably...
c     note: there's some trickery in here around the log terms

         ztemp0 = 0.0d0
         ztemp1 = 0.0d0
         ztempl0 = 0.0d0
         ztempl1 = 0.0d0
         zkrh2p = 1.0d0

c     note: can remove some of the redundancy here later...
         
         do p = 0,pmax
            ztemp0 = ztemp0 + cs0(p)*zkrh2p
            ztempl0 = ztempl0 + logcs0(p)*zkrh2p
            ztemp1 = ztemp1 + cs1(p)*zkrh2p
            ztempl1 = ztempl1 + logcs1(p)*zkrh2p
            zkrh2p = zkrh2p*zkrh2
         enddo

         zlogtemp = cdlog(zk/2.0d0)*z2iopi
         
         ztemp = cdlog(zkrh)*z2iopi
         diffs(0) = ztemp0 + zlogtemp + ztempl0*ztemp
         diffs(1) = ztemp1*zkrh + ztempl1*ztemp*zkrh

         diffs(1) = diffs(1)*rscale

         diffs2 = 0.0d0

         if (n .ge. 2 .or. ifders .eq. 1) then
c     recursion formula is bad for diffs(2)
c     explicitly cancel the log terms, reusing the above
            diffs2 = (ztemp1-ztemp0
     1           +(ztempl1-ztempl0-1.0d0)*ztemp)*rscale2
         endif
c     simple recursion for rest
         if (n .ge. 2) then
            diffs(2) = diffs2
            ztemp = diffs(2)
c     recursion formula pretty good for the rest
            do j = 2,n-1
               ztemp = (j*ztemp)*zkrhsc
     1              - hvec(j-1)*rscale2
               diffs(j+1) = ztemp
            enddo
         endif

         if (ifders .eq. 1) then
            
            ders(0) = -diffs(1)*zk/rscale
            
            zkh = 0.5d0*zk
            do j = 1,n-1
               ders(j) = zkh*(hvec(j-1)*rscale-diffs(j+1)/rscale)
            enddo

c     use recursion formula for diffs(n+1) implicitly
            ders(n) = zkh*(2*hvec(n-1)*rscale-diffs(n)*n/zkrh)
c     fix if n = 1 (recursion formula for diffs(2) is unstable in this regime)
            if (n .eq. 1) ders(1) = zkh*(2*hvec(0)*rscale-diffs2/rscale)
            
         endif

      else

c     directly evaluate

         diffs(0) = hvec(0) - dlog(r)*z2iopi

         ztemp = z2iopi/zkr
         diffs(1) = hvec(1) + ztemp*rscale
         do j = 1,n-1
            ztemp = (j*ztemp)*zkrhsc
            diffs(j+1) = hvec(j+1) + ztemp
         enddo

         if (ifders .eq. 1) then
            
            ders(0) = -diffs(1)*zk/rscale
            
            zkh = 0.5d0*zk
            do j = 1,n-1
               ders(j) = zkh*(hvec(j-1)*rscale-diffs(j+1)/rscale)
            enddo

c     use recursion formula for diffs(n+1) implicitly
            ders(n) = zkh*(2*hvec(n-1)*rscale-diffs(n)*n/zkrh)
            
         endif
         
      endif


      return
      end

      subroutine diffszkjfun(r,zk,rscale,diffs,ifders,ders,jvec,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     J_n(zk*r) - (zk*r)^n/(rscale^2*2^n*n!)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, diffs(0:n), ders(0:n), jvec(0:n)
      real *8 r, rscale
      integer n, ifders
c     local variables
      complex *16 eye
      integer ifder, p, j, ntop, ier, n1
      data eye / (0.0d0,1.0d0) /
      integer pmax, lwfjs
      parameter (pmax = 12)
      real *8 cs0(0:pmax), cs1(0:pmax), rscale2
      complex *16 zkr, zkrh, zkr2p, jder(0:1), ztemp
      complex *16, allocatable :: jvals(:)
      integer, allocatable :: iscale(:)

c     J0-1 and J1-zk*r/2 are given as power series
c     coefficients are precomputed

c     J0-1, coefficients of (zk*r)^(2*j)

      data cs0 /   0.000000000000000000D+000,
     1     -0.250000000000000000D+000,
     2     0.156250000000000000D-001,
     3     -0.434027777777777778D-003,
     4     0.678168402777777778D-005,
     5     -0.678168402777777778D-007,
     6     0.470950279706790123D-009,
     7     -0.240280754952443941D-011,
     8     0.938596699032984143D-014,
     9     -0.289690339207711155D-016,
     $     0.724225848019277888D-019,
     1     -0.149633439673404522D-021,
     2     0.259780277210771740D-024 /

c     J1-zk*r/2, coefficients of (zk*r)^(2*j+1)
      
      data cs1 /   0.000000000000000000D+000,
     1     -0.625000000000000000D-001,
     2     0.260416666666666667D-002,
     3     -0.542534722222222222D-004,
     4     0.678168402777777778D-006,
     5     -0.565140335648148148D-008,
     6     0.336393056933421517D-010,
     7     -0.150175471845277463D-012,
     8     0.521442610573880079D-015,
     9     -0.144845169603855578D-017,
     $     0.329193567281489949D-020,
     1     -0.623472665305852176D-023,
     2     0.999154912349122077D-026 /

      rscale2 = rscale**2
      zkr = zk*r
      zkrh = 0.5d0*zkr

c     get Bessel function values 

      lwfjs = n+10+4*n+100
      allocate(jvals(0:lwfjs),iscale(lwfjs))
      n1 = n+1
      ifder = 0
      call jfuns2d(ier,n1,zkr,rscale,jvals,ifder,jder,
     1     lwfjs,iscale,ntop)
      if (ier .ne. 0) then
         write(*,*) "WARNING, JFUNS BOMBED!!!", ier, zkr
         return
      endif
      do j = 0,n
         jvec(j) = jvals(j)
      enddo

      if (abs(zkr) .le. 1.0d0) then

c     small arg, use power series to evaluate stably...

         diffs(0) = 0.0d0
         diffs(1) = 0.0d0
         zkr2p = zkr**2

         do p = 1,pmax
            diffs(0) = diffs(0) + cs0(p)*zkr2p
            zkr2p = zkr2p*zkr
            diffs(1) = diffs(1) + cs1(p)*zkr2p
            zkr2p = zkr2p*zkr
         enddo

c     simple recursion for rest

         diffs(1) = diffs(1)/rscale
         ztemp = diffs(1)
         do j = 2,n
            ztemp = (zkrh*(ztemp+jvals(j+1)*rscale2))/(rscale*j)
            diffs(j) = ztemp
         enddo

         
      else

c     directly evaluate

         ztemp = 1.0d0
         diffs(0) = jvals(0) - ztemp
         do j = 1,n
            ztemp = (zkrh*ztemp)/(j*rscale)
            diffs(j) = jvals(j) - ztemp
         enddo
      endif


      if (ifders .eq. 1) then
         
      endif

      return
      end
