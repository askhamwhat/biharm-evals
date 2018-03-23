cc Copyright (C) 2017: Travis Askham
cc Contact: askhamwhat@gmail.com
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
      
      subroutine helmbhgreen_all(zk,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green function 
c
c     G(x,y) = (S_zk(r) - S_0(r))/zk^2 ,
c
c     where S_zk is the fundamental solution of the Yukawa equation
c     and S_0 is the fundamental solution of the Laplace equation and 
c     r = sqrt( (x_1-y_1)^2 + (x_2-y_2)^2),
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     zk - real *8, modified helmholtz parameter
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
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c     der3(4) - real *8, 3rd order derivatives, if requested.
c     der4(5) - real *8, 4th order derivatives, if requested.
c     der5(6) - real *8, 5th order derivatives, if requested.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk
      real *8 zx(2), zy(2), pot, grad(2), hess(3)
      real *8 der3(4), der4(5), der5(6)
      integer ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
c     local variables
      real *8 dx, dy, r, r2, r3, r4, r5, r6
      real *8 dx2, dx3, dx4, dx5, dy2, dy3, dy4, dy5
      real *8 g0, g1, g2, g3, g4, g5, dtemp
      real *8 g0temp, g1temp, g2temp, g3temp
      integer if1, if2, if3, if0, ifders
      real *8 rscale, pih
      integer nmax, ifder, j
      complex *16 z, eye, hvec(0:10), hder(0:10), diffs(0:5)
      complex *16 zkh      
      real *8 kvec(0:10), ders(0:5)
      data eye /(0.0d0,1.0d0)/

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

c     get values of difference kernel and higher modes
      
      rscale = 1.0d0
      ifders = 0
      call diffsloghank(r,zk,rscale,diffs,ifders,ders,kvec,nmax)

c     get values of modified Bessel functions

      ifder = 0
      rscale = 1.0d0


      pih = 2.0d0*atan(1.0d0)
      

c     combine to get values of derivatives of G 
c     with respect to r

      dtemp = 1.0d0/(4.0d0*pih*zk**2)
      zkh = 0.5d0*zk

      g0 = diffs(0)*dtemp
      dtemp = -dtemp*zk
      g1 = diffs(1)*dtemp
      dtemp = -dtemp*zkh
      g2 = (diffs(2)+hvec(0))*dtemp
      dtemp = -dtemp*zkh
      g3 = (diffs(3)+3.0d0*hvec(1))*dtemp
      dtemp = -dtemp*zkh
      g4 = (diffs(4)+3.0d0*hvec(0)+4.0d0*hvec(2))*dtemp
      dtemp = -dtemp*zkh
      g5 = (diffs(5)+10.0d0*hvec(1)+5.0d0*hvec(3))*dtemp
      
      if0 = 1
      if1 = 1
      if2 = 1
      if3 = 1

c     evaluate potential and derivatives

      if (ifpot .eq. 1) then
         pot = g0
      endif
      
      if (ifgrad .eq. 1) then
         grad(1) = dx*g1/r
         grad(2) = dy*g1/r
      endif

      if (ifhess .eq. 1) then
         hess(1) = dx2*g2/r2+g1*(1.0d0/r-dx2/r3)
         hess(2) = dx*dy*(g2/r2-g1/r3)
         hess(3) = dy2*g2/r2+g1*(1.0d0/r-dy2/r3)
      endif

      if (ifder3 .eq. 1) then
         der3(1) = (dx3*g3+3.0d0*dy2*dx*(g2/r-g1/r2))/r3
         der3(2) = dx2*dy*(g3/r3-3.0d0*(g2/r4-g1/r5)) + 
     1        dy*(g2/r2-g1/r3)
         der3(3) = dx*dy2*(g3/r3-3.0d0*(g2/r4-g1/r5)) + 
     1        dx*(g2/r2-g1/r3)
         der3(4) = (dy3*g3+3.0d0*dx2*dy*(g2/r-g1/r2))/r3
      endif

      if (ifder4 .eq. 1) then
         der4(1) = (dx4*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3)))/r4 +
     1        (dx2*6.0d0*(g3-3.0d0*(g2/r-g1/r2)))/r3 +
     2        (3.0d0*(g2-g1/r))/r2
         der4(2) = (dx3*dy*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3)))/r4 + 
     1        (3.0d0*dx*dy*(g3-3.0d0*(g2/r-g1/r2)))/r3
         der4(3) = dx2*dy2*(g4-6.0d0*g3/r+15.0d0*g2/r2-15.0d0*g1/r3)/r4
     1        + g3/r - 2.0d0*g2/r2 + 2.0d0*g1/r3
         der4(4) = dx*dy3*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3))/r4 + 
     1        3.0d0*dx*dy*(g3-3.0d0*(g2/r-g1/r2))/r3
         der4(5) = dy4*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3))/r4 +
     1        dy2*6.0d0*(g3-3.0d0*(g2/r-g1/r2))/r3 +
     2        3.0d0*(g2-g1/r)/r2        
      endif

      if (ifder5 .eq. 1) then
         der5(1) = (dx5*g5+10.0d0*dy2*dx3*g4/r +
     1        (15.0d0*dy4*dx-30.0d0*dy2*dx3)*g3/r2 +
     2        (60.0d0*dy2*dx3-45.0d0*dy4*dx)*g2/r3 +
     3        (45.0d0*dy4*dx-60.0d0*dy2*dx3)*g1/r4)/r5
         der5(2) = (dy*dx4*g5+(6.0d0*dy3*dx2-4.0d0*dy*dx4)*g4/r +
     1        (3.0d0*dy5+12.0d0*dy*dx4-30.0d0*dy3*dx2)*g3/r2 +
     2        (72.0d0*dy3*dx2-9.0d0*dy5-24.0d0*dy*dx4)*g2/r3 +
     3        (9.0d0*dy5-72.0d0*dy3*dx2+24.0d0*dy*dx4)*g1/r4)/r5
         der5(3) = (dy2*dx3*g5+(3.0d0*dy4*dx-6.0d0*dy2*dx3+dx5)*g4/r +
     1        (27.0d0*dy2*dx3-15.0d0*dy4*dx-3.0d0*dx5)*g3/r2 +
     2        (36.0d0*dy4*dx-63.0d0*dy2*dx3+6.0d0*dx5)*g2/r3 +
     3        (63.0d0*dy2*dx3-36.0d0*dy4*dx-6.0d0*dx5)*g1/r4)/r5
         der5(4) = (dx2*dy3*g5+(3.0d0*dx4*dy-6.0d0*dx2*dy3+dy5)*g4/r +
     1        (27.0d0*dx2*dy3-15.0d0*dx4*dy-3.0d0*dy5)*g3/r2 +
     2        (36.0d0*dx4*dy-63.0d0*dx2*dy3+6.0d0*dy5)*g2/r3 +
     3        (63.0d0*dx2*dy3-36.0d0*dx4*dy-6.0d0*dy5)*g1/r4)/r5
         der5(5) = (dx*dy4*g5+(6.0d0*dx3*dy2-4.0d0*dx*dy4)*g4/r +
     1        (3.0d0*dx5+12.0d0*dx*dy4-30.0d0*dx3*dy2)*g3/r2 +
     2        (72.0d0*dx3*dy2-9.0d0*dx5-24.0d0*dx*dy4)*g2/r3 +
     3        (9.0d0*dx5-72.0d0*dx3*dy2+24.0d0*dx*dy4)*g1/r4)/r5
         der5(6) = (dy5*g5+10.0d0*dx2*dy3*g4/r +
     1        (15.0d0*dx4*dy-30.0d0*dx2*dy3)*g3/r2 +
     2        (60.0d0*dx2*dy3-45.0d0*dx4*dy)*g2/r3 +
     3        (45.0d0*dx4*dy-60.0d0*dx2*dy3)*g1/r4)/r5
      endif


      return
      end

      subroutine helmbhgreen(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green's function (S_beta - S_0)/beta^2 
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 der3(4), der4(5), der5(6)
      integer ifder3, ifder4, ifder5

      ifder3 = 0
      ifder4 = 0
      ifder5 = 0
      
      call helmbhgreen_all(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)

      return
      end

      subroutine helmbhgreend1(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates a directional derivative of the Green's function 
c     G = (S_beta - S_0)/beta^2 at the source ("y" variable)
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     The potential is (nabla_y * dir1) G
c     The gradient is grad_x (nabla_y * dir1) G
c     The hessian is hess_x (nabla_y * dir1) G
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3), dir1(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 potloc, gradloc(2), hessloc(3), der3(4), der4(5)
      real *8 der5(6)
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

      call helmbhgreen_all(beta,zx,zy,ifpotloc,potloc,ifgradloc,gradloc,
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

      subroutine helmbhgreend2(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1,dir2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green's function (S_beta - S_0)/beta^2 
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3), dir1(2)
      real *8 dir2(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 r, r2, dx, dy, g0, g1, g2, g3, temp1, dx2, dy2
      real *8 r3, pi, r4, r5
      real *8 g, gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy
      real *8 gxxxx, gxxxy, gxxyy, gxyyy, gyyyy
      real *8 quadvec(3)
      integer if0, if1, if2, if3
      integer ifpot1,ifgrad1,ifhess1,ifder3,ifder4,ifder5
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)

      pi = 4.0d0*datan(1.0d0)

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

      call helmbhgreen_all(beta,zx,zy,ifpot1,potloc,ifgrad1,gradloc,
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
      complex *16 zkr, zkrh, zkrh2p, hder(0:1), z2iopi, ztemp
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

      data cs0 / (0.100000000000000000D+001,-0.738042951086872253D-001),
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

c     get Hankel function values 

      ifder = 0
      call h2dall(n,zkr,rscale,hvec,ifder,hder)

      if (abs(zkrh) .le. 1.0d0) then

c     small arg, use power series to evaluate stably...
c     note: there's some trickery in here around the log terms

         diffs(0) = log(zk)*z2iopi
         diffs(1) = 0.0d0
         ztemp = log(zkrh)*z2iopi
         zkrh2p = 1.0d0

         do p = 0,pmax
            diffs(0) = diffs(0) + logcs0(p)*zkrh2p*ztemp
            diffs(0) = diffs(0) + cs0(p)*zkrh2p
            zkrh2p = zkrh2p*zkrh
            diffs(1) = diffs(1) + logcs1(p)*zkrh2p*ztemp
            diffs(1) = diffs(1) + cs1(p)*zkrh2p
            zkrh2p = zkrh2p*zkrh
         enddo

c     simple recursion for rest

         diffs(1) = diffs(1)*rscale
         ztemp = diffs(1)
         do j = 1,n-1
            ztemp = ((j*rscale)*ztemp)/zkrh
     1           - hvec(j-1)*rscale2
            diffs(j+1) = ztemp
         enddo

         
      else

c     directly evaluate

         diffs(0) = hvec(0) - dlog(r)*z2iopi

         ztemp = z2iopi/zkr
         diffs(1) = hvec(1) + ztemp*rscale
         do j = 1,n-1
            ztemp = ((j*rscale)*ztemp)/zkrh
            diffs(j+1) = hvec(j+1) + ztemp
         enddo
      endif


      if (ifders .eq. 1) then


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
