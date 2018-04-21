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

      subroutine hbhpotgrad2dall_cdq(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk
      real *8 source(2,*)
      complex *16 charge(*), dipstr(*)
      complex *16 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      complex *16 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif

      call hbhpotgrad2dall_cdq_add(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)

      

      return
      end

      subroutine hbhpotgrad2dall_cdq_add(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk
      real *8 source(2,*)
      complex *16 charge(*), dipstr(*)
      complex *16 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      complex *16 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
c     local variables
      complex *16 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

      do i = 1,ns
         call helmbhgreen_all(zk,target,source(1,i),ifpot,potloc,
     1        ifgrad,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,der4,
     2        ifder5,der5)

         if (ifcharge .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + potloc*charge(i)
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradloc(1)*charge(i)
               grad(2) = grad(2) + gradloc(2)*charge(i)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hessloc(1)*charge(i)
               hess(2) = hess(2) + hessloc(2)*charge(i)
               hess(3) = hess(3) + hessloc(3)*charge(i)
            endif
         endif

         if (ifdipole .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1              + gradloc(2)*dipvec(2,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1              + hessloc(2)*dipvec(2,i))
               grad(2) = grad(2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1              + hessloc(3)*dipvec(2,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1              + der3(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1              + der3(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1              + der3(4)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3(1)*quadvec(1,i)
     1              + der3(2)*quadvec(2,i) + der3(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3(2)*quadvec(1,i)
     1              + der3(3)*quadvec(2,i) + der3(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
         endif
      enddo

      
      return
      end
  
  
      subroutine hbhpotgrad2dall_cdqo(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, charge(*), dipstr(*), quadstr(*), octstr(*)
      real *8 source(2,*)
      real *8 dipvec(2,*), quadvec(3,*)
      real *8 octvec(4,*)
      real *8 target(2)
      complex *16 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif

      call hbhpotgrad2dall_cdqo_add(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess)

      
      return
      end

      subroutine hbhpotgrad2dall_cdqo_add(zk,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 zk, charge(*), dipstr(*), quadstr(*), octstr(*)
      real *8 source(2,*)
      real *8 dipvec(2,*), quadvec(3,*)
      real *8 octvec(4,*)
      real *8 target(2)
      complex *16 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct
c     local variables
      complex *16 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 1

      do i = 1,ns
         call helmbhgreen_all(zk,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

c         write(*,*) i
c         write(*,*) charge(i), dipstr(i), quadstr(i), octstr(i)
c         write(*,*) octvec(1,i), octvec(2,i), octvec(3,i), octvec(4,i)


         if (ifcharge .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + potloc*charge(i)
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradloc(1)*charge(i)
               grad(2) = grad(2) + gradloc(2)*charge(i)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hessloc(1)*charge(i)
               hess(2) = hess(2) + hessloc(2)*charge(i)
               hess(3) = hess(3) + hessloc(3)*charge(i)
            endif
         endif

         if (ifdipole .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1              + gradloc(2)*dipvec(2,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1              + hessloc(2)*dipvec(2,i))
               grad(2) = grad(2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1              + hessloc(3)*dipvec(2,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1              + der3(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1              + der3(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1              + der3(4)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3(1)*quadvec(1,i)
     1              + der3(2)*quadvec(2,i) + der3(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3(2)*quadvec(1,i)
     1              + der3(3)*quadvec(2,i) + der3(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
         endif

         if (ifoct .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - octstr(i)*(der3(1)*octvec(1,i)
     1              + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1              + der3(4)*octvec(4,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - octstr(i)*(der4(1)*octvec(1,i)
     1              + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1              + der4(4)*octvec(4,i))
               grad(2) = grad(2) - octstr(i)*(der4(2)*octvec(1,i)
     1              + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1              + der4(5)*octvec(4,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - octstr(i)*(der5(1)*octvec(1,i)
     1              + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1              + der5(4)*octvec(4,i))
               hess(2) = hess(2) - octstr(i)*(der5(2)*octvec(1,i)
     1              + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1              + der5(5)*octvec(4,i))
               hess(3) = hess(3) - octstr(i)*(der5(3)*octvec(1,i)
     1              + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1              + der5(6)*octvec(4,i))
            endif
         endif

      enddo

      return
      end
