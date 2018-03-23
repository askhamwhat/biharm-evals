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

      subroutine hbhpotgrad2dall_cdq(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

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

      do i = 1,ns
         call helmbhgreen_all(beta,target,source(1,i),ifpot,potloc,
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

      subroutine hbhpotgrad2dall_cdq_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3)

      call hbhpotgrad2dall_cdq(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pottemp,ifgrad,gradtemp,ifhess,hesstemp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif

      return
      end
  
      subroutine hbhpotgrad2dall_cdq3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3loc(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3loc, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3loc = 1
      ifder4 = 1
      ifder5 = 1

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
      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0
      endif

      do i = 1,ns
         call helmbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3loc,der3loc,
     2        ifder4,der4,ifder5,der5)

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
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) + der3loc(1)*charge(i)
               der3(2) = der3(2) + der3loc(2)*charge(i)
               der3(3) = der3(3) + der3loc(3)*charge(i)
               der3(4) = der3(4) + der3loc(4)*charge(i)
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
               hess(1) = hess(1) - dipstr(i)*(der3loc(1)*dipvec(1,i)
     1              + der3loc(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3loc(2)*dipvec(1,i)
     1              + der3loc(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3loc(3)*dipvec(1,i)
     1              + der3loc(4)*dipvec(2,i))
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) - dipstr(i)*(der4(1)*dipvec(1,i)
     1              + der4(2)*dipvec(2,i))
               der3(2) = der3(2) - dipstr(i)*(der4(2)*dipvec(1,i)
     1              + der4(3)*dipvec(2,i))
               der3(3) = der3(3) - dipstr(i)*(der4(3)*dipvec(1,i)
     1              + der4(4)*dipvec(2,i))
               der3(4) = der3(4) - dipstr(i)*(der4(4)*dipvec(1,i)
     1              + der4(5)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3loc(1)*quadvec(1,i)
     1              + der3loc(2)*quadvec(2,i) + der3loc(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3loc(2)*quadvec(1,i)
     1              + der3loc(3)*quadvec(2,i) + der3loc(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) + quadstr(i)*(der5(1)*quadvec(1,i)
     1              + der5(2)*quadvec(2,i) + der5(3)*quadvec(3,i))
               der3(2) = der3(2) + quadstr(i)*(der5(2)*quadvec(1,i)
     1              + der5(3)*quadvec(2,i) + der5(4)*quadvec(3,i))
               der3(3) = der3(3) + quadstr(i)*(der5(3)*quadvec(1,i)
     1              + der5(4)*quadvec(2,i) + der5(5)*quadvec(3,i))
               der3(4) = der3(4) + quadstr(i)*(der5(4)*quadvec(1,i)
     1              + der5(5)*quadvec(2,i) + der5(6)*quadvec(3,i))
            endif
         endif
      enddo
      
      return
      end

      subroutine hbhpotgrad2dall_cdq3_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)

      call hbhpotgrad2dall_cdq3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pottemp,ifgrad,gradtemp,ifhess,hesstemp,
     3     ifder3,der3temp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif
      if (ifder3 .eq. 1) then
         der3(1) = der3(1) + der3temp(1)
         der3(2) = der3(2) + der3temp(2)
         der3(3) = der3(3) + der3temp(3)
         der3(4) = der3(4) + der3temp(4)
      endif

      return
      end
  
      subroutine hbhpotgrad2dall_cdqo(beta,source,ns,ifcharge,
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
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 octstr(*), octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 1

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

      do i = 1,ns
         call helmbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

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

      subroutine hbhpotgrad2dall_cdqo_add(beta,source,ns,ifcharge,
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
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3)

      call hbhpotgrad2dall_cdqo(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pottemp,ifgrad,gradtemp,
     3     ifhess,hesstemp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif

      return
      end
  
      subroutine hbhpotgrad2dall_cdqo3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      integer i, ntermstemp
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)
      real *8 rscale, rs2, pi, xtemp, ytemp, rtemp, tiny
      complex *16 ima, hexphbh(0:3), hexpy(0:3)
      data ima /(0.0d0,1.0d0)/
      data tiny / 1.0d-200 /

      pi = 4.0d0*datan(1.0d0)

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
      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0
      endif

      ntermstemp = 0
      if (ifdipole .eq. 1) ntermstemp = 1
      if (ifquad .eq. 1) ntermstemp = 2
      if (ifoct .eq. 1) ntermstemp = 3

      do i = 1,ns
         
         xtemp = target(1)-source(1,i)
         ytemp = target(2)-source(2,i)

         rtemp = dsqrt(xtemp**2+ytemp**2)

         rscale = min(rtemp*beta,1.0d0)
         if (rscale .lt. tiny) return

         rs2 = rscale*rscale

         hexphbh(0) = 0.0d0
         hexphbh(1) = 0.0d0
         hexphbh(2) = 0.0d0
         hexphbh(3) = 0.0d0

         hexpy(0) = 0.0d0
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0
         hexpy(3) = 0.0d0

         if (ifcharge .eq. 1) then
            hexphbh(0) = hexphbh(0) + charge(i)/(2.0d0*pi*beta**2)
         endif

         if (ifdipole .eq. 1) then
            hexphbh(1) = hexphbh(1) + dipstr(i)*(dipvec(1,i)+
     1           ima*dipvec(2,i))/(rscale*2.0d0*pi*beta)
         endif

         if (ifquad .eq. 1) then
            hexphbh(2) = hexphbh(2) + quadstr(i)*(quadvec(1,i)
     1           +quadvec(2,i)*ima
     2           -quadvec(3,i))/(rs2*4.0d0*pi)

            hexpy(0) = quadstr(i)*(quadvec(1,i)
     1           +quadvec(3,i))/(4.0d0*pi)
         endif

         if (ifoct .eq. 1) then
            hexphbh(3) = hexphbh(3) + beta*octstr(i)*(octvec(1,i) 
     1           + ima*octvec(2,i) - octvec(3,i) 
     2           -ima*octvec(4,i))/(rs2*rscale*8.0d0*pi)
            hexpy(1) = hexpy(1) + beta*octstr(i)*(3.0d0*octvec(1,i)
     1           + ima*octvec(2,i) + octvec(3,i) 
     2           + 3.0d0*ima*octvec(4,i))/(rscale*8.0d0*pi)
         endif

         call hbh2dmpeval3(beta,rscale,source(1,i),hexphbh,hexpy,
     1        ntermstemp,target,pottemp,ifgrad,gradtemp,ifhess,hesstemp,
     2        ifder3,der3temp)

         if (ifpot .eq. 1) pot = pot+pottemp
         if (ifgrad .eq. 1) then
            grad(1) = grad(1)+gradtemp(1)
            grad(2) = grad(2)+gradtemp(2)
         endif
         if (ifhess .eq. 1) then
            hess(1) = hess(1)+hesstemp(1)
            hess(2) = hess(2)+hesstemp(2)
            hess(3) = hess(3)+hesstemp(3)
         endif
         if (ifder3 .eq. 1) then
            der3(1) = der3(1)+der3temp(1)
            der3(2) = der3(2)+der3temp(2)
            der3(3) = der3(3)+der3temp(3)
            der3(4) = der3(4)+der3temp(4)
         endif
         
      enddo
      
      return
      end

      subroutine hbhpotgrad2dall_cdqo3_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 target(2), octvec(4,*)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)

      call hbhpotgrad2dall_cdqo3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pottemp,ifgrad,gradtemp,
     3     ifhess,hesstemp,ifder3,der3temp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif
      if (ifder3 .eq. 1) then
         der3(1) = der3(1) + der3temp(1)
         der3(2) = der3(2) + der3temp(2)
         der3(3) = der3(3) + der3temp(3)
         der3(4) = der3(4) + der3temp(4)
      endif

      return
      end
  
