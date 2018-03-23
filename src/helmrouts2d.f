cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-04-25 22:45:52 -0400 (Mon, 25 Apr 2011) $
c    $Revision: 1889 $
c
c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	hvec(n)= H_n(z)*rscale^(n)
c       	jvec(n)= J_n(z)/rscale^(n)
c
c          rscale should be of the order of |z| if |z| < 1. Otherwise,
c          rscale should be set to 1.
c
c      2) The potential scaling is that required of the delta function
c      response:  pot = H_0(k*r)*(eye/4)
c
c-----------------------------------------------------------------------
c
c      H2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      H2DFORMMP: creates multipole expansion (outgoing) due to 
c                 a collection of charge sources.
c
c      H2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      H2DFORMTA: creates local expansion due to 
c                 a collection of charge sources.
c
c      H2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c
c      H2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c
c      H2DMPMP:     Converts multipole expansion to a multipole expansion.
c      H2DMPLOC:     Converts multipole expansion to a local expansion.
c      H2DLOCLOC:     Converts local expansion to a local expansion.
c
c      HPOTGRAD2DALL:  direct calculation for a collection of charge sources
c      HPOTGRAD2D : direct calculation for a single charge source
c
c
c
c      H2DFORMMP_DP: creates multipole expansion (outgoing) due to 
c                 a collection of dipole sources.
c
c      H2DFORMTA_DP: creates local expansion due to 
c                 a collection of dipole sources.
c
c      HPOTGRAD2DALL_DP:  direct calculation for a collection of dipoles
c      HPOTGRAD2D_DP : direct calculation for a single dipole
c
c
c      HPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      HPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c
c
c
c**********************************************************************
      subroutine h2dmpeval(zk,rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n H_n(k r) exp(i n theta)
c              n=-nterms  
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,pot,grad(2),hess(3),mpole(-nterms:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
c
        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
        call h2cart2polar(zdiff,r,theta)
c
        z=zk*r
        ifder=0
        call h2dall(nterms+3,z,rscale,hval,ifder,hder)
c
        allocate(mptemp(-nterms-2:nterms+2), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        mptemp(0)=hval(0)
        do n=1,nterms+2
        mptemp(n)=hval(n)*exp(ima*n*theta)
        mptemp(-n)=hval(n)*exp(ima*(-n)*theta)*(-1)**n
        enddo
c
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
c
           allocate(mpolex(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-1,nterms+1
              mpolex(i)=0
              mpoley(i)=0
           enddo
           do i=-nterms,nterms
             if (i .le. 0) then
             mpolex(i-1)=mpolex(i-1)+zk/2/rscale*mpole(i)
             mpoley(i-1)=mpoley(i-1)+zk/2*(ima)/rscale*mpole(i)
             endif
             if( i .gt. 0 ) then
                mpolex(i-1)=mpolex(i-1)+zk/2*rscale*mpole(i)
                mpoley(i-1)=mpoley(i-1)+zk/2*(ima)*rscale*mpole(i)
             endif
             if( i .ge. 0 ) then
                mpolex(i+1)=mpolex(i+1)-zk/2/rscale*mpole(i)
                mpoley(i+1)=mpoley(i+1)+zk/2*(ima)/rscale*mpole(i)
             endif
             if( i .lt. 0 ) then
                mpolex(i+1)=mpolex(i+1)-zk/2*rscale*mpole(i)
                mpoley(i+1)=mpoley(i+1)+zk/2*(ima)*rscale*mpole(i)
             endif
           enddo
        endif
c
        if( ifhess .eq. 1 ) then
           allocate(mpolexx(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpolexy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoleyy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-2,nterms+2
              mpolexx(i)=0
              mpolexy(i)=0
              mpoleyy(i)=0
           enddo
           do i=-nterms+1,nterms+1
              if( i .le. 0 ) then 
                 mpolexx(i-1)=mpolexx(i-1)+zk/2/rscale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+zk/2*(ima)/rscale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+zk/2*(ima)/rscale*mpoley(i)
              endif
              if( i .gt. 0 ) then
                 mpolexx(i-1)=mpolexx(i-1)+zk/2*rscale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+zk/2*(ima)*rscale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+zk/2*(ima)*rscale*mpoley(i)
              endif
              if( i .ge. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-zk/2/rscale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+zk/2*(ima)/rscale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+zk/2*(ima)/rscale*mpoley(i)
              endif
              if( i .lt. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-zk/2*rscale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+zk/2*(ima)*rscale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+zk/2*(ima)*rscale*mpoley(i)
              endif
           enddo
        endif
c
c
        pot=hval(0)*mpole(0)
        do n=1,nterms
           pot=pot+mpole(n)*mptemp(n)+mpole(-n)*mptemp(-n)
        enddo
        pot=pot/ima4
c
        if( ifgrad .eq. 1 ) then
           grad(1)=0
           grad(2)=0
c
           if( 1 .eq. 2 ) then
              rx=zdiff(1)
              ry=zdiff(2)
              ctheta=rx/r
              stheta=ry/r
              grad(1)=-hval(1)/rscale*zk*ctheta *mpole(0)
              grad(2)=-hval(1)/rscale*zk*stheta *mpole(0)
              do n=1,nterms
                 grad(1)=grad(1)+(-hval(n+1)/rscale*exp(ima*theta)+
     $              (exp(-ima*theta))*hval(n-1)*rscale)
     $              *(mpole(n)*exp(ima*n*theta))*zk/2
                 grad(1)=grad(1)+(-hval(n+1)/rscale*exp(-ima*theta) +
     $              hval(n-1)*rscale*exp(ima*theta))
     $              *(mpole(-n)*exp(ima*(-n)*theta)*(-1)**n)*zk/2
                 grad(2)=grad(2)+(-hval(n+1)/rscale*exp(ima*theta)-
     $              hval(n-1)*rscale*exp(-ima*theta))
     $              *(mpole(n)*exp(ima*n*theta))*zk*(-ima)/2
                 grad(2)=grad(2)+(-hval(n+1)/rscale*exp(-ima*theta)-
     $              hval(n-1)*rscale*exp(ima*theta))
     $              *(mpole(-n)*exp(ima*(-n)*theta)*(-1)**n)*zk*(ima)/2
              enddo
              grad(1)=grad(1)/ima4
              grad(2)=grad(2)/ima4
              write(*,*) grad(1), grad(2)
           endif
c
           grad(1)=hval(0) *mpolex(0)
           grad(2)=hval(0) *mpoley(0)
           do n=1,nterms+1
              grad(1)=grad(1)+mpolex(n)*mptemp(n)+mpolex(-n)*mptemp(-n)
              grad(2)=grad(2)+mpoley(n)*mptemp(n)+mpoley(-n)*mptemp(-n)
           enddo
           grad(1)=grad(1)/ima4
           grad(2)=grad(2)/ima4
        endif
c
c
        if( ifhess .eq. 1 ) then
           hess(1)=0
           hess(2)=0
           hess(3)=0
c
           hess(1)=hval(0) *mpolexx(0)
           hess(2)=hval(0) *mpolexy(0)
           hess(3)=hval(0) *mpoleyy(0)
           do n=1,nterms+2
              hess(1)=hess(1)+
     $                mpolexx(n)*mptemp(n)+mpolexx(-n)*mptemp(-n)
              hess(2)=hess(2)+
     $                mpolexy(n)*mptemp(n)+mpolexy(-n)*mptemp(-n)
              hess(3)=hess(3)+
     $                mpoleyy(n)*mptemp(n)+mpoleyy(-n)*mptemp(-n)
           enddo
           hess(1)=hess(1)/ima4
           hess(2)=hess(2)/ima4
           hess(3)=hess(3)/ima4
        endif
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine h2dtaeval(zk,rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n J_n(k r) exp(i n theta)
c              n=-nterms  
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad  :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,pot,grad(2),hess(3),mpole(-nterms:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        integer itt,iisc
        parameter (itt=30000)
        parameter (iisc=1000)
        complex *16, allocatable :: cw1(:)
        integer, allocatable :: iscale1(:)
        complex *16 ima,ima4,ima4inv,z,zmull,zmullinv,zfac
        complex *16 cw(0:itt)
        integer iscale(0:iisc)
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        ima4inv=ima/4
c
        lwfjs = nterms+5 + 4*nterms + 100
        ijval = 0
        ijder = ijval + lwfjs+4
        imptemp = ijder + lwfjs+4
        impolex = imptemp + 2*nterms+5
        impoley = impolex + 2*nterms+3
        impolexx = impoley + 2*nterms+3
        impolexy = impolexx + 2*nterms+5
        impoleyy = impolexy + 2*nterms+5
        itot = impoleyy + 2*nterms+5
c
        ialloc = 0
        if ((itot.gt.itt).or.(lwfjs+10 .gt. iisc)) then
           allocate(cw1(0:itot))
           allocate(iscale1(0:lwfjs+10))
           ialloc = 1
        endif
ccc        write(6,*) ' ialloc is ',ialloc
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
        call h2cart2polar(zdiff,r,theta)
c
        z=zk*r
        ifder=0
        if (ialloc.eq.0) then 
           call jfuns2d(ier,nterms+3,z,rscale,cw(ijval),ifder,cw(ijder),
     1        lwfjs,iscale,ntop)
           call mkmptemp(cw(imptemp),theta,cw(ijval),nterms)
           if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
              call mkmpole12(cw(impolex),cw(impoley),
     1             ifhess,cw(impolexx),cw(impolexy),
     2             cw(impoleyy),mpole,zk,rscale,nterms)
           endif
           call taevals(mpole,cw(impolex),cw(impoley),cw(impolexx),
     1             cw(impolexy),cw(impoleyy),rscale,nterms,
     2             cw(ijval),cw(imptemp),pot,ifgrad,grad,
     3             ifhess,hess,ima4inv)
        else
           call jfuns2d(ier,nterms+3,z,rscale,cw1(ijval),ifder,
     1        cw1(ijder),lwfjs,iscale,ntop)
           call mkmptemp(cw1(imptemp),theta,cw1(ijval),nterms)
           if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
              call mkmpole12(cw1(impolex),cw1(impoley),
     1             ifhess,cw1(impolexx),cw1(impolexy),
     2             cw1(impoleyy),mpole,zk,rscale,nterms)
           endif
           call taevals(mpole,cw1(impolex),cw1(impoley),cw1(impolexx),
     1             cw1(impolexy),cw1(impoleyy),rscale,nterms,
     2             cw1(ijval),cw1(imptemp),pot,ifgrad,grad,
     3             ifhess,hess,ima4inv)
        endif
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine h2dmpevalall
     $     (zk,rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot1,ifgrad,grad1,ifhess,hess1)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n H_n(k r) exp(i n theta)
c              n=-nterms  
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,pot,grad(2),hess(3),mpole(-nterms:nterms)
        complex *16 pot1(*),grad1(2,*),hess1(3,*)
        real *8 center(2),ztarg(2,*),zdiff(2)
c
        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
        complex *16 zmul,zinv,ztemp1,ztemp2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        ima4inv=ima/4
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
           allocate(mpolex(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-1,nterms+1
              mpolex(i)=0
              mpoley(i)=0
           enddo
           z1scale = zk/2/rscale
           z2scale = zk/2*rscale
           z3scale = zk/2/rscale * ima
           z4scale = zk/2*rscale * ima
           do i=-nterms,nterms
              if( i .le. 0 ) then
                 mpolex(i-1)=mpolex(i-1)+z1scale*mpole(i)
                 mpoley(i-1)=mpoley(i-1)+z3scale*mpole(i)
              endif
              if( i .gt. 0 ) then
                 mpolex(i-1)=mpolex(i-1)+z2scale*mpole(i)
                 mpoley(i-1)=mpoley(i-1)+z4scale*mpole(i)
              endif
              if( i .ge. 0 ) then
                 mpolex(i+1)=mpolex(i+1)-z1scale*mpole(i)
                 mpoley(i+1)=mpoley(i+1)+z3scale*mpole(i)
              endif
              if( i .lt. 0 ) then
                 mpolex(i+1)=mpolex(i+1)-z2scale*mpole(i)
                 mpoley(i+1)=mpoley(i+1)+z4scale*mpole(i)
              endif
           enddo
        endif
c
        if( ifhess .eq. 1 ) then
           allocate(mpolexx(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpolexy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoleyy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-2,nterms+2
              mpolexx(i)=0
              mpolexy(i)=0
              mpoleyy(i)=0
           enddo
           z1scale = zk/2/rscale
           z2scale = zk/2*rscale
           z3scale = zk/2/rscale * ima
           z4scale = zk/2*rscale * ima
           do i=-nterms+1,nterms+1
              if( i .le. 0 ) then
                 mpolexx(i-1)=mpolexx(i-1)+z1scale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+z3scale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+z3scale*mpoley(i)
              endif
              if( i .gt. 0 ) then
                 mpolexx(i-1)=mpolexx(i-1)+z2scale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+z4scale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+z4scale*mpoley(i)
              endif
              if( i .ge. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-z1scale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+z3scale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+z3scale*mpoley(i)
              endif
              if( i .lt. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-z2scale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+z4scale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+z4scale*mpoley(i)
              endif
           enddo
        endif
c
        allocate(mptemp(-nterms-2:nterms+2), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do k=1,ntarg
           zdiff(1)=ztarg(1,k)-center(1)
           zdiff(2)=ztarg(2,k)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call h2dall(nterms+3,z,rscale,hval,ifder,hder)
c
           mptemp(0)=hval(0)
           zmul=exp(ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul
           ztemp2=-zinv
           do j = 1,nterms+2
              mptemp( j) = ztemp1*hval(j)
              mptemp(-j) = ztemp2*hval(j)
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
           if( ifpot .eq. 1 ) then
              pot=0
              pot=hval(0)*mpole(0)
              do n=1,nterms
                 pot=pot+mpole(n)*mptemp(n)+mpole(-n)*mptemp(-n)
              enddo
              pot=pot*ima4inv
              pot1(k)=pot1(k)+pot
           endif
           if( ifgrad .eq. 1 ) then
              grad(1)=0
              grad(2)=0
              grad(1)=hval(0) *mpolex(0)
              grad(2)=hval(0) *mpoley(0)
              do n=1,nterms+1
                 grad(1)=grad(1)+
     $                   mpolex(n)*mptemp(n)+mpolex(-n)*mptemp(-n)
                 grad(2)=grad(2)+
     $                   mpoley(n)*mptemp(n)+mpoley(-n)*mptemp(-n)
              enddo
              grad(1)=grad(1)*ima4inv
              grad(2)=grad(2)*ima4inv
              grad1(1,k)=grad1(1,k)+grad(1)
              grad1(2,k)=grad1(2,k)+grad(2)
c
           endif
           if( ifhess .eq. 1 ) then
              hess(1)=0
              hess(2)=0
              hess(3)=0
              hess(1)=hval(0) *mpolexx(0)
              hess(2)=hval(0) *mpolexy(0)
              hess(3)=hval(0) *mpoleyy(0)
              do n=1,nterms+2
                 hess(1)=hess(1)+
     $                   mpolexx(n)*mptemp(n)+mpolexx(-n)*mptemp(-n)
                 hess(2)=hess(2)+
     $                   mpolexy(n)*mptemp(n)+mpolexy(-n)*mptemp(-n)
                 hess(3)=hess(3)+
     $                   mpoleyy(n)*mptemp(n)+mpoleyy(-n)*mptemp(-n)
              enddo
              hess(1)=hess(1)*ima4inv
              hess(2)=hess(2)*ima4inv
              hess(3)=hess(3)*ima4inv
              hess1(1,k)=hess1(1,k)+hess(1)
              hess1(2,k)=hess1(2,k)+hess(2)
              hess1(3,k)=hess1(3,k)+hess(3)
           endif
        enddo
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine h2dtaevalall
     $     (zk,rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot1,ifgrad,grad1,ifhess,hess1)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n J_n(k r) exp(i n theta)
c              n=-nterms  
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of targets
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier    :    error return code
c                     ier=0  successful execution
c                     ier=8  insufficient work space (not implemented yet)
c                     ier=16 ztarg too close to center
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,pot,grad(2),hess(3),mpole(-nterms:nterms)
        complex *16 pot1(*),grad1(2,*),hess1(3,*)
        real *8 center(2),ztarg(2,1),zdiff(2)
c
        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
        complex *16 zmul,zinv,ztemp1,ztemp2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        ima4inv=ima/4
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
           allocate(mpolex(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(-nterms-1:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-1,nterms+1
              mpolex(i)=0
              mpoley(i)=0
           enddo
           z1scale = zk/2*rscale
           z2scale = zk/2/rscale
           z3scale = zk/2*rscale * ima
           z4scale = zk/2/rscale * ima
           do i=-nterms,nterms
              if( i .le. 0 ) then
                 mpolex(i-1)=mpolex(i-1)+z1scale*mpole(i)
                 mpoley(i-1)=mpoley(i-1)+z3scale*mpole(i)
              endif
              if( i .gt. 0 ) then
                 mpolex(i-1)=mpolex(i-1)+z2scale*mpole(i)
                 mpoley(i-1)=mpoley(i-1)+z4scale*mpole(i)
              endif
              if( i .ge. 0 ) then
                 mpolex(i+1)=mpolex(i+1)-z1scale*mpole(i)
                 mpoley(i+1)=mpoley(i+1)+z3scale*mpole(i)
              endif
              if( i .lt. 0 ) then
                 mpolex(i+1)=mpolex(i+1)-z2scale*mpole(i)
                 mpoley(i+1)=mpoley(i+1)+z4scale*mpole(i)
              endif
           enddo
        endif
c
c
        if( ifhess .eq. 1 ) then
           allocate(mpolexx(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpolexy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoleyy(-nterms-2:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=-nterms-2,nterms+2
              mpolexx(i)=0
              mpolexy(i)=0
              mpoleyy(i)=0
           enddo
           z1scale = zk/2*rscale
           z2scale = zk/2/rscale
           z3scale = zk/2*rscale * ima
           z4scale = zk/2/rscale * ima
           do i=-nterms+1,nterms+1
              if( i .le. 0 ) then
                 mpolexx(i-1)=mpolexx(i-1)+z1scale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+z3scale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+z3scale*mpoley(i)
              endif
              if( i .gt. 0 ) then
                 mpolexx(i-1)=mpolexx(i-1)+z2scale*mpolex(i)
                 mpolexy(i-1)=mpolexy(i-1)+z4scale*mpolex(i)
                 mpoleyy(i-1)=mpoleyy(i-1)+z4scale*mpoley(i)
              endif
              if( i .ge. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-z1scale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+z3scale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+z3scale*mpoley(i)
              endif
              if( i .lt. 0 ) then
                 mpolexx(i+1)=mpolexx(i+1)-z2scale*mpolex(i)
                 mpolexy(i+1)=mpolexy(i+1)+z4scale*mpolex(i)
                 mpoleyy(i+1)=mpoleyy(i+1)+z4scale*mpoley(i)
              endif
          enddo
        endif
c
        allocate(mptemp(-nterms-2:nterms+2), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do k=1,ntarg
           zdiff(1)=ztarg(1,k)-center(1)
           zdiff(2)=ztarg(2,k)-center(2)
           call h2cart2polar(zdiff,r,theta)
c
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+3,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)
c
           mptemp(0)=jval(0)
           zmul=exp(ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul
           ztemp2=-zinv
           do j = 1,nterms+2
              mptemp( j) = ztemp1*jval(j)
              mptemp(-j) = ztemp2*jval(j)
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
           if( ifpot .eq. 1 ) then
              pot=jval(0)*mpole(0)
              do n=1,nterms
                 pot=pot+mpole(n)*mptemp(n)+mpole(-n)*mptemp(-n)
              enddo
              pot=pot*ima4inv
              pot1(k)=pot1(k)+pot
           endif
           if( ifgrad .eq. 1 ) then
              grad(1)=jval(0) *mpolex(0)
              grad(2)=jval(0) *mpoley(0)
              do n=1,nterms+1
               grad(1)=grad(1)+mpolex(n)*mptemp(n)+mpolex(-n)*mptemp(-n)
               grad(2)=grad(2)+mpoley(n)*mptemp(n)+mpoley(-n)*mptemp(-n)
              enddo
              grad(1)=grad(1)*ima4inv
              grad(2)=grad(2)*ima4inv
              grad1(1,k)=grad1(1,k)+grad(1)
              grad1(2,k)=grad1(2,k)+grad(2)
           endif
           if( ifhess .eq. 1 ) then
              hess(1)=jval(0) *mpolexx(0)
              hess(2)=jval(0) *mpolexy(0)
              hess(3)=jval(0) *mpoleyy(0)
              do n=1,nterms+2
                 hess(1)=hess(1)+
     $                   mpolexx(n)*mptemp(n)+mpolexx(-n)*mptemp(-n)
                 hess(2)=hess(2)+
     $                   mpolexy(n)*mptemp(n)+mpolexy(-n)*mptemp(-n)
                 hess(3)=hess(3)+
     $                   mpoleyy(n)*mptemp(n)+mpoleyy(-n)*mptemp(-n)
              enddo
              hess(1)=hess(1)*ima4inv
              hess(2)=hess(2)*ima4inv
              hess(3)=hess(3)*ima4inv
              hess1(1,k)=hess1(1,k)+hess(1)
              hess1(2,k)=hess1(2,k)+hess(2)
              hess1(3,k)=hess1(3,k)+hess(3)
           endif
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformmp(ier,zk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum charge_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        complex *16 zk,mpole(-nterms:nterms),charge(*)
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do n=-nterms,nterms
           mpole(n)=0
        enddo

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)
           mpole(0)=mpole(0)+charge(j)*jval(0)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul*charge(j)
           ztemp2=-zinv*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+jval(n)*ztemp1
              mpole(-n)=mpole(-n)+jval(n)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformmp_add(ier,zk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum charge_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        complex *16 zk,mpole(-nterms:nterms),charge(*)
        real *8 center(2),source(2,1),zdiff(2)
        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)

           mpole(0)=mpole(0)+charge(j)*jval(0)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul*charge(j)
           ztemp2=-zinv*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+jval(n)*ztemp1
              mpole(-n)=mpole(-n)+jval(n)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformta(ier,zk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum charge_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(-nterms:nterms),charge(*)
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do n=-nterms,nterms
        mpole(n)=0
        enddo

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call h2dall(nterms+1,z,rscale,hval,ifder,hder)
           mpole(0)=mpole(0)+charge(j)*hval(0)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul*charge(j)
           ztemp2=-zinv*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+hval(n)*ztemp1
              mpole(-n)=mpole(-n)+hval(n)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformta_add(ier,zk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum charge_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(-nterms:nterms),charge(*)
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call h2dall(nterms+1,z,rscale,hval,ifder,hder)
           mpole(0)=mpole(0)+charge(j)*hval(0)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1= zmul*charge(j)
           ztemp2=-zinv*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+hval(n)*ztemp1
              mpole(-n)=mpole(-n)+hval(n)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmpmp(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
        theta=theta-pi
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
c
        jtemp(0) = jval(0)
        zmul=exp(-ima*theta)
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           jtemp( j) = ztemp1*jval(j)
           jtemp(-j) = ztemp2*jval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*jtemp(0)
        rsj=rscale1
        do j = 1,nterms1
           jexp(0) = jexp(0)+(hexp(+j)*jtemp(-j))*rsj**2
           jexp(0) = jexp(0)+(hexp(-j)*jtemp(+j))*rsj**2
           rsj=rsj*rscale1
        enddo
c
        rsi=rscale1
        rsi7=rscale2
        rsi5=rscale1/rscale2
        do i = 1,nterms2
           jexp(i) = jexp(i) + hexp(0)*jtemp(i)*rsi5
           jexp(-i) = jexp(-i) + hexp(0)*jtemp(-i)*rsi5
           rsj=rscale1
        do j = 1,min(nterms1,i)
           jexp(i) = jexp(i)+(hexp(+j)*jtemp(i-j))*rsi5
           jexp(i) = jexp(i)+(hexp(-j)*jtemp(i+j))*rsj**2*rsi5
           jexp(-i) = jexp(-i)+(hexp(+j)*jtemp(-i-j))*rsj**2*rsi5
           jexp(-i) = jexp(-i)+(hexp(-j)*jtemp(-i+j))*rsi5
           rsj=rsj*rscale1
        enddo
        rsj=rscale1**(i+1)
        fs2=rsi5
        do j = i+1,nterms1
           jexp(i) = jexp(i)+(hexp(+j)*jtemp(i-j))*rscale1**2*fs2
           jexp(i) = jexp(i)+(hexp(-j)*jtemp(i+j))*rsj**2*rsi5
           jexp(-i) = jexp(-i)+(hexp(+j)*jtemp(-i-j))*rsj**2*rsi5
           jexp(-i) = jexp(-i)+(hexp(-j)*jtemp(-i+j))*rscale1**2*fs2
           rsj=rsj*rscale1
           fs2=fs2*rscale1**2
        enddo
        rsi=rsi*rscale1
        rsi7=rsi7*rscale2
        rsi5=rsi5*rscale1/rscale2
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dlocloc(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for local expansion
C           center1 = center of original local expansion
C           hexp    = coefficients of original local expansion
C           nterms1 = order of original local expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
c
        jtemp(0) = jval(0)
        zmul=exp(-ima*theta)
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           jtemp( j) = ztemp1*jval(j)
           jtemp(-j) = ztemp2*jval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*jtemp(0)
        do j = 1,nterms1
           jexp(0) = jexp(0)+(hexp(+j)*jtemp(-j))
           jexp(0) = jexp(0)+(hexp(-j)*jtemp(+j))
        enddo
c
        rsi=rscale1
        rsi7=rscale2
        rsi5=rscale2/rscale1
        do i = 1,nterms2
           jexp(i) = jexp(i) + hexp(0)*jtemp(i)*rsi7*rsi
           jexp(-i) = jexp(-i) + hexp(0)*jtemp(-i)*rsi7*rsi
           fs2=rsi5
           if( nterms1 .le. i ) fs2=fs2*rscale1**(2*(i-nterms1))
           do j = min(nterms1,i),1,-1
              jexp(i)=jexp(i)+(hexp(+j)*jtemp(i-j))*fs2
              jexp(i)=jexp(i)+(hexp(-j)*jtemp(i+j))*rsi7*rsi
              jexp(-i)=jexp(-i)+(hexp(+j)*jtemp(-i-j))*rsi7*rsi
              jexp(-i)=jexp(-i)+(hexp(-j)*jtemp(-i+j))*fs2
              fs2=fs2*rscale1**2
           enddo
           rsj=rscale1**(i+1)
           do j = i+1,nterms1
              jexp(i)=jexp(i)+(hexp(+j)*jtemp(i-j))*rsi5
              jexp(i)=jexp(i)+(hexp(-j)*jtemp(i+j))*rsi7*rsi
              jexp(-i)=jexp(-i)+(hexp(+j)*jtemp(-i-j))*rsi7*rsi
              jexp(-i)=jexp(-i)+(hexp(-j)*jtemp(-i+j))*rsi5
              rsj=rsj*rscale1
           enddo
           rsi=rsi*rscale1
           rsi7=rsi7*rscale2
           rsi5=rsi5*rscale2/rscale1
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmploc(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call h2dall(nterms+1,z,rscale1,hval,ifder,hder)
c        
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
c
        htemp(0) = hval(0)
        zmul=exp(-ima*theta)
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           htemp( j) = ztemp1*hval(j)
           htemp(-j) = ztemp2*hval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*htemp(0)
        do j = 1,nterms1
           jexp(0) = jexp(0)+(hexp(+j)*htemp(-j))
           jexp(0) = jexp(0)+(hexp(-j)*htemp(+j))
        enddo
c
        rsi=rscale1
        rsi2=rscale1**2
        do i = 1,nterms2
           jexp(i) = jexp(i) + hexp(0)*htemp(i)
           jexp(-i) = jexp(-i) + hexp(0)*htemp(-i)
           rsj=rscale1
           rsj2=rscale1**2
           do j = 1,min(nterms1,i)
              jexp(i) = jexp(i)+(hexp(+j)*htemp(i-j))*rsj2
              jexp(i) = jexp(i)+(hexp(-j)*htemp(i+j))
              jexp(-i) = jexp(-i)+(hexp(+j)*htemp(-i-j))
              jexp(-i) = jexp(-i)+(hexp(-j)*htemp(-i+j))*rsj2
              rsj=rsj*rscale1
              rsj2=rsj2*rscale1**2
           enddo
           do j = i+1,nterms1
              jexp(i) = jexp(i)+(hexp(+j)*htemp(i-j))*rsi2
              jexp(i) = jexp(i)+(hexp(-j)*htemp(i+j))
              jexp(-i) = jexp(-i)+(hexp(+j)*htemp(-i-j))
              jexp(-i) = jexp(-i)+(hexp(-j)*htemp(-i+j))*rsi2
           enddo
           rsi=rsi*rscale1
           rsi2=rsi2*rscale1**2
        enddo
        rsi=rscale2/rscale1
        do i = 1,nterms2
           jexp(+i) = jexp(+i)*rsi
           jexp(-i) = jexp(-i)*rsi
           rsi=rsi*rscale2/rscale1
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmploc_add(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call h2dall(nterms+1,z,rscale1,hval,ifder,hder)
c        
        htemp(0) = hval(0)
        zmul=exp(-ima*theta)
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           htemp( j) = ztemp1*hval(j)
           htemp(-j) = ztemp2*hval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*htemp(0)
        do j = 1,nterms1
           jexp(0) = jexp(0)+(hexp(+j)*htemp(-j))
           jexp(0) = jexp(0)+(hexp(-j)*htemp(+j))
        enddo
c
        rsi=rscale1
        rsi2=rscale1**2
        do i = 1,nterms2
           jexp(i) = jexp(i) + hexp(0)*htemp(i)
           jexp(-i) = jexp(-i) + hexp(0)*htemp(-i)
           rsj=rscale1
           rsj2=rscale1**2
           do j = 1,min(nterms1,i)
              jexp(i) = jexp(i)+(hexp(+j)*htemp(i-j))*rsj2
              jexp(i) = jexp(i)+(hexp(-j)*htemp(i+j))
              jexp(-i) = jexp(-i)+(hexp(+j)*htemp(-i-j))
              jexp(-i) = jexp(-i)+(hexp(-j)*htemp(-i+j))*rsj2
              rsj=rsj*rscale1
              rsj2=rsj2*rscale1**2
           enddo
           do j = i+1,nterms1
              jexp(i) = jexp(i)+(hexp(+j)*htemp(i-j))*rsi2
              jexp(i) = jexp(i)+(hexp(-j)*htemp(i+j))
              jexp(-i) = jexp(-i)+(hexp(+j)*htemp(-i-j))
              jexp(-i) = jexp(-i)+(hexp(-j)*htemp(-i+j))*rsi2
           enddo
           rsi=rsi*rscale1
           rsi2=rsi2*rscale1**2
        enddo
c
        rsi=rscale2/rscale1
        do i = 1,nterms2
           jexp(+i) = jexp(+i)*rsi
           jexp(-i) = jexp(-i)*rsi
           rsi=rsi*rscale2/rscale1
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmploc_add_trunc(zk,
     $     rscale1,center1,hexp,nterms1,nterms1_trunc,
     $     rscale2,center2,jexp,nterms2,nterms2_trunc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion
C           with truncation, and add to jexp.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           nterms1_trunc = number of terms used in the multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C           nterms2_trunc = number of terms used in the shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
c
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2)
        complex *16, allocatable :: hexp1(:), jexp1(:)

        allocate(hexp1(-nterms1_trunc:nterms1_trunc), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(jexp1(-nterms2_trunc:nterms2_trunc), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do i = -nterms1_trunc,nterms1_trunc
           hexp1(i) = hexp(i)
        enddo
c
        call h2dmploc(zk,
     $     rscale1,center1,hexp1,nterms1_trunc,
     $     rscale2,center2,jexp1,nterms2_trunc)
c
        do i = -nterms2_trunc,nterms2_trunc
           jexp(i) = jexp(i) + jexp1(i)
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmpmp_old(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c
c
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
c
        jtemp(0) = jval(0)
        zmul=exp(-ima*theta)*rscale1
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           jtemp( j) = ztemp1*jval(j)
           jtemp(-j) = ztemp2*jval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        do i = -nterms2,nterms2
           jexp(i) = jexp(i) + hexp(0)*jtemp(i)
           rsj=rscale1
           do j = 1,nterms1
              jexp(i) = jexp(i)+
     $          (hexp(+j)*jtemp(i-j)+hexp(-j)*jtemp(i+j))*rsj
              rsj=rsj*rscale1
           enddo
        enddo
        rsi=rscale2
        do i = 1,nterms2
           jexp(+i) = jexp(+i)/rsi
           jexp(-i) = jexp(-i)/rsi
           rsi=rsi*rscale2
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dlocloc_old(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
        jtemp(0) = jval(0)
        zmul=exp(-ima*theta)*rscale1
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           jtemp( j) = ztemp1*jval(j)
           jtemp(-j) = ztemp2*jval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        rscale1inv=1/rscale1
        do i = -nterms2,nterms2
           jexp(i) = jexp(i) + hexp(0)*jtemp(i)
           rsj=rscale1inv
           do j = 1,nterms1
              jexp(i) = jexp(i)+
     $           (hexp(+j)*jtemp(i-j)+hexp(-j)*jtemp(i+j))*rsj
              rsj=rsj*rscale1inv
           enddo
        enddo
        rsi=rscale2
        do i = 1,nterms2
           jexp(+i) = jexp(+i)*rsi
           jexp(-i) = jexp(-i)*rsi
           rsi=rsi*rscale2
        enddo
        return
        end
c
c
c
c
c
        subroutine h2dmploc_old(zk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
        complex *16 zk,hexp(-nterms1:nterms1),jexp(-nterms2:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(-nterms-5:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call h2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        z=zk*r
        ifder=0
        call h2dall(nterms+1,z,rscale1,hval,ifder,hder)
c
        do i = -nterms2,nterms2
           jexp(i) = 0
        enddo
c
        htemp(0) = hval(0)
        zmul=exp(-ima*theta)/rscale1
        zinv=conjg(zmul)
        ztemp1= zmul
        ztemp2=-zinv
        do j = 1,nterms
           htemp( j) = ztemp1*hval(j)
           htemp(-j) = ztemp2*hval(j)
           ztemp1= ztemp1*zmul
           ztemp2=-ztemp2*zinv
        enddo
c
        do i = -nterms2,nterms2
           jexp(i) = jexp(i) + hexp(0)*htemp(i)
           rsj=rscale1
           do j = 1,nterms1
              jexp(i) = jexp(i)+
     $           (hexp(+j)*htemp(i-j)+hexp(-j)*htemp(i+j))*rsj
              rsj=rsj*rscale1
           enddo
        enddo
        rsi=rscale2
        do i = 1,nterms2
           jexp(+i) = jexp(+i)*rsi
           jexp(-i) = jexp(-i)*rsi
           rsi=rsi*rscale2
        enddo
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine hpotgrad2dall(ifgrad,ifhess,sources,charge,ns,
     1                   target,wavek,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad         : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2)
      complex *16 wavek,pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call hpotgrad2d(ifgrad,ifhess,sources(1,i),charge(i),target,
     1        wavek,potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine hpotgrad2d(ifgrad,ifhess,source,charge,target,wavek,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a charge at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 wavek,pot,grad(2),hess(3)
      complex *16 charge
      complex *16 z, h0, h1, h2, cd, h2z, zk, ima, ima4, ima4inv
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4=-4*ima
      ima4inv=ima/4
c
      call prin2(' source is *',source,2)
      call prin2(' target is *',target,2)
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
      call prin2(' z is *',z,2)
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
      call prin2(' h0 is *',h0,2)
      pot = h0*charge*ima4inv
      call prin2(' pot is *',pot,2)
c
      if (ifgrad.eq.1) then
         cd = -h1*(wavek*charge*ima4inv/r)
         grad(1) = cd*xdiff
         grad(2) = cd*ydiff
c         ctheta=xdiff/r
c         stheta=ydiff/r
c         cd = -h1*(wavek*charge*ima4inv)
c         grad(1) = cd*ctheta
c         grad(2) = cd*stheta
      endif
c
      if (ifhess.eq.1) then
         cd = (wavek*charge*ima4inv/r)/rr
         h2z=(-z*h0+2*h1)
         hess(1) = cd*(h2z*xdiff*xdiff-rr*h1)
         hess(2) = cd*(h2z*xdiff*ydiff      )
         hess(3) = cd*(h2z*ydiff*ydiff-rr*h1)
c         ctheta=xdiff/r
c         stheta=ydiff/r
c         cd = (wavek*charge*ima4inv/r)
c         h2z=(-z*h0+2*h1)
c         hess(1) = cd*(h2z*ctheta*ctheta-h1)
c         hess(2) = cd*(h2z*ctheta*stheta   )
c         hess(3) = cd*(h2z*stheta*stheta-h1)
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h2cart2polar(zat,r,theta)
c**********************************************************************
c
c     Convert form Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r      :  |zat|
c       theta  : angle of (zat(1),zat(2)) subtended with 
c                 respect to x-axis
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 zat(2)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2)
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      theta = 0
      else
      theta = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c**********************************************************************
      subroutine h2dall(nterms,z,rscale,hvec,ifder,hder)
c**********************************************************************
c
c     This subroutine computes scaled versions of the Hankel 
c     functions H_n of orders 0 to nterms.
c
c       	hvec(n)= H_n(z)*rscale^(n)
c
c     The parameter SCALE is useful when |z| < 1, in which case
c     it damps out the rapid growth of H_n as n increases. In such 
c     cases, we recommend setting 
c                                 
c               rscale approx |z|
c
c     or something close. If |z| > 1, set scale = 1.
c
c     If the flag IFDER is set to one, it also computes the 
c     derivatives of h_n.
c
c		hder(n)= H_n'(z)*rscale^(n)
c
c     NOTE: If |z| < 1.0d-200, the subroutine returns zero.
c     
c-----------------------------------------------------------------------
c     INPUT:
c
c     nterms  : highest order of the Hankel functions to be computed.
c     z       : argument of the Hankel functions.
c     rscale   : scaling parameter discussed above
c     ifder   : flag indcating whether derivatives should be computed.
c		ifder = 1   ==> compute 
c		ifder = 0   ==> do not compute
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     hvec    : the vector of Hankel functions 
c     hder    : the derivatives of the Hankel functions 
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 hvec(0:*),hder(0:*)
      complex *16 zk2,z,zinv,ztmp,fhextra,h0,h1
c
      data thresh/1.0d-200/,done/1.0d0/
c
c     If |z| < thresh, return zeros.
c
      if (abs(z).lt.thresh) then
         do i=0,nterms
            hvec(i)=0
            hder(i)=0
         enddo
         return
      endif
c
c     Otherwise, get H_0 and H_1 via hank103 and the rest via
c     recursion.
c
      ifexpon=1
      call hank103(z,h0,h1,ifexpon)
      hvec(0)=h0
      hvec(1)=h1*rscale
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c       H_{n-1}(z) + H_{n+1}(z) = 2*n/z * H_n(z)
c
c     With scaling:
c
c       H_{n-1}(z) *rscale + H_{n+1}(z)/rscale = 2*n/z * H_n(z)
c       H_{n+1}(z) = rscale*(2*n/z*H_n(z) - H_{n-1}(z)*rscale)
c
      scal2=rscale*rscale
      zinv=rscale/z
      do i=1,nterms-1
         dtmp=2*i
         ztmp=zinv*dtmp
         hvec(i+1)=ztmp*hvec(i)-scal2*hvec(i-1)
      enddo
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c     H_{n}'(z)= H_{n-1}(z) - (n)/z * H_n(z)
c
c     With scaling:
c
c     hder(n)=scale* hvec(n-1) - n/z * hvec(n)
c
      if (ifder.eq.1) then
c
         hder(0)=-hvec(1)/rscale
         zinv=1.0d0/z
         do i=1,nterms
            dtmp=(i)
            ztmp=zinv*dtmp
            hder(i)=rscale*hvec(i-1)-ztmp*hvec(i)
         enddo
c
      endif
c
      return
      end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformmp_dp(ier,zk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        complex *16 zk,mpole(-nterms:nterms),dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do n=-nterms,nterms
           mpole(n)=0
        enddo

        do j=1,ns
c
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)

           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1=-zmul*dipstr(j)*zk/2
           ztemp2=+zinv*dipstr(j)*zk/2
           ztemp3= zinv/rscale
           ztemp4= zmul*rscale
           ztemp5= zmul/rscale
           ztemp6= zinv*rscale
           ztemp3 = ztemp3*(-dipvec(1,j)+ima*dipvec(2,j))
           ztemp4 = ztemp4*(+dipvec(1,j)+ima*dipvec(2,j))
           ztemp5 = ztemp5*(-dipvec(1,j)-ima*dipvec(2,j))
           ztemp6 = ztemp6*(+dipvec(1,j)-ima*dipvec(2,j))
           mpole(0)=mpole(0)-dipstr(j)*jval(1)*zk/2*rscale*
     $       (
     $       (zmul+zinv)*dipvec(1,j)+
     $       (zmul-zinv)*ima*dipvec(2,j)
     $       )
c
           do n=1,nterms
              mpole(n)=mpole(n)+
     $                  (jval(n-1)*ztemp3+jval(n+1)*ztemp4)*ztemp1
              mpole(-n)=mpole(-n)+
     $                  (jval(n-1)*ztemp5+jval(n+1)*ztemp6)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformmp_dp_add(ier,zk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        complex *16 zk,mpole(-nterms:nterms),dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1          lwfjs,iscale,ntop)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1=-zmul*dipstr(j)*zk/2
           ztemp2=+zinv*dipstr(j)*zk/2
           ztemp3= zinv/rscale
           ztemp4= zmul*rscale
           ztemp5= zmul/rscale
           ztemp6= zinv*rscale
           ztemp3 = ztemp3*(-dipvec(1,j)+ima*dipvec(2,j))
           ztemp4 = ztemp4*(+dipvec(1,j)+ima*dipvec(2,j))
           ztemp5 = ztemp5*(-dipvec(1,j)-ima*dipvec(2,j))
           ztemp6 = ztemp6*(+dipvec(1,j)-ima*dipvec(2,j))
           mpole(0)=mpole(0)-dipstr(j)*jval(1)*zk/2*rscale*
     $        (
     $        (zmul+zinv)*dipvec(1,j)+
     $        (zmul-zinv)*ima*dipvec(2,j)
     $        )
           do n=1,nterms
              mpole(n)=mpole(n)+
     $                 (jval(n-1)*ztemp3+jval(n+1)*ztemp4)*ztemp1
              mpole(-n)=mpole(-n)+
     $                  (jval(n-1)*ztemp5+jval(n+1)*ztemp6)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformta_dp(ier,zk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(-nterms:nterms),dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do n=-nterms,nterms
           mpole(n)=0
        enddo

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call h2dall(nterms+2,z,rscale,hval,ifder,hder)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1=-zmul*dipstr(j)*zk/2
           ztemp2=+zinv*dipstr(j)*zk/2
           ztemp3= zinv*rscale
           ztemp4= zmul/rscale
           ztemp5= zmul*rscale
           ztemp6= zinv/rscale
           ztemp3 = ztemp3*(-dipvec(1,j)+ima*dipvec(2,j))
           ztemp4 = ztemp4*(+dipvec(1,j)+ima*dipvec(2,j))
           ztemp5 = ztemp5*(-dipvec(1,j)-ima*dipvec(2,j))
           ztemp6 = ztemp6*(+dipvec(1,j)-ima*dipvec(2,j))
c
           mpole(0)=mpole(0)-dipstr(j)*hval(1)*zk/2/rscale*
     $        (
     $        (zmul+zinv)*dipvec(1,j)+
     $        (zmul-zinv)*ima*dipvec(2,j)
     $        )
           do n=1,nterms
              mpole( n)=mpole( n)+
     $                  (hval(n-1)*ztemp3+hval(n+1)*ztemp4)*ztemp1
              mpole(-n)=mpole(-n)+
     $                  (hval(n-1)*ztemp5+hval(n+1)*ztemp6)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine h2dformta_dp_add(ier,zk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(-nterms:nterms),dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call h2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call h2dall(nterms+2,z,rscale,hval,ifder,hder)
           zmul=exp(-ima*theta)
           zinv=conjg(zmul)
           ztemp1=-zmul*dipstr(j)*zk/2
           ztemp2=+zinv*dipstr(j)*zk/2
           ztemp3= zinv*rscale
           ztemp4= zmul/rscale
           ztemp5= zmul*rscale
           ztemp6= zinv/rscale
           ztemp3 = ztemp3*(-dipvec(1,j)+ima*dipvec(2,j))
           ztemp4 = ztemp4*(+dipvec(1,j)+ima*dipvec(2,j))
           ztemp5 = ztemp5*(-dipvec(1,j)-ima*dipvec(2,j))
           ztemp6 = ztemp6*(+dipvec(1,j)-ima*dipvec(2,j))
           mpole(0)=mpole(0)-dipstr(j)*hval(1)*zk/2/rscale*
     $        (
     $        (zmul+zinv)*dipvec(1,j)+
     $        (zmul-zinv)*ima*dipvec(2,j)
     $        )
           do n=1,nterms
              mpole(n)=mpole(n)+
     $                  (hval(n-1)*ztemp3+hval(n+1)*ztemp4)*ztemp1
              mpole(-n)=mpole(-n)+
     $                  (hval(n-1)*ztemp5+hval(n+1)*ztemp6)*ztemp2
              ztemp1= ztemp1*zmul
              ztemp2=-ztemp2*zinv
           enddo
        enddo
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine hpotgrad2dall_dp(ifgrad,ifhess,sources,dipstr,dipvec,
     1                   ns,target,wavek,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot  = H_0(k*r)*(eye/4)
c		grad  = gradient(pot)
c		hess = Hessian
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad         : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     dipstr        : dipole strengths
c     dipvec        : dipole vectors
c     ns            : number of sources
c     target        : location of the target
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2),dipvec(2,ns)
      complex *16 wavek,pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call hpotgrad2d_dp(ifgrad,ifhess,sources(1,i),dipstr(i),
     1        dipvec(1,i),target,wavek,potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine hpotgrad2d_dp(ifgrad,ifhess,source,dipstr,dipvec,
     $     target,wavek,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     dipstr    : dipole strength
c     dipvec    : dipole orientation
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 wavek,pot,grad(2),hess(3)
      complex *16 dipstr
      real *8 dipvec(2)
      complex *16 z, h0, h1, h2, h3, h4, cd, h2z, zk, ima, ima4, ima4inv
      complex *16 hx,hy
      complex *16 hxx,hxy,hyy
      complex *16 hxxx,hxxy,hxyy,hyyy
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4=-4*ima
      ima4inv=ima/4
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
      ctheta=xdiff/r
      stheta=ydiff/r

      h2 = 2*h1/z-h0
      h3 = 4*h2/z-h1
c
      cd = h1/r*wavek*dipstr*ima4inv
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)
      pot = cd*dotprod
c
      if (ifgrad.eq.1) then
         cd = -wavek**2*dipstr*ima4inv
         hxx = cd*(h2*(ctheta*ctheta-0.5d0)-h0/2)
         hxy = cd*(h2*ctheta*stheta             )
         hyy = cd*(h2*(stheta*stheta-0.5d0)-h0/2)
         grad(1) = hxx*dipvec(1)+hxy*dipvec(2)
         grad(2) = hxy*dipvec(1)+hyy*dipvec(2)
      endif
c
      if (ifhess.eq.1) then
c         hess(1) = 0
c         hess(2) = 0
c         hess(3) = 0
         cd = -wavek**3*dipstr*ima4inv
         hxxx = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $      )*ctheta
         hxxy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $      )*stheta
         hxyy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $      )*ctheta
         hyyy = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $      )*stheta
         hess(1) = hxxx*dipvec(1)+hxxy*dipvec(2)
         hess(2) = hxxy*dipvec(1)+hxyy*dipvec(2)
         hess(3) = hxyy*dipvec(1)+hyyy*dipvec(2)
      endif
c
      return
      end
c
c
c
c**********************************************************************
c
c     ... direct calculation for a collection of charges and dipoles
c
c**********************************************************************
      subroutine hpotgrad2dall_sdp(wavek,sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot  = H_0(k*r)*(eye/4)
c		grad  = gradient(pot)
c		hess = Hessian
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing 
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     dipstr        : dipole strengths
c     dipvec        : dipole vectors
c     ns            : number of sources
c     target        : location of the target
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 sources(2,ns),target(2),dipvec(2,ns)
      complex *16 wavek,pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      complex *16 dipstr(ns),charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call hpotgrad2d_sdp(wavek,sources(1,i),
     1   ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     $   target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
         if (ifpot.eq.1) pot = pot + potloc
         if (ifgrad.eq.1) then
            grad(1) = grad(1) + gradloc(1)
            grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
            hess(1) = hess(1) + hessloc(1)
            hess(2) = hess(2) + hessloc(2)
            hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine hpotgrad2d_sdp(wavek,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing 
c	                	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 wavek,pot,grad(2),hess(3)
      complex *16 charge,dipstr
      real *8 dipvec(2)
      complex *16 z, h0, h1, h2, h3, h4, cd, h2z, zk, ima, ima4, ima4inv
      complex *16 hx,hy
      complex *16 hxx,hxy,hyy
      complex *16 hxxx,hxxy,hxyy,hyyy
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4=-4*ima
      ima4inv=ima/4
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
c
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif
c
c
      if( ifcharge .eq. 1 ) then
c
      if( ifpot.eq.1) pot = h0*charge*ima4inv
c
      if (ifgrad.eq.1) then
         cd = -h1*(wavek*charge*ima4inv/r)
         grad(1) = cd*xdiff
         grad(2) = cd*ydiff
      endif
c
      if (ifhess.eq.1) then
         cd = (wavek*charge*ima4inv/r)/rr
         h2z=(-z*h0+2*h1)
         hess(1) = cd*(h2z*xdiff*xdiff-rr*h1)
         hess(2) = cd*(h2z*xdiff*ydiff      )
         hess(3) = cd*(h2z*ydiff*ydiff-rr*h1)
      endif
c
      endif
c
c
      if( ifdipole .eq. 1 ) then
      
      ctheta=xdiff/r
      stheta=ydiff/r

      h2 = 2*h1/z-h0
      h3 = 4*h2/z-h1
c
      cd = h1/r*wavek*dipstr*ima4inv
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)

      if(ifpot.eq.1) pot = pot+cd*dotprod
c
      if (ifgrad.eq.1) then
         cd = -wavek**2*dipstr*ima4inv
         hxx = cd*(h2*(ctheta*ctheta-0.5d0)-h0/2)
         hxy = cd*(h2*ctheta*stheta             )
         hyy = cd*(h2*(stheta*stheta-0.5d0)-h0/2)
         grad(1) = grad(1) + hxx*dipvec(1)+hxy*dipvec(2)
         grad(2) = grad(2) + hxy*dipvec(1)+hyy*dipvec(2)
      endif
c
      if (ifhess.eq.1) then
         cd = -wavek**3*dipstr*ima4inv
         hxxx = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $      )*ctheta
         hxxy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $      )*stheta
         hxyy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $      )*ctheta
         hyyy = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $      )*stheta
         hess(1) = hess(1) + hxxx*dipvec(1)+hxxy*dipvec(2)
         hess(2) = hess(2) + hxxy*dipvec(1)+hxyy*dipvec(2)
         hess(3) = hess(3) + hxyy*dipvec(1)+hyyy*dipvec(2)
      endif
c
      endif
c
      return
      end
c
c
c 
c**********************************************************************
      subroutine mkmptemp(mptemp,theta,jval,nterms)
c**********************************************************************
c
c     This subroutine computes the complex terms in a partial
c     wave expansion:
c                       J_n(k r) exp(i n theta)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     theta   :    angular argument
c     jval    :    array of J_n values
c     nterms  :    order of expansion is [-nterms-2:nterms+2]
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mptemp  :    J_n(kr)*exp(i n theta)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 mptemp(-nterms-2:nterms+2)
      complex *16 jval(0:nterms+2)
      real *8 theta
      complex *16 ima,zfac,zfac2,zmull,zmullinv
      data ima/(0.0d0,1.0d0)/
c
      mptemp(0)=jval(0)
      zfac = exp(ima*theta)
      zfac2 = -dconjg(zfac)
      zmull = zfac
      zmullinv = zfac2
      do n=1,nterms+2
         mptemp(n)=jval(n)*zmull
         mptemp(-n)=jval(n)*zmullinv
         zmull = zmull*zfac
         zmullinv = zmullinv*zfac2
      enddo
      return
      end
c
c
c
       subroutine mkmpole12(mpolex,mpoley,ifhess,mpolexx,
     1               mpolexy,mpoleyy,mpole,zk,rscale,nterms)
c
c     This subroutine converts multipole expansion for
c     potential into expansions for various derivatives.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     ifhess  :    flag determines whether 2nd derivs are desired.
c     mpole   :    multipole expansion 
c     zk      :    Helmholtz parameter
c     rscale  :    scaling parameter for expansion
c     nterms  :    order of expansion
c                  for mpole [-nterms:nterms] 
c                  for mpolex,mpoley [-nterms-1:nterms+1] 
c                  for 2nd derivs [-nterms-2:nterms+2] 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpolex,mpoley   : expansions for x,y derivatives
c     mpolexx,mpolexy,
c     mpoleyy         : expansions for 2nd derivatives if requested
c
c
c
       implicit real *8 (a-h,o-z)
       complex *16 mpolex(-nterms-1:nterms+1)
       complex *16 mpoley(-nterms-1:nterms+1)
       complex *16 mpolexx(-nterms-2:nterms+2)
       complex *16 mpolexy(-nterms-2:nterms+2)
       complex *16 mpoleyy(-nterms-2:nterms+2)
       complex *16 mpole(-nterms:nterms)
       complex *16 ima,zk
       data ima/(0.0d0,1.0d0)/
c
       do i=-nterms-1,nterms+1
          mpolex(i)=0
          mpoley(i)=0
       enddo
       do i=-nterms,nterms
          if (i.le.0) then
             mpolex(i-1)=mpolex(i-1)+zk/2*rscale*mpole(i)
             mpoley(i-1)=mpoley(i-1)+zk/2*(ima)*rscale*mpole(i)
          else
             mpolex(i-1)=mpolex(i-1)+zk/2/rscale*mpole(i)
             mpoley(i-1)=mpoley(i-1)+zk/2*(ima)/rscale*mpole(i)
          endif
          if( i .ge. 0 ) then
             mpolex(i+1)=mpolex(i+1)-zk/2*rscale*mpole(i)
             mpoley(i+1)=mpoley(i+1)+zk/2*(ima)*rscale*mpole(i)
          else
             mpolex(i+1)=mpolex(i+1)-zk/2/rscale*mpole(i)
             mpoley(i+1)=mpoley(i+1)+zk/2*(ima)/rscale*mpole(i)
          endif
       enddo
c
       if (ifhess.eq.1) then
          do i=-nterms-2,nterms+2
             mpolexx(i)=0
             mpolexy(i)=0
             mpoleyy(i)=0
          enddo
          do i=-nterms+1,nterms+1
             if (i.le.0) then
                mpolexx(i-1)=mpolexx(i-1)+zk/2*rscale*mpolex(i)
                mpolexy(i-1)=mpolexy(i-1)+zk/2*(ima)*rscale*mpolex(i)
                mpoleyy(i-1)=mpoleyy(i-1)+zk/2*(ima)*rscale*mpoley(i)
             else
                mpolexx(i-1)=mpolexx(i-1)+zk/2/rscale*mpolex(i)
                mpolexy(i-1)=mpolexy(i-1)+zk/2*(ima)/rscale*mpolex(i)
                mpoleyy(i-1)=mpoleyy(i-1)+zk/2*(ima)/rscale*mpoley(i)
             endif
             if (i.ge.0) then
                mpolexx(i+1)=mpolexx(i+1)-zk/2*rscale*mpolex(i)
                mpolexy(i+1)=mpolexy(i+1)+zk/2*(ima)*rscale*mpolex(i)
                mpoleyy(i+1)=mpoleyy(i+1)+zk/2*(ima)*rscale*mpoley(i)
             else
                mpolexx(i+1)=mpolexx(i+1)-zk/2/rscale*mpolex(i)
                mpolexy(i+1)=mpolexy(i+1)+zk/2*(ima)/rscale*mpolex(i)
                mpoleyy(i+1)=mpoleyy(i+1)+zk/2*(ima)/rscale*mpoley(i)
             endif
          enddo
       endif
       return
       end

      subroutine taevals(mpole,mpolex,mpoley,mpolexx,mpolexy,mpoleyy,
     1                   rscale,nterms,jval,mptemp,pot,
     2                   ifgrad,grad,ifhess,hess,ima4inv)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     mpole           : multipole expansion 
c     mpolex,mpoley   : expansions for x,y derivatives if requested
c     mpolexx,mpolexy,
c     mpoleyy         : expansions for 2nd derivatives if requested
c     rscale          :    scaling parameter for expansion
c     nterms          :    order of expansion
c                          for mpole [-nterms:nterms] 
c                          for mpolex,mpoley [-nterms-1:nterms+1] 
c                          for 2nd derivs [-nterms-2:nterms+2] 
c     jval            : Bessel expansion values
c     mptemp          : J_n exp(in theta) values for target from
c                       prior call to mkmptemp
c     ifgrad          : flag to request gradient computation
c     ifhess          : flag to request Hessian computation
c     ima4inv         : scaling factor (1/4i)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot             : potential
c     grad            : gradient (if requested)  (pot_x,pot_y)
c     hess            : Hessian (if requested) (pot_xx,pot_xy,pot_yy)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 mpole(-nterms:nterms)
      complex *16 mptemp(-nterms-2:nterms+2)
      complex *16 mpolex(-nterms-1:nterms+1)
      complex *16 mpoley(-nterms-1:nterms+1)
      complex *16 mpolexx(-nterms-2:nterms+2)
      complex *16 mpolexy(-nterms-2:nterms+2)
      complex *16 mpoleyy(-nterms-2:nterms+2)
      complex *16 jval(0:nterms+3)
      complex *16 pot,grad(2),hess(3)
      complex *16 ima4inv
c
      pot=jval(0)*mpole(0)
      do n=1,nterms
         pot=pot+mpole(n)*mptemp(n)+mpole(-n)*mptemp(-n)
      enddo
      pot=pot*ima4inv
c
      if (ifgrad .eq. 1) then
         grad(1)=0
         grad(2)=0
         grad(1)=jval(0) *mpolex(0)
         grad(2)=jval(0) *mpoley(0)
         do n=1,nterms+1
            grad(1)=grad(1)+mpolex(n)*mptemp(n)+mpolex(-n)*mptemp(-n)
            grad(2)=grad(2)+mpoley(n)*mptemp(n)+mpoley(-n)*mptemp(-n)
         enddo
         grad(1)=grad(1)*ima4inv
         grad(2)=grad(2)*ima4inv
      endif
c
      if (ifhess .eq. 1) then
         hess(1)=0
         hess(2)=0
         hess(3)=0
         hess(1)=jval(0) *mpolexx(0)
         hess(2)=jval(0) *mpolexy(0)
         hess(3)=jval(0) *mpoleyy(0)
         do n=1,nterms+2
            hess(1)=hess(1)+mpolexx(n)*mptemp(n)+mpolexx(-n)*mptemp(-n)
            hess(2)=hess(2)+mpolexy(n)*mptemp(n)+mpolexy(-n)*mptemp(-n)
            hess(3)=hess(3)+mpoleyy(n)*mptemp(n)+mpoleyy(-n)*mptemp(-n)
         enddo
         hess(1)=hess(1)*ima4inv
         hess(2)=hess(2)*ima4inv
         hess(3)=hess(3)*ima4inv
      endif
      return
      end
c
c
c
