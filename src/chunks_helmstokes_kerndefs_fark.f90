
!***************************************************************
!
! Kernels and Green's functions for implementing the
! integral equations for the helmholtz stokes 
! The integral equations are the extension of
! one's presented in Peter Farkas' thesis and the Jiang Quaife Kropinski
! paper
!
!***************************************************************

subroutine zfark_kernel(par0, src, targ, src_normal, &
     targ_normal, q1, q2, fgreen, pars1, pars2, val)
  implicit real *8 (a-h,o-z)
  real *8 src(2), targ(2), src_normal(2)
  real *8 targ_normal(2), pars1(*), pars2(*)
  complex *16 par0, q1, q2, g, dnyg, dnxg, dnxdnyg
  complex *16 val

  complex *16 val2, grad2(10), hess2(2,2), vec(10)

  call fgreen(par0, src, targ, src_normal, pars2, val2, &
       grad2, hess2)

  val = val2

  return
end subroutine zfark_kernel
!--------------------------------
!
subroutine zfark_kernelp(par0, src, targ, src_normal, &
     targ_normal, q1, q2, fgreen, pars1, pars2, val)
  implicit real *8 (a-h,o-z)
  real *8 src(2), targ(2), src_normal(2)
  real *8 targ_normal(2), pars1(*), pars2(*)
  complex *16 par0, q1, q2, g, dnyg, dnxg, dnxdnyg
  complex *16 val

  complex *16 val2, grad2(2), hess2(2,2)

  call fgreen(par0, src, targ, src_normal, pars2, val2, &
       grad2, hess2)

  val = targ_normal(1)*grad2(1) + targ_normal(2)*grad2(2)

  return
end subroutine zfark_kernelp

!----------------------------------
subroutine zhfark_ck2(zk,src,targ,src_normal, &
     pars2,val,grad,hess)
  implicit real *8 (a-h,o-z)
  ! global
  complex *16, intent(in) :: zk
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: pars2(*)
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: octvec(4),quadvec(3), dipvec(2)
  complex *16 :: octstr(1),quadstr(1), dipstr(1)
  complex *16 :: dpot,dgrad(2),dhess(4) 
  real *8 :: charge(1)
  real *8 :: tautemp(2), rnormtemp(2)
  integer :: ifcharge, ifdipole, ifquad, ns

  tautemp(1) = -src_normal(2)
  tautemp(2) = src_normal(1)

  rnormtemp(1) = src_normal(1)
  rnormtemp(2) = src_normal(2)

  ifcharge = 0
  ifdipole = 0
  ifquad = 1
  ifoct = 0

  ! quadrupole corresponding to G_tt - G_nn

  quadstr(1) = 1.0d0
  quadvec(1) = tautemp(1)**2 - rnormtemp(1)**2
  quadvec(2) = 2*tautemp(1)*tautemp(2) &
       - 2*rnormtemp(1)*rnormtemp(2)
  quadvec(3) = tautemp(2)**2 - rnormtemp(2)**2

  ns = 1

  ifpot = 1
  ifgrad = 1
  ifhess = 0

!  call prin2('zk=*',zk,2)
!  call prin2('src=*',src,2)
!  call prinf('ifcharge=*',ifcharge,1)
!  call prinf('ifquad=*',ifquad,1)
!  call prin2('quadstr=*',quadstr,2)
!  call prin2('quadvec=*',quadvec,3)
!  call prinf('ifdipole=*',ifdipole,1)
!  call prinf('ifoct=*',ifoct,1)
!  call prinf('ifpot=*',ifpot,1)



  call hbhpotgrad2dall_cdqo(zk,src,ns,ifcharge, &
       charge,ifdipole,dipstr,dipvec,ifquad,quadstr, &
       quadvec,ifoct,octstr,octvec,targ,ifpot,dpot,ifgrad, &
       dgrad,ifhess,dhess)

!  call prin2('dpot=*',dpot,2)



  val = dpot
  grad(1) = dgrad(1)
  grad(2) = dgrad(2)
  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  return
end subroutine zhfark_ck2

subroutine zhfark_ck1(zk,src,targ,src_normal, &
     pars2,val,grad,hess)
  implicit real *8 (a-h,o-z)
  ! global
  complex *16, intent(in) :: zk
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: pars2(*)
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: octvec(4),quadvec(3), dipvec(2)
  complex *16 :: octstr(1),quadstr(1), dipstr(1)
  complex *16 :: dpot,dgrad(2),dhess(4)
  real *8 :: charge(1)
  real *8 :: tautemp(2), rnormtemp(2)
  integer :: ifcharge, ifdipole, ifquad, ns

  tautemp(1) = -src_normal(2)
  tautemp(2) = src_normal(1)

  rnormtemp(1) = src_normal(1)
  rnormtemp(2) = src_normal(2)

  ifcharge = 0
  ifdipole = 1
  ifquad = 0
  ifoct = 1

  ! dipole  corresponding to \lambda^2 G_n
  
  dipstr(1) = abs(zk)**2
  dipvec(1) = rnormtemp(1)
  dipvec(2) = rnormtemp(2)


  ! octopole  corresponding to 3G_ttn + G_nnn 

  octstr(1) = 1.0d0
  octvec(1) = 3*rnormtemp(1)*tautemp(1)**2 + rnormtemp(1)**3
  octvec(2) = 6*rnormtemp(1)*tautemp(1)*tautemp(2) + &
       3*rnormtemp(2)*tautemp(1)**2 + 3*rnormtemp(1)**2*rnormtemp(2) 
  octvec(3) = 6*rnormtemp(2)*tautemp(1)*tautemp(2) + &
       3*rnormtemp(1)*tautemp(2)**2 + 3*rnormtemp(1)*rnormtemp(2)**2 
  octvec(4) = 3*rnormtemp(2)*tautemp(2)**2 + rnormtemp(2)**3 
   

  ns = 1

  ifpot = 1
  ifgrad = 1
  ifhess = 0

!  call prinf('ifdipole=*',ifdipole,1)
!  call prinf('ifquad=*',ifquad,1)
!  call prinf('ifoct=*',ifoct,1)
!  call prin2('octvec=*',octvec,4)

  call hbhpotgrad2dall_cdqo(zk,src,ns,ifcharge, &
       charge,ifdipole,dipstr,dipvec,ifquad,quadstr, &
       quadvec,ifoct,octstr,octvec,targ,ifpot,dpot,ifgrad, &
       dgrad,ifhess,dhess)

!  call prin2('dpot=*',dpot,2)

  val = dpot
  grad(1) = dgrad(1)
  grad(2) = dgrad(2)
  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  return
end subroutine zhfark_ck1

