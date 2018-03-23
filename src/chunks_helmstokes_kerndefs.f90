
subroutine zchker_helmstokes_s1(zk,src,targ,src_normal, &
     pars2,val,grad,hess)
  implicit real *8 (a-h,o-z)
  ! global
  complex *16, intent(in) :: zk
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: pars2(*)
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: quadvec(3), dipvec(2)
  complex *16 :: quadstr(1), dipstr(1)
  real *8 :: charge(1)
  real *8 :: tautemp(2), rnormtemp(2)
  integer :: ifcharge, ifdipole, ifquad, ns

  tautemp(1) = -src_normal(2)
  tautemp(2) = src_normal(1)

  rnormtemp(1) = src_normal(1)
  rnormtemp(2) = src_normal(2)

  beta = dreal(par0)

  ifcharge = 0
  ifdipole = 0
  ifquad = 1

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

  call mbhpotgrad2dall_cdq(beta,src,ns,ifcharge, &
       charge,ifdipole,dipstr,dipvec,ifquad,quadstr, &
       quadvec,targ,ifpot,dpot,ifgrad,dgrad,ifhess,dhess)

  val = dpot
  grad(1) = dgrad(1)
  grad(2) = dgrad(2)
  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  return
end subroutine zchker_helmstokes_s


subroutine zchcsm_modstokes_g1_nh2(par0,src,targ,src_normal, &
     pars2,val,grad,hess)
  implicit real *8 (a-h,o-z)
  ! global
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: pars2(*)
  complex *16, intent(in) :: par0
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: beta, quadvec(3), quadstr(1), dipstr(1), dipvec(2)
  real *8 :: charge(1)
  real *8 :: dpot, dgrad(2), dhess(3)
  real *8 :: tautemp(2), rnormtemp(2)
  integer :: ifcharge, ifdipole, ifquad, ns

  tautemp(1) = -src_normal(2)
  tautemp(2) = src_normal(1)

  rnormtemp(1) = src_normal(1)
  rnormtemp(2) = src_normal(2)

  beta = dreal(par0)

  ifcharge = 0
  ifdipole = 0
  ifquad = 1

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

  call mbhpotgrad2dall_cdq(beta,src,ns,ifcharge, &
       charge,ifdipole,dipstr,dipvec,ifquad,quadstr, &
       quadvec,targ,ifpot,dpot,ifgrad,dgrad,ifhess,dhess)

  val = dpot
  grad(1) = dgrad(1)
  grad(2) = dgrad(2)
  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  return
end subroutine zchcsm_modstokes_g1_nh2

subroutine zchcsm_modstokes_g2_nh2(par0,src,targ,src_normal, &
     pars2,val,grad,hess)
  implicit real *8 (a-h,o-z)
  ! global
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: pars2(*)
  complex *16, intent(in) :: par0
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: beta, quadvec(3), quadstr(1), dipstr(1), dipvec(2)
  real *8 :: charge(1)
  real *8 :: dpot, dgrad(2), dhess(3)
  real *8 :: tautemp(2), rnormtemp(2)
  real *8 :: dx, dy, pi, pi2, pi4, rr2, rr4
  integer :: ifcharge, ifdipole, ifquad, ns

  ! -grad log(r)/2pi + 2 grad perp G_nt

  tautemp(1) = -src_normal(2)
  tautemp(2) = src_normal(1)

  rnormtemp(1) = src_normal(1)
  rnormtemp(2) = src_normal(2)

  beta = dreal(par0)

  ifcharge = 0
  ifdipole = 0
  ifquad = 1

  ! 

  dx = targ(1)-src(1)
  dy = targ(2)-src(2)

  pi = 4*datan(1.d0)
  pi2 = 2.0d0*pi
  pi4 = 4.0d0*pi
  rr2 = dx**2 + dy**2
  rr4 = rr2**2

  if (rr2 .eq. 0.0d0) then
     val = 0.0d0
     grad(1) = 0.0d0
     grad(2) = 0.0d0
     hess(1,1) = 0.0d0
     hess(1,2) = 0.0d0
     hess(2,1) = 0.0d0
     hess(2,2) = 0.0d0
     return
  endif

  val = -log(rr2)/(pi4)      

  grad(1) = -dx/rr2/(pi2)
  grad(2) = -dy/rr2/(pi2)

  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  ! quadrupole corresponding to 2 G_nt

  quadstr(1) = 2.0d0
  quadvec(1) = tautemp(1)*rnormtemp(1)
  quadvec(2) = tautemp(1)*rnormtemp(2) &
       + rnormtemp(1)*tautemp(2)
  quadvec(3) = tautemp(2)*rnormtemp(2)

  ns = 1

  ifpot = 1
  ifgrad = 1
  ifhess = 0

  call mbhpotgrad2dall_cdq(beta,src,ns,ifcharge, &
       charge,ifdipole,dipstr,dipvec,ifquad,quadstr, &
       quadvec,targ,ifpot,dpot,ifgrad,dgrad,ifhess,dhess)

  val = dpot
  grad(1) = grad(1) + dgrad(2)
  grad(2) = grad(2) - dgrad(1)
  hess(1,1) = 0.0d0
  hess(1,2) = 0.0d0
  hess(2,1) = 0.0d0
  hess(2,2) = 0.0d0

  return
end subroutine zchcsm_modstokes_g2_nh2
