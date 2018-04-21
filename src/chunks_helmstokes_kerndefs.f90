
subroutine fgreensdummy(par0,src,targ, &
     pars1,pars2,val,grad,hess)
  implicit none
  complex *16, intent(in) :: par0
  complex *16, intent(in) :: pars1(*), pars2(*)
  real *8, intent(in) :: src(2), targ(2)
  complex *16, intent(out) :: val, grad(2), hess(2,2)

  return
end subroutine fgreensdummy
  
subroutine zhelmstokes_kern(zk,src,targ,&
     src_normal,targ_normal,q1,q2,fgreens,&
     pars1,pars2,val)
  !
  ! Note that this routine returns a 2 x 2
  ! matrix (as opposed to scalar output)
  ! fgreens is ignored
  !
  ! The output is a linear combo of the
  ! single and double layer potentials
  ! with coefficients cs = q1(1) and cd = q2(1)
  !
  ! If either coefficient is smaller than
  ! 1d-16, it is ignored
  !
  implicit none
  complex *16, intent(in) :: zk
  complex *16, intent(in) :: q1(*), q2(*)
  complex *16, intent(in) :: pars1(*), pars2(*)  
  real *8, intent(in) :: src(2), targ(2)
  real *8, intent(in) :: src_normal(2), targ_normal(2)  
  external fgreens
  complex *16, intent(out) :: val(2,2)
  ! local
  complex *16 :: mats(2,2), matd(2,2), cs, cd
  complex *16 :: stressmat(2,2,2), nu(2)
  integer :: ifs, ifd, ifstress
  real *8 :: small
  character :: trans
  data small / 1d-16 /

  ifs = 1
  ifd = 1
  
  mats(:,:) = 0.0d0
  matd(:,:) = 0.0d0

  cs = q1(1)
  cd = q2(1)

  if (abs(cs) .lt. small) ifs = 0
  if (abs(cd) .lt. small) ifd = 0

  if (ifs .eq. 1) then
     ifstress = 0
     call zhelmstokesletmat(zk,src,targ, &
          mats,ifstress,stressmat)
  end if
  if (ifd .eq. 1) then
     trans = 'T'
     nu(1) = src_normal(1)
     nu(2) = src_normal(2)
     call zhelmstressletmat(zk,src,targ,nu, &
          matd,trans)
  end if

  val = cs*mats + cd*matd

  return
end subroutine zhelmstokes_kern


subroutine zhelmstokeslet_p(zk,src,targ, &
     mu,p)

  implicit none
  ! global 
  complex *16, intent(in) :: zk, mu(2)
  real *8, intent(in) :: src(2), targ(2)
  complex *16, intent(out) :: p
  ! local
  real *8 rx, ry, r2, pi2

  pi2 = 8.0d0*datan(1.0d0)
  

  rx = targ(1)-src(1)
  ry = targ(2)-src(2)

  r2 = rx**2 + ry**2

  p = (-rx*mu(1)-ry*mu(2))/(r2*pi2)

  return
end subroutine zhelmstokeslet_p

subroutine zhelmstokeslet_pvec(zk,src,targ, &
     pvec)

  implicit none
  ! global 
  complex *16, intent(in) :: zk
  real *8, intent(in) :: src(2), targ(2)
  complex *16, intent(out) :: pvec(2)
  ! local
  real *8 rx, ry, r2, pi2

  pi2 = 8.0d0*datan(1.0d0)
  

  rx = targ(1)-src(1)
  ry = targ(2)-src(2)

  r2 = rx**2 + ry**2

  pvec(1) = -rx/(r2*pi2)
  pvec(2) = -ry/(r2*pi2)

  return
end subroutine zhelmstokeslet_pvec

subroutine zhelmstokeslet(zk,src,targ,mu, &
     zvel,ifstress,stress)
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  ! This function evaluates the velocity field
  ! and its derivatives (if requested)
  ! of a stokeslet with density mu
  !
  ! Let gradperp = (-dx2,dx1), and
  ! tau = (-rnorm(2),rnorm(1))
  !
  ! This function returns
  !
  ! zvel = gradperp (gradperp dot mu) G
  !
  ! and zvelder(i,j) = dxi zvel(j) (if requested)
  !
  ! Where G is the Bi-Helmholtz Green's function
  !
  ! G(x,y) = -(i*H_0^(1)(zk*r)/4 + log(r)/(2*pi))/zk^2
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  ! global
  complex *16, intent(in) :: zk, mu(2)
  real *8, intent(in) :: src(2), targ(2)
  integer, intent(in) :: ifstress
  complex *16, intent(out) :: zvel(2), stress(2,2)
  ! local
  complex *16 mat(2,2), stressmat(2,2,2)
  
  call zhelmstokesletmat(zk,src,targ, &
     mat,ifstress,stressmat)

  zvel(1) = mat(1,1)*mu(1) + mat(1,2)*mu(2)
  zvel(2) = mat(2,1)*mu(1) + mat(2,2)*mu(2)

  if (ifstress .eq. 1) then
     stress(1,1) = stressmat(1,1,1)*mu(1)+stressmat(1,1,2)*mu(2)
     stress(2,1) = stressmat(2,1,1)*mu(1)+stressmat(2,1,2)*mu(2)
     stress(1,2) = stressmat(1,2,1)*mu(1)+stressmat(1,2,2)*mu(2)
     stress(2,2) = stressmat(2,2,1)*mu(1)+stressmat(2,2,2)*mu(2)
  end if
  return
end subroutine zhelmstokeslet

subroutine zhelmstokesletmat(zk,src,targ, &
     mat,ifstress,stressmat)
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  ! global
  complex *16, intent(in) :: zk
  real *8, intent(in) :: src(2), targ(2)
  integer, intent(in) :: ifstress
  complex *16, intent(out) :: mat(2,2), stressmat(2,2,2)
  ! local
  integer :: ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
  complex *16 :: pot, grad(2), hess(3), der3(4), der4(5), der5(6)
  complex *16 :: pvec(2)
  
  ifpot = 0
  ifgrad = 0
  ifhess = 1
  ifder3 = 0
  ifder4 = 0
  ifder5 = 0
  if(ifstress.eq. 1) ifder3=1

  ! grab some derivatives of G

  call helmbhgreen_all(zk,targ,src,ifpot,pot,ifgrad,grad, &
       ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)


  ! zvel = gradperp (gradperp dot mu) G

  mat(1,1) = hess(3)
  mat(2,1) = -hess(2)
  mat(1,2) = -hess(2)
  mat(2,2) = hess(1)

  if (ifstress .eq. 1) then
     call zhelmstokeslet_pvec(zk,src,targ,pvec)

     ! stress(i,j) = -p delta(i,j) + dxi*zvel(j) + dxj*zvel(i)

     stressmat(1,1,1) = 2*der3(3) - pvec(1)
     stressmat(1,1,2) = -2*der3(2) - pvec(2)
     stressmat(2,1,1) = der3(4)-der3(2)
     stressmat(2,1,2) = der3(1)-der3(3)
     stressmat(1,2,1) = der3(4)-der3(2)
     stressmat(1,2,2) = der3(1)-der3(3)
     stressmat(2,2,1) = -2*der3(3) - pvec(1)
     stressmat(2,2,2) = 2*der3(2) - pvec(2)
  endif

  return
end subroutine zhelmstokesletmat

subroutine zhelmstresslet(zk,src,targ,mu,nu, &
     zvel,trans)
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  ! This function evaluates the velocity field
  ! of a stresslet with densities mu and nu
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  ! global
  complex *16, intent(in) :: zk, mu(2), nu(2)
  real *8, intent(in) :: src(2), targ(2)
  character, intent(in) :: trans  
  complex *16, intent(out) :: zvel(2)
  ! local
  complex *16 :: mat(2,2)

  call zhelmstressletmat(zk,src,targ,nu, &
     mat,trans)

  zvel(1) = mat(1,1)*mu(1) + mat(1,2)*mu(2)
  zvel(2) = mat(2,1)*mu(1) + mat(2,2)*mu(2)

  return
end subroutine zhelmstresslet

subroutine zhelmstressletmat(zk,src,targ,nu, &
     mat,trans)
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  ! This function produces a matrix such that
  ! mat*mu is a stresslet. trans='T' gives the
  ! transpose of the matrix
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  ! global
  complex *16, intent(in) :: zk, nu(2)
  real *8, intent(in) :: src(2), targ(2)
  character, intent(in) :: trans  
  complex *16, intent(out) :: mat(2,2)
  ! local
  integer :: ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
  complex *16 :: pot, grad(2), hess(3), der3(4), der4(5), der5(6)
  complex *16 :: mat1(2,2), mat2(2,2), tau(2), matl(2,2)
  real *8 :: gradgl(2), pi2, rx, ry, r2, c1,c2,c3

  pi2 = 8.0d0*datan(1.0d0)
  rx = targ(1)-src(1)
  ry = targ(2)-src(2)

  r2 = rx**2 + ry**2

  gradgl(1) = -rx/(r2*pi2)
  gradgl(2) = -ry/(r2*pi2)

  tau(1) = -nu(2)
  tau(2) = nu(1)
  
  ifpot = 0
  ifgrad = 0
  ifhess = 0
  ifder3 = 1
  ifder4 = 0
  ifder5 = 0

  ! grab some derivatives of G

  call helmbhgreen_all(zk,targ,src,ifpot,pot,ifgrad,grad, &
       ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)

  ! stresslet = -nu outer grad G_L - mat1 + mat2

  ! matl = -nu outer gradgl

  matl(1,1) = -nu(1)*gradgl(1)
  matl(2,1) = -nu(2)*gradgl(1)
  matl(1,2) = -nu(1)*gradgl(2)
  matl(2,2) = -nu(2)*gradgl(2)
  
  ! mat1
  ! -gradperp outer gradperp dnu G
  ! gradperp outer gradperp = (dx2dx2,-dx1dx2;-dx1dx2,dx1dx1)
  mat1(1,1) = -der3(3)*nu(1)-der3(4)*nu(2)
  mat1(2,1) = der3(2)*nu(1)+der3(3)*nu(2)
  mat1(1,2) = mat1(2,1)
  mat1(2,2) = -der3(1)*nu(1)-der3(2)*nu(2)

  ! mat2
  ! grad outer gradperp dtau G
  ! grad outer gradperp = (-dx1dx2,dx1dx1;-dx2dx2,dx1dx2)
  mat2(1,1) = -der3(2)*tau(1)-der3(3)*tau(2)
  mat2(2,1) = -der3(3)*tau(1)-der3(4)*tau(2)
  mat2(1,2) = der3(1)*tau(1)+der3(2)*tau(2)
  mat2(2,2) = der3(2)*tau(1)+der3(3)*tau(2)

  c1 = 1.0d0
  c2 = -1.0d0
  c3 = -1.0d0
  
  if (trans .eq. 'T') then
     mat(1,1) = c1*matl(1,1)+c2*mat1(1,1)+c3*mat2(1,1)
     mat(1,2) = c1*matl(2,1)+c2*mat1(2,1)+c3*mat2(2,1)
     mat(2,1) = c1*matl(1,2)+c2*mat1(1,2)+c3*mat2(1,2)
     mat(2,2) = c1*matl(2,2)+c2*mat1(2,2)+c3*mat2(2,2)
  else
     mat(1,1) = c1*matl(1,1)+c2*mat1(1,1)+c3*mat2(1,1)
     mat(1,2) = c1*matl(1,2)+c2*mat1(1,2)+c3*mat2(1,2)
     mat(2,1) = c1*matl(2,1)+c2*mat1(2,1)+c3*mat2(2,1)
     mat(2,2) = c1*matl(2,2)+c2*mat1(2,2)+c3*mat2(2,2)
  end if
     

  return
end subroutine zhelmstressletmat
