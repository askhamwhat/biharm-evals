
subroutine zhbh_stokes_matbuild(zk,wgeo,ncomp,nchs,ccs, &
  q1,q2,ntot,sysmat,ier)

  !
  ! This routine builds the matrix for solving
  ! a Dirichlet problem of the form
  !
  ! \Delta (\Delta + zk^2) u = 0 in the domain
  !
  ! with boundary conditions:
  !
  ! u = f
  ! d/dn u = g
  !
  ! on a bounded, multiply connected domain.
  ! The integral representation is based on the
  ! Stokes-type layer potentials
  !
  ! INPUT:
  !
  ! zk - COMPLEX *16, wavenumber of PDE
  ! wgeo - REAL *8 array, packed chunks description
  !        of boundary
  ! ncomp - INTEGER, number of boundary components
  !         the first boundary component should be
  !         the OUTER boundary
  ! nchs - INTEGER array, nchs(i) is the number of
  !        chunks on the ith boundary component. the
  !        chunks in wgeo should be ordered such that
  !        the first nchs(1) chunks comprise the 1st
  !        boundary component, etc.
  ! ccs - REAL *8 array, ccs(1:2,2:ncomp) describe
  !       point inside each of the connected
  !       boundary components
  ! q1 - COMPLEX *16, weight of single layer part
  !      of stokes rep
  ! q2 - COMPLEX *16, weight of double layer part
  !      of stokes rep
  ! ntot - INTEGER, size of sysmat. if ntot=-1
  !        ntot is the size ntot should be on
  !        return. Note that ntot=2*k*nch+ncomp
  !
  ! RETURN
  !
  ! ntot - INTEGER if ntot=-1 on entry, then ntot
  !        is the system size on return
  ! sysmat - COMPLEX *16, if ier = 0, then this
  !          should contain the system matrix as
  !          described in detail below.
  ! ier - INTEGER, error flag
  !        ier = 0 --> normal operation
  !        ier = 4 --> insufficient size matrix
  !        ier = 1024 --> query mode (NTOT=-1 on entry)
  !
  ! SYSTEM DESCRIPTION:
  !
  ! on simply connected q1 = 0, q2 = 1 are good
  ! choices. on multiply connected q1 = q2 = 1
  ! are good choices.
  !
  ! sysmat = ( A B )
  !          ( C D )
  !
  ! A = q1*S + q2*(D - 1/2), is the system matrix for
  ! a stokes layer potential. It is 2*number of
  ! boundary points-by-2*number of boundary points
  ! in dimension. The x and y directions of the
  ! vector field are interleaved, i.e. entries in
  ! the 2*j+1 columns correspond to the x directoin
  ! and entries in the 2*j columns correspond to the
  ! y directions.
  !
  ! B = a matrix corresponding to other parts of
  ! the representation. The first column is zero
  ! the remaining ncomp-1 columns correspond to
  ! the gradient of log charges in the holes of
  ! the domain.
  !
  ! C = a matrix corresponding to the integral of
  ! the stream function of the stokes field
  ! on each boundary component. This matrix requies some
  ! precomputation because of the potential field
  ! in the Stokes double layer representation
  !
  ! D = the integral of the constant added to
  ! the representation and the log charges along
  ! each boundary component.
  !
  ! The right hand side for this system matrix
  ! should be of the form
  !
  ! (G; I)
  !
  ! where G is the gradient of u and I is the
  ! integral of u along each boundary component
  !
  implicit none
  ! global
  complex *16 :: zk, q1, q2, sysmat(ntot,ntot)
  integer :: ncomp, nchs(*), ntot, ier
  real *8 :: wgeo(*), ccs(2,*)
  ! local
  complex *16, allocatable :: staut(:,:), sneut(:,:), &
       slay(:,:), xx(:,:), yy(:,:), stokesmat(:,:), &
       onesmat(:,:)
  real *8, allocatable :: rnorms(:,:,:), whts(:)
  integer k, nch, ichunks, iadjs, iders, iders2, ihs
  integer npts, i, j
  complex *16 :: zero
  data zero / (0.0d0,0.0d0) /

  

  return
end subroutine zhbh_stokes_matbuild


subroutine zkernel_sprimet(par0, src, targ, src_normal, &
     targ_normal, q1, q2, fgreen, pars1, pars2, val)
  implicit real *8 (a-h,o-z)
  real *8, intent(in) :: src(2), targ(2), src_normal(2)
  real *8, intent(in) :: targ_normal(2), pars1(*), pars2(*)
  complex *16, intent(in) :: par0, q1, q2
  complex *16, intent(out) :: val

  complex *16 :: val2, grad2(2), hess2(2,2)
  !
  ! this is just a wrapper routine for the kernel of the 
  ! normal derivative of the single layer (transpose)
  !

  call fgreen(par0, targ, src, pars1, pars2, val2, &
       grad2, hess2)
  val = src_normal(1)*grad2(1) + src_normal(2)*grad2(2)

  return
end subroutine zkernel_sprimet


subroutine fgreenlap(par0,src,targ,pars1,pars2,val, &
     grad,hess)
  implicit real *8 (a-h,o-z)
  real *8, intent(in) :: src(2), targ(2)
  complex *16, intent(in) :: par0, pars1(*), pars2(*)
  complex *16, intent(out) :: val, grad(2), hess(2,2)
  ! local
  real *8 :: pi2, rx, ry, pi, r2, r4, r, rx2, ry2

  pi2 = 8.0d0*atan(1.0d0)
  pi = pi2/2.0d0
  rx = targ(1)-src(1)
  ry = targ(2)-src(2)

  rx2 = rx*rx
  ry2 = ry*ry

  r2 = rx2 + ry2
  r4 = r2*r2
  r = sqrt(r2)

  val = -log(r)/pi2
  grad(1) = -rx/(r2*pi2)
  grad(2) = -ry/(r2*pi2)

  hess(1,1) = rx2/(pi*r4)-1.0d0/(pi2*r2)
  hess(2,1) = rx*ry/(pi*r4)
  hess(1,2) = hess(2,1)
  hess(2,2) = ry2/(pi*r4)-1.0d0/(pi2*r2)

  return
end subroutine fgreenlap
  
