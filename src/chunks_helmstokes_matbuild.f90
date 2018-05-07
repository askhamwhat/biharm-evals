
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
  complex *16, allocatable :: temp1(:,:),temp2(:,:), &
       xx(:,:), yy(:,:), stokesmat(:,:), work(:), xxlong(:,:), &
       yylongtrans(:,:), logval(:,:), logvalperp(:,:,:), brmat(:,:)
  real *8, allocatable :: rnorms(:,:), whts(:)
  real *8 :: eps, errs(1000), p1,p2,p3,p4
  integer job, lw, ndim, ngmrec, numit, niter, j1
  integer k, nch, ichunks, iadjs, iders, iders2, ihs
  integer npts, i, j, iii, incx, incy, ind, npts2, jj, ii
  complex *16 :: zero, alpha, beta, pars1, pars2, one
  complex *16 :: val, grad(2), hess(2,2)
  data zero, one / (0.0d0,0.0d0), (1.0d0,0.0d0) /
  external zkernel_slp, zkernel_sprime
  external fgreenlap, chs_multatrans, fgreensdummy
  external zhelmstokes_stream_kern, zhelmstokes_kern

  ! get dimensions
  call chunkunpack1(wgeo,k,nch,ichunks,iadjs, &
       iders,iders2,ihs)
  
  ! normals and smooth integration weights
  allocate(rnorms(2,k*nch),whts(k*nch))
  call chunknormals(wgeo,rnorms)
  call chunkwhts(k,nch,wgeo(ichunks),wgeo(iders), &
       wgeo(ihs),whts)
  
  npts = k*nch

  call prinf('npts*',npts,1)

  allocate(xx(npts,ncomp),yy(npts,ncomp))

  do i = 1,ncomp
     do j = 1,npts
        xx(j,i) = zero
        yy(j,i) = zero
     enddo
  enddo

  ! set up integration weights on separate components

  ind = 1
  do iii = 1,ncomp
     call prinf('ind *',ind,1)
     do i = 1,nchs(iii)
        do j = 1,k
           xx(ind,iii) = whts(ind)
           ind = ind+1
        enddo
     enddo
  enddo

  call prinf('ind *',ind,1)

  ! charges and constants

  allocate(logval(npts,ncomp),logvalperp(2,npts,ncomp), &
       brmat(ncomp,ncomp))

  ! constant
  
  do j = 1,npts
     logval(j,1) = one
     logvalperp(1,j,1) = zero
     logvalperp(2,j,1) = zero
  enddo

  ! log charge
  
  do i = 2,ncomp
     do j = 1,npts
        call fgreenlap(zk,ccs(1,i),wgeo(ichunks+2*(j-1)), &
             pars1,pars2,val,grad,hess)
        logval(j,i) = val
        logvalperp(1,j,i) = -grad(2)
        logvalperp(2,j,i) = grad(1)
     enddo
  enddo

  ! bottom right matrix is integral of these charges
  ! on each component

  beta = 0.0d0
  incx = 1
  incy = 1
  alpha = 1.0d0

  call zgemm('T','N',ncomp,ncomp,npts,alpha,xx, &
       npts,logval,npts,beta,brmat,ncomp)

  ! copy into system matrix
  ! (bottom right)

  do i = 1,ncomp
     do j = 1,ncomp
        sysmat(2*npts+j,2*npts+i) = brmat(j,i)
     enddo
  enddo

  ! copy velocity due to charges into system matrix
  ! (top right)

  do i = 1,ncomp
     ii = i+2*npts
     do j = 1,npts
        jj = 2*j-1
        sysmat(jj,ii) = logvalperp(1,j,i)
        jj = 2*j
        sysmat(jj,ii) = logvalperp(2,j,i)
     enddo
  enddo

  ! single layer evaluation matrix in temp1
  allocate(temp1(npts,npts),temp2(npts,npts))
  call zbuildmat(k,wgeo,zkernel_slp,q1,q2, &
       fgreenlap,zk,pars1,pars2,npts,temp1)

  beta = 0.0d0
  incx = 1
  incy = 1
  alpha = 1.0d0

  ! apply transpose to weights vectors (1 per
  ! boundary component)
  call zgemm('T','N',npts,ncomp,npts,alpha,temp1, &
       npts,xx,npts,beta,yy,npts)

  ! matrix for neumann problem in temp2

  call zbuildmat(k,wgeo,zkernel_sprime,q1,q2, &
       fgreenlap,zk,pars1,pars2,npts,temp2)
  do i = 1,npts
     temp2(i,i) = temp2(i,i) + 0.5d0
  enddo
  do j = 1,npts
     do i = 1,npts
        temp2(i,j) = temp2(i,j) + whts(j)
     enddo
  enddo

  ! solve transpose problem for each yy vector
  
  ngmrec = 200
  numit = 200
  lw = (ngmrec*2 + 4)*npts
  allocate(work(lw))
  job = 0
  ier = 0
  eps = 1.0d-14

  do i = 1,ncomp
     call prinf('solving neumann problem, comp *',i,1)

     call cgmres(ier,npts,temp2,chs_multatrans, &
          p1,p2,p3,p4,yy(1,i), &
          eps,numit,xx(1,i),niter,errs,ngmrec,work)
     call prin2('errs *',errs,niter)
  enddo

  ! form S_tau

  call prinf('forming stau ...*',i,0)
  do i = 1,npts
     ndim = 2
     call chunkderf(temp1(1,i),temp2(1,i),ndim,k,nch, &
          wgeo(ichunks),wgeo(iders),wgeo(ihs))
  enddo
  call prinf('done*',i,0)

  ! apply Stau transpose to solutions of transpose problem
  beta = 0.0d0
  incx = 1
  incy = 1
  alpha = 1.0d0

  call zgemm('T','N',npts,ncomp,npts,alpha,temp2, &
       npts,xx,npts,beta,yy,npts)

  ! get whts^T*stream matrix

  npts2 = 2*npts
  allocate(xxlong(npts2,ncomp),yylongtrans(ncomp,npts2))
  allocate(stokesmat(npts2,npts2))

  ndim = 2
  call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_stream_kern, &
       q1,q2,fgreensdummy,zk,pars1,pars2,npts2, &
       stokesmat)

  do i = 1,ncomp
     do j = 1,npts2
        xxlong(j,i) = zero
     enddo
  enddo

  do i = 1,npts2
     do j = 1,ncomp
        yylongtrans(j,i) = zero
     enddo
  enddo

  ! set up integration weights on separate components

  ind = 1
  do iii = 1,ncomp
     call prinf('ind *',ind,1)
     do i = 1,nchs(iii)
        do j = 1,k
           xxlong(2*ind-1,iii) = whts(ind)
           ind = ind+1
        enddo
     enddo
  enddo

  ! integrate stream function matrix
  
  call zgemm('T','N',ncomp,npts2,npts2,alpha,xxlong, &
       npts2,stokesmat,npts2,beta,yylongtrans,ncomp)
  
  ! add these components together and store in
  ! system matrix (bottom left)

  do j = 1,npts
     do j1 = 1,2
        do i = 1,ncomp
           jj = 2*j-2+j1
           sysmat(npts2+i,jj) = yylongtrans(i,jj) - &
                q2*yy(j,i)*rnorms(j1,j)
        enddo
     enddo
  enddo

  ! get stokes system matrix (top left)

  ndim = 2
  call zbuildmat_vec(ndim,k,wgeo,zhelmstokes_kern, &
       q1,q2,fgreensdummy,zk,pars1,pars2,npts2, &
       stokesmat)
  
  ! copy in (top left)

  do i = 1,npts2
     do j = 1,npts2
        sysmat(j,i) = stokesmat(j,i)
     enddo
  enddo

  ! add onesmat modification
  
  call normalonesmat(stokesmat,whts,rnorms,k,nch)
  do j = 1,npts2
     do i = 1,npts2
        sysmat(i,j) = sysmat(i,j) + stokesmat(i,j)
        if (i .eq. j) sysmat(i,j) = sysmat(i,j)-0.5d0*q2
     enddo
  enddo
  
  
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
  
subroutine chs_multatrans(a,p1,p2,p3,p4,x,y,n)
  implicit real *8 (a-h,o-z)
  !
  !       computes the product a^T*x = y
  !
  complex *16 a(n,*),p1,p2,p3,p4,x(*),y(*)
  complex *16 alpha, beta

  beta = 0.0d0
  incx = 1
  incy = 1
  alpha = 1.0d0

  call zgemv ('T', n, n, alpha, a, n, x, incx, beta, y, incy)

  return
end subroutine chs_multatrans
