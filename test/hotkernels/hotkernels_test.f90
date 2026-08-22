! Standalone regression tests for the hot wave advection kernels
! (advecxho, advecyho, advecthetaho) from wave_functions_module.
!
! Each test builds a small deterministic field at production-scale
! magnitudes (J/m^2, m/s, 1 m spacing, dt ~ 0.16 s) with a wet/dry mask
! that exercises both upwind branches, then compares the kernel output
! bit-for-bit against an independent reference written out in this file.
! The reference mirrors the production arithmetic operation-by-operation,
! including the Warming-Beam correction term (whose dsz division was the
! source of the 100x GPU-chapter bug - keeping it in a CPU test guards
! against any future "simplification").
!
! Plain program, no test-framework dependency. Exit code 0 = pass.

program hotkernel_tests

   use paramsconst
   use spaceparamsdef
   use wave_functions_module

   implicit none

   integer :: failures, i, j, k

   failures = 0

   ! Production-scale grid (DELILAH uses 1 m spacing, dt ~ 0.16 s)
   call test_advecxho_upwind1
   call test_advecyho_upwind1
   call test_advecthetaho_upwind1
   call test_advecxho_warmbeam

   if (failures == 0) then
      write(*,'(a)') "ALL HOTKERNEL TESTS PASSED"
   else
      write(*,'(a,i0,a)') "FAILURES: ", failures, " test(s) failed"
      stop 1
   end if

contains

   subroutine check(name, got, want, n1, n2, n3)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: got(n1,n2,n3), want(n1,n2,n3)
      integer, intent(in) :: n1,n2,n3
      integer :: i,j,k, first_bad

      first_bad = -1
      do k=1,n3
         do j=1,n2
            do i=1,n1
               if (got(i,j,k) /= want(i,j,k)) then
                  if (first_bad < 0) first_bad = 1
                  write(*,'(a,a,a)') "  FAIL ", name, " at ("// &
                       trim(adjustl(itoa(i)))//","// &
                       trim(adjustl(itoa(j)))//","// &
                       trim(adjustl(itoa(k)))//")"
                  write(*,'(a,es18.11,a,es18.11)') "    got ", got(i,j,k), &
                       "  want ", want(i,j,k)
                  failures = failures + 1
                  return
               end if
            end do
         end do
      end do
      write(*,'(a,a,a)') "  pass ", name, " (bit-exact)"
   end subroutine check

   function itoa(n) result(s)
      integer, intent(in) :: n
      character(len=8) :: s
      write(s,'(i0)') n
   end function itoa

! ------------------------------------------------------------------
! advecxho, SCHEME_UPWIND_1
! ------------------------------------------------------------------
   subroutine test_advecxho_upwind1
      integer, parameter :: nx=30, ny=12, ntheta=6
      real*8 :: ee(nx+1,ny+1,ntheta), cgx(nx+1,ny+1,ntheta)
      real*8 :: xadvec(nx+1,ny+1,ntheta), ref(nx+1,ny+1,ntheta)
      real*8 :: dnu(nx+1,ny+1), dsu(nx+1,ny+1)
      real*8 :: dsz(nx+1,ny+1), dsdnzi(nx+1,ny+1), flux(nx+1,ny+1)
      integer :: wete(nx+1,ny+1)
      integer :: i,j,k
      real*8 :: cgxu_local

      do k=1,ntheta
         do j=1,ny+1
            do i=1,nx+1
               ee(i,j,k)  = 5.0d0 + 0.1d0*i + 0.05d0*j + 0.01d0*k
               ! mix of + and - cgx so both upwind branches run
               cgx(i,j,k) = 1.5d0 + 0.2d0*mod(i+k,3) - 0.4d0*mod(i-j,2)
            end do
         end do
      end do
      dnu    = 1.0d0/0.16d0          ! 1/dt-like coefficient
      dsu    = 1.0d0
      dsz    = 1.0d0
      dsdnzi = 1.0d0
      wete   = 1
      ! some dry cells, including a dry column so the flux difference sees 0
      wete(5,4) = 0
      wete(6,4) = 0
      wete(10,7) = 0
      wete(1,2) = 0

      jmin_ee = 1
      jmax_ee = ny+1

      call advecxho(ee,cgx,xadvec,nx,ny,ntheta,dnu,dsu,dsdnzi, &
                    SCHEME_UPWIND_1, wete, 0.16d0, dsz)

      ! independent reference (same operation order as production)
      flux = 0.0d0
      ref  = 0.0d0
      do k=1,ntheta
         do j=1,ny+1
            do i=1,nx
               if (wete(i,j)==1) then
                  cgxu_local = 0.5d0*(cgx(i+1,j,k)+cgx(i,j,k))
                  if (cgxu_local>0) then
                     flux(i,j) = ee(i,j,k)*cgxu_local*dnu(i,j)
                  else
                     flux(i,j) = ee(i+1,j,k)*cgxu_local*dnu(i,j)
                  end if
               end if
            end do
         end do
         do j=jmin_ee,jmax_ee
            do i=2,nx
               if (wete(i,j)==1) then
                  ref(i,j,k) = (flux(i,j)-flux(i-1,j))*dsdnzi(i,j)
               end if
            end do
         end do
      end do

      call check("advecxho upwind1", xadvec, ref, nx+1, ny+1, ntheta)
   end subroutine test_advecxho_upwind1

! ------------------------------------------------------------------
! advecyho, SCHEME_UPWIND_1
! ------------------------------------------------------------------
   subroutine test_advecyho_upwind1
      integer, parameter :: nx=30, ny=12, ntheta=6
      real*8 :: ee(nx+1,ny+1,ntheta), cgy(nx+1,ny+1,ntheta)
      real*8 :: yadvec(nx+1,ny+1,ntheta), ref(nx+1,ny+1,ntheta)
      real*8 :: dsv(nx+1,ny+1), dnv(nx+1,ny+1)
      real*8 :: dnz(nx+1,ny+1), dsdnzi(nx+1,ny+1), flux(nx+1,ny+1)
      integer :: wete(nx+1,ny+1)
      integer :: i,j,k
      real*8 :: cgyv_local

      do k=1,ntheta
         do j=1,ny+1
            do i=1,nx+1
               ee(i,j,k)  = 4.0d0 + 0.15d0*j + 0.02d0*i + 0.03d0*k
               cgy(i,j,k) = -1.2d0 + 0.3d0*mod(j+k,3)
            end do
         end do
      end do
      dsv    = 1.0d0
      dnv    = 1.0d0/0.16d0
      dnz    = 1.0d0
      dsdnzi = 1.0d0
      wete   = 1
      wete(:,6) = 0      ! dry band
      wete(:,7) = 0

      call advecyho(ee,cgy,yadvec,nx,ny,ntheta,dsv,dnv,dsdnzi, &
                    SCHEME_UPWIND_1, wete, 0.16d0, dnz)

      flux = 0.0d0
      ref  = 0.0d0
      do k=1,ntheta
         do j=1,ny
            do i=1,nx+1
               if (wete(i,j)==1) then
                  cgyv_local = 0.5d0*(cgy(i,j+1,k)+cgy(i,j,k))
                  if (cgyv_local>0) then
                     flux(i,j) = ee(i,j,k)*cgyv_local*dsv(i,j)
                  else
                     flux(i,j) = ee(i,j+1,k)*cgyv_local*dsv(i,j)
                  end if
               end if
            end do
         end do
         do j=2,ny
            do i=1,nx+1
               if (wete(i,j)==1) then
                  ref(i,j,k) = (flux(i,j)-flux(i,j-1))*dsdnzi(i,j)
               end if
            end do
         end do
      end do

      call check("advecyho upwind1", yadvec, ref, nx+1, ny+1, ntheta)
   end subroutine test_advecyho_upwind1

! ------------------------------------------------------------------
! advecthetaho, SCHEME_UPWIND_1
! ------------------------------------------------------------------
   subroutine test_advecthetaho_upwind1
      integer, parameter :: nx=20, ny=10, ntheta=8
      real*8 :: ee(nx+1,ny+1,ntheta), cth(nx+1,ny+1,ntheta)
      real*8 :: tadvec(nx+1,ny+1,ntheta), ref(nx+1,ny+1,ntheta)
      integer :: wete(nx+1,ny+1)
      real*8 :: fluxt(ntheta)
      integer :: i,j,k
      real*8 :: ctb_local
      real*8, parameter :: dtheta = 0.1d0

      do k=1,ntheta
         do j=1,ny+1
            do i=1,nx+1
               ee(i,j,k)  = 2.0d0 + 0.1d0*k + 0.05d0*mod(i,5)
               cth(i,j,k) = 0.3d0 + 0.1d0*mod(i+k,3) - 0.2d0*mod(j,2)
            end do
         end do
      end do
      wete = 1
      wete(4,4) = 0
      wete(12,8) = 0

      call advecthetaho(ee,cth,tadvec,nx,ny,ntheta,dtheta, &
                        SCHEME_UPWIND_1, wete)

      ref = 0.0d0
      do j=1,ny+1
         do i=1,nx+1
            if (wete(i,j)==1) then
               fluxt = 0.0d0
               do k=1,ntheta-1
                  ctb_local = 0.5d0*(cth(i,j,k+1)+cth(i,j,k))
                  if (ctb_local>0) then
                     fluxt(k) = ee(i,j,k)*ctb_local
                  else
                     fluxt(k) = ee(i,j,k+1)*ctb_local
                  end if
               end do
               ref(i,j,1)      = (fluxt(1)-0.0d0)/dtheta
               do k=2,ntheta-1
                  ref(i,j,k)   = (fluxt(k)-fluxt(k-1))/dtheta
               end do
               ref(i,j,ntheta) = (0.0d0-fluxt(ntheta-1))/dtheta
            end if
         end do
      end do

      call check("advecthetaho upwind1", tadvec, ref, nx+1, ny+1, ntheta)
   end subroutine test_advecthetaho_upwind1

! ------------------------------------------------------------------
! advecxho, SCHEME_WARMBEAM (includes the Warming-Beam correction)
! ------------------------------------------------------------------
   subroutine test_advecxho_warmbeam
      integer, parameter :: nx=30, ny=12, ntheta=6
      real*8 :: ee(nx+1,ny+1,ntheta), cgx(nx+1,ny+1,ntheta)
      real*8 :: xadvec(nx+1,ny+1,ntheta), ref(nx+1,ny+1,ntheta)
      real*8 :: dnu(nx+1,ny+1), dsu(nx+1,ny+1)
      real*8 :: dsz(nx+1,ny+1), dsdnzi(nx+1,ny+1), flux(nx+1,ny+1)
      integer :: wete(nx+1,ny+1)
      integer :: i,j,k
      real*8 :: cgxu_local, eupw_local
      real*8, parameter :: dt = 0.16d0

      do k=1,ntheta
         do j=1,ny+1
            do i=1,nx+1
               ee(i,j,k)  = 3.0d0 + 0.2d0*i + 0.07d0*j + 0.011d0*k
               cgx(i,j,k) = 2.0d0 + 0.3d0*mod(i,3) - 0.5d0*mod(j+k,2)
            end do
         end do
      end do
      dsu    = 1.0d0
      dnu    = 1.0d0/dt
      dsz    = 1.0d0
      dsdnzi = 1.0d0
      wete   = 1
      wete(3,3) = 0
      wete(20,9) = 0

      jmin_ee = 1
      jmax_ee = ny+1

      call advecxho(ee,cgx,xadvec,nx,ny,ntheta,dnu,dsu,dsdnzi, &
                    SCHEME_WARMBEAM, wete, dt, dsz)

      ! independent reference. Serial build: xmpi_istop/xmpi_isbot are
      ! .true., so boundary columns i=1 and i=nx are computed too.
      flux = 0.0d0
      ref  = 0.0d0
      do k=1,ntheta
         do j=1,ny+1
            do i=2,nx-1
               if (wete(i,j)==1) then
                  cgxu_local = 0.5d0*(cgx(i+1,j,k)+cgx(i,j,k))
                  if (cgxu_local>0) then
                     eupw_local = ((dsu(i-1,j)+0.5d0*dsu(i,j))*ee(i,j,k) &
                                   -0.5d0*dsu(i,j)*ee(i-1,j,k))/dsu(i-1,j)
                     if (eupw_local<0.0d0) eupw_local = ee(i,j,k)
                     flux(i,j) = eupw_local*cgxu_local*dnu(i,j)
                  else
                     eupw_local = ((dsu(i+1,j)+0.5d0*dsu(i,j))*ee(i+1,j,k) &
                                   -0.5d0*dsu(i,j)*ee(i+2,j,k))/dsu(i+1,j)
                     if (eupw_local<0.0d0) eupw_local = ee(i+1,j,k)
                     flux(i,j) = eupw_local*cgxu_local*dnu(i,j)
                  end if
               end if
            end do
            i = 1
            if (wete(i,j)==1) then
               cgxu_local = 0.5d0*(cgx(i+1,j,k)+cgx(i,j,k))
               if (cgxu_local>0) then
                  flux(i,j) = ee(i,j,k)*cgxu_local*dnu(i,j)
               else
                  eupw_local = ((dsu(i+1,j)+0.5d0*dsu(i,j))*ee(i+1,j,k) &
                                -0.5d0*dsu(i,j)*ee(i+2,j,k))/dsu(i+1,j)
                  if (eupw_local<0.0d0) eupw_local = ee(i+1,j,k)
                  flux(i,j) = eupw_local*cgxu_local*dnu(i,j)
               end if
            end if
            i = nx
            if (wete(i,j)==1) then
               cgxu_local = 0.5d0*(cgx(i+1,j,k)+cgx(i,j,k))
               if (cgxu_local>0) then
                  eupw_local = ((dsu(i-1,j)+0.5d0*dsu(i,j))*ee(i,j,k) &
                                -0.5d0*dsu(i,j)*ee(i-1,j,k))/dsu(i-1,j)
                  if (eupw_local<0.0d0) eupw_local = ee(i,j,k)
                  flux(i,j) = eupw_local*cgxu_local*dnu(i,j)
               else
                  flux(i,j) = ee(i+1,j,k)*cgxu_local*dnu(i,j)
               end if
            end if
         end do
         do j=jmin_ee,jmax_ee
            do i=2,nx
               if (wete(i,j)==1) then
                  ref(i,j,k) = (flux(i,j)-flux(i-1,j))*dsdnzi(i,j)
               end if
            end do
         end do
      end do
      ! Warming-Beam correction (the term the GPU chapter got wrong 100x)
      do k=1,ntheta
         do j=jmin_ee,jmax_ee
            do i=2,nx
               if (wete(i,j)==1) then
                  ref(i,j,k) = ref(i,j,k) &
                     -((ee(i+1,j,k)-ee(i  ,j,k))/dsu(i  ,j) &
                       -(ee(i  ,j,k)-ee(i-1,j,k))/dsu(i-1,j))/ &
                       dsz(i,j)*dt/2*cgx(i,j,k)**2
               end if
            end do
         end do
      end do

      call check("advecxho warmbeam", xadvec, ref, nx+1, ny+1, ntheta)
   end subroutine test_advecxho_warmbeam

end program hotkernel_tests
