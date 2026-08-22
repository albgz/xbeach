program vsm_test

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use vsm_interfaces, only: vsm_u_XB
   implicit none

   integer, parameter :: n = 5
   double precision :: sigz(n), uz(n), vz(n), ustz(n), nutz(n)
   double precision :: ue_sed, ve_sed
   integer :: swglm

   do swglm = 0, 1
      call vsm_u_XB(0.2d0, 0.1d0, 5.d0, 0.01d0, 1025.d0, &
                    0.1d0, 0.05d0, 0.2d0, 0.5d0, 0.8d0, &
                    0.1d0, 10.d0, 100.d0, 10.d0, 0.1d0, &
                    5.d0, 0.01d0, 0.02d0, 1.d0, swglm, n, &
                    sigz, uz, vz, ustz, nutz, ue_sed, ve_sed)

      call require(all(ieee_is_finite(sigz)), "finite sigma grid")
      call require(all(ieee_is_finite(uz)), "finite x velocity")
      call require(all(ieee_is_finite(vz)), "finite y velocity")
      call require(all(ieee_is_finite(ustz)), "finite Stokes drift")
      call require(all(ieee_is_finite(nutz)), "finite viscosity")
      call require(ieee_is_finite(ue_sed), "finite x sediment velocity")
      call require(ieee_is_finite(ve_sed), "finite y sediment velocity")
      call require(sigz(1) > 0.d0 .and. sigz(n) == 1.d0, "sigma endpoints")
      call require(all(sigz(2:n) > sigz(1:n-1)), "monotonic sigma grid")
      call require(all(nutz > 0.d0), "positive viscosity")

      if (swglm == 0) then
         call require(all(ustz == 0.d0), "defined zero Stokes output without GLM")
      else
         call require(any(ustz > 0.d0), "active Stokes profile with GLM")
      end if
   end do

   print *, "VSM path-active tests passed"

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print *, "FAIL: ", trim(label)
         error stop 1
      end if
   end subroutine require

end program vsm_test
