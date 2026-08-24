program dispersion_batch_test
   use iso_fortran_env, only: int64
   use wave_functions_module, only: iteratedispersion
   implicit none

   integer, parameter :: count = 12
   double precision, parameter :: px = 3.141592653589793238462643383279502884197d0
   double precision, parameter :: aphi = 1.d0/(((1.0d0 + sqrt(5.0d0))/2)+1)
   double precision, parameter :: bphi = ((1.0d0 + sqrt(5.0d0))/2)/(((1.0d0 + sqrt(5.0d0))/2)+1)
   double precision :: l0(count), estimate(count), depth(count)
   double precision :: reference(count), batched(count), lnext(count)
   logical :: active(count), wet(count)
   integer :: i, iter

   depth = [-10.d0, -1.d0, -0.1d0, 0.d0, 0.001d0, 0.01d0, &
            -1.d0, 0.25d0, 1.d0, 5.d0, 25.d0, 100.d0]
   wet = .true.
   wet(2) = .false.
   do i = 1,count
      l0(i) = 60.d0 + 2.d0*i
      estimate(i) = l0(i)
      reference(i) = estimate(i)
      if (wet(i)) then
         reference(i) = iteratedispersion(l0(i),estimate(i),px,depth(i))
         if (reference(i) < 0.d0) reference(i) = -reference(i)
      endif
   enddo

   batched = estimate
   active = wet
   iter = 0
   do while (any(active) .and. iter < 150)
      iter = iter+1
      where(active)
         lnext = l0*tanh(2*px*depth/batched)
         batched = batched*aphi + lnext*bphi
         active = abs(lnext-batched) > 0.00001d0
      endwhere
   enddo

   if (.not. any(wet .and. batched < 0.d0)) error stop 'negative-depth path not exercised'
   where(wet .and. batched < 0.d0) batched = -batched

   if (storage_size(reference(1)) /= storage_size(0_int64)) &
      error stop 'double precision is not 64 bit'
   do i = 1,count
      if (transfer(reference(i),0_int64) /= transfer(batched(i),0_int64)) then
         error stop 'batched dispersion changed a scalar result'
      endif
   enddo
end program dispersion_batch_test
