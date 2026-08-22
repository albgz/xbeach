! Standalone regression tests, part 2: the state-carrying hot kernels.
!
!   1) iteratedispersion        - dispersion-relation fixed-point
!      (elemental function; tested against the fixed-point residual and
!       deep/shallow-water limits)
!   2) compute_wave_direction_velocities - wave-current velocity split,
!      refraction velocity and the Dano clip (elementwise; exact
!      hand-computed reference, bit-for-bit)
!   3) wave_dispersion (wci = 0) - linear dispersion time step:
!      sigma invariant, L1 fixed-point closure, k = 2px/L1, lateral
!      boundary copies (all-true in the serial build), determinism
!
! Serial build: xmpi_isleft/isright/istop/isbot are all .true.
! Build/run: see the Makefile in this directory.

program hotkernels2

   use iso_fortran_env, only: int64
   use params
   use spaceparams
   use wave_functions_module

   implicit none

   integer, save :: nfail = 0

   ! shared state instance, wired up by init_state
   type(spacepars) :: s

   ! --- scalars (targets for the spacepars scalar pointers) --------------
   integer, target :: t_nx = 3, t_ny = 2, t_nth = 3, t_nts = 2

   ! --- dummy storage for the spacepars pointer components --------------
   integer,             allocatable,target :: wete(:,:)
   real*8, allocatable,target :: cg(:,:), c(:,:), sigm(:,:), hhw(:,:)
   real*8, allocatable,target :: u(:,:), v(:,:), k(:,:), nn(:,:), l1(:,:)
   real*8, allocatable,target :: thetamean(:,:), alfaz(:,:), wm(:,:), hbar(:,:)
   real*8, allocatable,target :: costh(:,:,:), sinth(:,:,:), sigt(:,:,:)
   real*8, allocatable,target :: cgx(:,:,:), cgy(:,:,:), ctheta(:,:,:)
   real*8, allocatable,target :: cx(:,:,:), cy(:,:,:)
   real*8, allocatable :: ref_cgx(:,:,:), ref_cgy(:,:,:), ref_cx(:,:,:)
   real*8, allocatable :: ref_cy(:,:,:), ref_cth(:,:,:)
   real*8, allocatable :: sigm_copy(:,:)
   real*8, allocatable :: dhdx2(:,:), dhdy2(:,:), dudx2(:,:), dudy2(:,:)
   real*8, allocatable :: dvdx2(:,:), dvdy2(:,:), s2k2(:,:)
   real*8 :: s2kv
   real*8 :: L0v, L, resid, deepL, shalL
   integer :: i, j, kk
   real*8 :: clip

   if (storage_size(0.0d0) /= 64) error stop "tests require 64-bit double precision"

   call hdr("1) iteratedispersion")

   ! deep water: tanh(2*pi*h/L0) -> 1, so L -> L0
   L0v = 100.d0
   L = iteratedispersion(L0v, 400.d0, 4.d0*atan(1.d0), 500.d0)
   deepL = L
   resid = abs(L - L0v*tanh(2.d0*(4.d0*atan(1.d0))*500.d0/L))
   call check("deep water: fixed-point residual < 1e-5", resid < 1.0d-5)
   call check("deep water: L ~ L0 (within 0.1%)", abs(L-L0v) < 1.0d-3*L0v)

   ! shallow water: L << L0
   L = iteratedispersion(L0v, 10.d0, 4.d0*atan(1.d0), 6.d0)
   shalL = L
   resid = abs(L - L0v*tanh(2.d0*(4.d0*atan(1.d0))*6.d0/L))
   call check("shallow water: fixed-point residual < 1e-5", resid < 1.0d-5)
   call check("shallow water: L < L0", shalL < L0v)
   call check("shallow water: L > 2*pi*h*0.9", shalL > 0.9d0*(4.d0*atan(1.d0))*6.d0)

   ! determinism
   call check("deterministic (bit-exact repeat)", bit_equal( &
        iteratedispersion(100.d0, 400.d0, 4.d0*atan(1.d0), 500.d0), deepL))

   call hdr("2) compute_wave_direction_velocities (wci = 0)")

   call init_state
   ! allocate the reference arrays (the state itself is wired in init_state)
   allocate(ref_cgx(t_nx+1,t_ny+1,t_nth)); allocate(ref_cgy(t_nx+1,t_ny+1,t_nth))
   allocate(ref_cx(t_nx+1,t_ny+1,t_nth));  allocate(ref_cy(t_nx+1,t_ny+1,t_nth))
   allocate(ref_cth(t_nx+1,t_ny+1,t_nth))
   ! wet everywhere except (2,2)
   wete(:,:) = 1
   wete(2,2) = 0
   do kk = 1, t_nth
      do j = 1, t_ny+1
         do i = 1, t_nx+1
            costh(i,j,kk) = 0.2d0*i + 0.1d0*j - 0.05d0*kk
            sinth(i,j,kk) = 0.1d0*kk + 0.05d0*j - 0.02d0*i
         enddo
      enddo
   enddo
   do j = 1, t_ny+1
      do i = 1, t_nx+1
         cg(i,j)   = 1.5d0 + 0.1d0*i
         c(i,j)    = 1.0d0 + 0.05d0*j
         sigm(i,j) = 2.0d0 + 0.25d0*i
      enddo
   enddo

   allocate(dhdx2(t_nx+1,t_ny+1)); allocate(dhdy2(t_nx+1,t_ny+1))
   allocate(dudx2(t_nx+1,t_ny+1)); allocate(dudy2(t_nx+1,t_ny+1))
   allocate(dvdx2(t_nx+1,t_ny+1)); allocate(dvdy2(t_nx+1,t_ny+1))
   allocate(s2k2(t_nx+1,t_ny+1))
   dhdx2(:,:) = 0.1d0;  dhdy2(:,:) = -0.2d0
   dudx2(:,:) = 0.d0;   dudy2(:,:) = 0.d0
   dvdx2(:,:) = 0.d0;   dvdy2(:,:) = 0.d0
   s2kv = 0.2d0                 ! small sinh(2kh): forces the Dano clip on
   s2k2(:,:) = s2kv
   s2k2(2,2) = 0.d0             ! dry-cell denominator must never be evaluated
   clip = 0.5d0*(4.d0*atan(1.d0))/10.d0    ! Trep = 10

   ! exact reference (same operation order as the production code)
   do kk = 1, t_nth
      do j = 1, t_ny+1
         do i = 1, t_nx+1
            if (wete(i,j) == 1) then
               ref_cgx(i,j,kk) = cg(i,j)*costh(i,j,kk)
               ref_cgy(i,j,kk) = cg(i,j)*sinth(i,j,kk)
               ref_cx(i,j,kk)  = c(i,j)*costh(i,j,kk)
               ref_cy(i,j,kk)  = c(i,j)*sinth(i,j,kk)
               L = sigm(i,j)/s2k2(i,j)*(dhdx2(i,j)*sinth(i,j,kk)-dhdy2(i,j)*costh(i,j,kk))
               ref_cth(i,j,kk) = sign(1.d0,L)*min(abs(L),clip)
            else
               ref_cgx(i,j,kk) = 0.d0
               ref_cgy(i,j,kk) = 0.d0
               ref_cx(i,j,kk)  = 0.d0
               ref_cy(i,j,kk)  = 0.d0
               ref_cth(i,j,kk) = 0.d0
            endif
         enddo
      enddo
   enddo

   call compute_wave_direction_velocities(s, par_wci(0), 0, dhdx2, dhdy2, &
        dudx2, dudy2, dvdx2, dvdy2, s2k2)

   call cmp3("cgx",  cgx,  ref_cgx)
   call cmp3("cgy",  cgy,  ref_cgy)
   call cmp3("cx",   cx,   ref_cx)
   call cmp3("cy",   cy,   ref_cy)
   call cmp3("ctheta (Dano-clipped)", ctheta, ref_cth)

   ! the Dano clip must actually be active somewhere (data above gives
   ! raw refraction > clip at e.g. (4,3,k=1))
   call check("Dano clip active (some wet point reaches the limit)", &
        maxval(abs(ctheta(:,:,:))) >= clip - 1.d-15)
   do kk = 1, t_nth
      call check("dry cell zero: cgx(2,2,"//trim(it2(kk))//")", bit_equal(cgx(2,2,kk), 0.d0))
      call check("dry cell zero: ctheta(2,2,"//trim(it2(kk))//")", bit_equal(ctheta(2,2,kk), 0.d0))
   enddo

   ! --- wci = 1 with zero currents must equal wci = 0 (bit-for-bit) -----
   ! (the wci = 0 results are still in cgx/cgy/cx/cy/ctheta)
   call hdr("2b) wci = 1, u = v = 0 (equivalence with wci = 0)")
   ref_cgx(:,:,:) = cgx; ref_cgy(:,:,:) = cgy; ref_cx(:,:,:) = cx
   ref_cy(:,:,:)  = cy;  ref_cth(:,:,:) = ctheta
   call compute_wave_direction_velocities(s, par_wci(1), 0, dhdx2, dhdy2, &
        dudx2, dudy2, dvdx2, dvdy2, s2k2)
   call cmp3("cgx   wci1(u=v=0) == wci0", cgx,  ref_cgx)
   call cmp3("cgy   wci1(u=v=0) == wci0", cgy,  ref_cgy)
   call cmp3("cx    wci1(u=v=0) == wci0", cx,   ref_cx)
   call cmp3("cy    wci1(u=v=0) == wci0", cy,   ref_cy)
   call cmp3("ctheta wci1(u=v=0) == wci0", ctheta, ref_cth)

   ! --- wci = 1 with non-zero currents and velocity gradients -----------
   ! This is the branch most likely to drift during the elementwise rewrite.
   call hdr("2c) wci = 1, non-zero currents (exact scalar reference)")
   do j = 1, t_ny+1
      do i = 1, t_nx+1
         u(i,j) = 0.03d0*i - 0.02d0*j
         v(i,j) = -0.01d0*i + 0.04d0*j
         dudx2(i,j) = 0.001d0*i
         dudy2(i,j) = -0.002d0*j
         dvdx2(i,j) = 0.003d0*j
         dvdy2(i,j) = -0.004d0*i
      enddo
   enddo
   do kk = 1, t_nth
      do j = 1, t_ny+1
         do i = 1, t_nx+1
            if (wete(i,j) == 1) then
               ref_cgx(i,j,kk) = cg(i,j)*costh(i,j,kk)+u(i,j)
               ref_cgy(i,j,kk) = cg(i,j)*sinth(i,j,kk)+v(i,j)
               ref_cx(i,j,kk)  = c(i,j)*costh(i,j,kk)+u(i,j)
               ref_cy(i,j,kk)  = c(i,j)*sinth(i,j,kk)+v(i,j)
               L = sigm(i,j)/s2k2(i,j)*(dhdx2(i,j)*sinth(i,j,kk)-dhdy2(i,j)*costh(i,j,kk)) + &
                    (costh(i,j,kk)*(sinth(i,j,kk)*dudx2(i,j)-costh(i,j,kk)*dudy2(i,j)) + &
                     sinth(i,j,kk)*(sinth(i,j,kk)*dvdx2(i,j)-costh(i,j,kk)*dvdy2(i,j)))
               ref_cth(i,j,kk) = sign(1.d0,L)*min(abs(L),clip)
            else
               ref_cgx(i,j,kk) = 0.d0
               ref_cgy(i,j,kk) = 0.d0
               ref_cx(i,j,kk)  = 0.d0
               ref_cy(i,j,kk)  = 0.d0
               ref_cth(i,j,kk) = 0.d0
            endif
         enddo
      enddo
   enddo
   call compute_wave_direction_velocities(s, par_wci(1), 0, dhdx2, dhdy2, &
        dudx2, dudy2, dvdx2, dvdy2, s2k2)
   call cmp3("cgx   wci1(non-zero)", cgx, ref_cgx)
   call cmp3("cgy   wci1(non-zero)", cgy, ref_cgy)
   call cmp3("cx    wci1(non-zero)", cx, ref_cx)
   call cmp3("cy    wci1(non-zero)", cy, ref_cy)
   call cmp3("ctheta wci1(non-zero)", ctheta, ref_cth)

   call hdr("3) wave_dispersion (wci = 0), two steps")

   ! reset the state fields dispersion touches (keep the wci fields)
   l1(:,:)  = 0.d0
   k(:,:)   = 0.d0
   nn(:,:)  = 0.d0
   cg(:,:)  = 0.d0
   sigm(:,:) = 0.d0
   sigt(:,:,:) = 0.d0
   thetamean(:,:) = 0.d0
   alfaz(:,:) = 0.d0
   wm(:,:) = 0.d0
   do j = 1, t_ny+1
      do i = 1, t_nx+1
         wete(i,j) = 1
         hhw(i,j)  = 10.d0 + 0.5d0*i + 0.25d0*j     ! gently varying, all wet
      enddo
   enddo

   call wave_dispersion(s, par_disp())
   call wave_dispersion(s, par_disp())      ! second step: save-state steady
   allocate(sigm_copy(size(l1,1),size(l1,2)))

   ! fixed-Trep invariant: sigt = 2px/Trep (wci = 0 branch, Trep changed
   ! from 0 to 10 on the first step)
   call check("sigt == 2*px/Trep (maxdev = "// &
        trim(es2(maxval(abs(sigt - 2.d0*(4.d0*atan(1.d0))/10.d0))))//")", &
        maxval(abs(sigt - 2.d0*(4.d0*atan(1.d0))/10.d0)) < 1.0d-12)

   ! L1 lateral (y) boundary copies in the serial build (all flags true)
   call check("L1 lateral bc: L1(:,1) == L1(:,2)",  all(bit_equal(l1(:,1), l1(:,2))))
   call check("L1 lateral bc: L1(:,"//trim(it2(t_ny+1))//") == L1(:,"//trim(it2(t_ny))//")", &
        all(bit_equal(l1(:,t_ny+1), l1(:,t_ny))))

   ! interior-row fixed-point closure: L1 = L0 tanh(2 pi h / L1),
   ! L0 = g Trep^2 / (2 px). Only interior y-rows (2..ny) are pure fixed
   ! points; rows 1 and ny+1 were overwritten by the lateral copies.
   resid = 0.d0
   do i = 1, t_nx+1
      L = (9.81d0*10.d0**2)/(2.d0*(4.d0*atan(1.d0)))
      resid = max(resid, abs(l1(i,2) - L*tanh(2.d0*(4.d0*atan(1.d0))*hhw(i,2)/l1(i,2))))
   enddo
   call check("L1 fixed-point closure interior row (max resid = "//trim(es2(resid))//")", &
        resid < 1.0d-5)

   ! k = 2px/L1 (same formula, bit-for-bit; all cells wet here)
   call check("k == 2*px/L1 bit-exact", &
        all(bit_equal(k(:,:), 2.d0*(4.d0*atan(1.d0))/l1(:,:))))

   ! determinism: L1 is recomputed fresh every step from the same inputs,
   ! so a third step must reproduce it bit-for-bit
   l1 = l1 + 0.d0
   sigm_copy(:,:) = l1
   call wave_dispersion(s, par_disp())
   call check("step-to-step determinism (L1 bit-identical)", all(bit_equal(l1, sigm_copy)))

   call summary
   stop 0

contains

   pure elemental logical function bit_equal(a, b)
   real*8, intent(in) :: a, b
   integer(int64) :: a_bits, b_bits
      a_bits = transfer(a, a_bits)
      b_bits = transfer(b, b_bits)
      bit_equal = a_bits == b_bits
   end function bit_equal

   !--------------------------------------------------------------------
   function par_wci(wci) result(p)
   integer, intent(in) :: wci
   type(parameters) :: p
      p%wci        = wci
      p%single_dir = 0
      p%px         = 4.d0*atan(1.d0)
      p%Trep       = 10.d0
      p%g          = 9.81d0
      p%eps        = 0.01d0
      p%t          = 0.d0
      p%dt         = 1.d0
   end function par_wci

   function par_disp() result(p)
   type(parameters) :: p
      p%wci         = 0
      p%single_dir  = 0
      p%shoaldelay  = 0
      p%px          = 4.d0*atan(1.d0)
      p%Trep        = 10.d0
      p%g           = 9.81d0
      p%eps         = 0.01d0
      p%t           = 0.d0
      p%dt          = 0.5d0
   end function par_disp

   subroutine init_state
   ! (uses host-associated s, t_nx, t_ny, t_nth, t_nts)
      s%nx       => t_nx
      s%ny       => t_ny
      s%ntheta   => t_nth
      s%ntheta_s => t_nts

      allocate(wete(t_nx+1,t_ny+1));        s%wete => wete
      allocate(cg(t_nx+1,t_ny+1));          s%cg   => cg
      allocate(c(t_nx+1,t_ny+1));           s%c    => c
      allocate(sigm(t_nx+1,t_ny+1));        s%sigm => sigm
      allocate(hhw(t_nx+1,t_ny+1));         s%hhw  => hhw
      allocate(u(t_nx+1,t_ny+1));           s%u    => u
      allocate(v(t_nx+1,t_ny+1));           s%v    => v
      allocate(k(t_nx+1,t_ny+1));           s%k    => k
      allocate(nn(t_nx+1,t_ny+1));          s%n    => nn
      allocate(l1(t_nx+1,t_ny+1));          s%L1   => l1
      allocate(thetamean(t_nx+1,t_ny+1));   s%thetamean => thetamean
      allocate(alfaz(t_nx+1,t_ny+1));       s%alfaz => alfaz
      allocate(wm(t_nx+1,t_ny+1));          s%wm   => wm
      allocate(hbar(t_nx+1,t_ny+1));        s%H    => hbar
      allocate(costh(t_nx+1,t_ny+1,t_nth)); s%costh => costh
      allocate(sinth(t_nx+1,t_ny+1,t_nth)); s%sinth => sinth
      allocate(sigt(t_nx+1,t_ny+1,t_nth));  s%sigt  => sigt
      allocate(cgx(t_nx+1,t_ny+1,t_nth));   s%cgx   => cgx
      allocate(cgy(t_nx+1,t_ny+1,t_nth));   s%cgy   => cgy
      allocate(ctheta(t_nx+1,t_ny+1,t_nth));s%ctheta=> ctheta
      allocate(cx(t_nx+1,t_ny+1,t_nth));    s%cx    => cx
      allocate(cy(t_nx+1,t_ny+1,t_nth));    s%cy    => cy
      u = 0.d0; v = 0.d0; hbar = 0.d0

   end subroutine init_state

   function ones2d() result(r)
   real*8 :: r(t_nx+1,t_ny+1)
      r = 1.d0
   end function ones2d

   subroutine cmp3(label, a, b)
   character(len=*), intent(in) :: label
   real*8, intent(in) :: a(:,:,:), b(:,:,:)

      call check(label//" bit-exact", all(bit_equal(a, b)))

   end subroutine cmp3

   subroutine check(label, cond)
   character(len=*), intent(in) :: label
   logical, intent(in) :: cond
   if (cond) then
      write(*,'(A,A)') 'PASS ', label
   else
      write(*,'(A,A)') 'FAIL ', label
      nfail = nfail + 1
   endif
   end subroutine check

   subroutine hdr(label)
   character(len=*), intent(in) :: label
      write(*,'(A,A,A)') '=== ', label, ' ==='
   end subroutine hdr

   function it2(n) result(r)
   integer, intent(in) :: n
   character(len=4) :: r
      write(r,'(I4)') n
      r = adjustl(r)
   end function it2

   function es2(x) result(r)
   real*8, intent(in) :: x
   character(len=24) :: r
      write(r,'(ES24.16)') x
      r = adjustl(r)
   end function es2

   subroutine summary
   if (nfail == 0) then
      write(*,'(A)') 'ALL HOT-KERNEL-2 TESTS PASSED'
   else
      write(*,'(I0,A)') nfail, ' HOT-KERNEL-2 TESTS FAILED'
      stop 1
   endif
   end subroutine summary

end program hotkernels2
