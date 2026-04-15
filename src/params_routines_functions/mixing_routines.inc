
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Andrassy et al. 2023 
! Rees & Izzard 2024
!
! Adjust temperature gradient in the penetration zone to be adiabatic.
! Only works if conv_premix = .true. since the last iteration in a time
! step has NaNs in D_mix if conv_premix = .false.
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine adjust_mlt_gradT_fraction_conv_pen(id, ierr)
    use utils_lib, only: is_bad
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    integer :: k, i, max_mt, idx

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (.not. seed_phase_over) then
        s% do_conv_premix = .false.
        return
    end if

    ! if (is_bad(s% D_mix(1))) return

    ! Apparently required
    s% do_conv_premix = .true.

    idx = -1
    do i = 1, s% nz, 1
        if (idx .gt. 0) exit
        ! Only care about core
        if ((s% lnT(i) .gt. 6.75d0) .and. (s% lnd(i) .gt. 1d0) .and. (s% eps_nuc(i) .gt. 1)) then
            idx = i
        end if
    end do

    ! write(*,*) 'idx:', idx

    if (s% overshoot_scheme(1) == 'step') then
        do k = s%nz, idx, -1
            if (s%D_mix(k) .lt. s% min_D_mix) exit
            ! f is the fraction such that grad_T = f*grad_ad + (1-f)*grad_rad
            s% adjust_mlt_gradT_fraction(k) = 1d0
        end do
    end if

    ! write(*,*) 'fixed gradT'
    ! write(*,*) 'power:', s% power_nuc_burn

end subroutine adjust_mlt_gradT_fraction_conv_pen

! subroutine my_other_adjust_mlt_gradT_fraction_Gumbel(id, ierr)
        ! Michielsen+ 2019 
!     integer, intent(in) :: id
!     integer, intent(out) :: ierr
!     type(star_info), pointer :: s
!     real(dp) :: mass_size_exp_ovr, f      ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
!     integer :: k, loc_core_boundary, loc_exp_boundary
!     real(dp) :: f_above, Hp, r0, D0, fctr

!     ierr = 0
!     call star_ptr(id, s, ierr)
!     if (ierr /= 0) return

!     f_above = s% overshoot_f_above_burn_h_core
!     Hp = s% core_overshoot_Hp
!     r0 = s% core_overshoot_r0
!     D0 = s% x_ctrl(11)

!     do k= s% nz, 1, -1                          ! From center of core to core bounndary
!         if (s% r(k) < r0) then
!             s% adjust_mlt_gradT_fraction(k) = 1
!         else
!             loc_core_boundary = k
!             exit
!         endif
!     end do

!     do k = loc_core_boundary, 1, -1             ! From core bounndary to end exponential
!         fctr = 1-exp(-exp(-abs((s% r(k)-r0)/(f_above* Hp))))    ! calculate if mixing would be higher than min_D_mix
!         if (D0*fctr < s% min_D_mix) then
!             loc_exp_boundary = k
!             exit
!         end if
!     end do

!     mass_size_exp_ovr = DBLE(s% q(loc_exp_boundary) - s% q(loc_core_boundary))
!     do k=loc_core_boundary, loc_exp_boundary, -1
!         f = DBLE(s% q(loc_exp_boundary)- s% q(k))/mass_size_exp_ovr ! zero at k=loc_exp_boundary and one at k= loc_step_boundary
!         s% adjust_mlt_gradT_fraction(k) = f
!     end do

! end subroutine my_other_adjust_mlt_gradT_fraction_Gumbel

! subroutine other_adjust_mlt_gradT_fraction_overshoot(id, ierr)
        ! Chris Lindsay 2023
!         integer, intent(in) :: id
!         integer, intent(out) :: ierr
!         type(star_info), pointer :: s
!         real(dp) :: fraction, f_ov, Hp, m0, dr       ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
!         integer :: i_bdy, k
!         logical, parameter :: DEBUG = .FALSE.

!         ierr = 0
!         call star_ptr(id, s, ierr)
!         if (ierr /= 0) return

!         if (.not. s% x_logical_ctrl(1)) return

!         f_ov = s% x_ctrl(1)

!         m0 = s% cz_bot_mass(1)
        
!         if (m0 <= 0 ) return

!         do k = 1, s% nz
!           if ( s % m(k) < m0 ) then
!             i_bdy = k - 1
!             exit
!           endif
!         end do
        
!         do k= i_bdy, s%nz
!               call evaluate_Hp (s, k, s%r(k), Hp, ierr)
!               dr = s%r(i_bdy) - s%r(k)
!               if ( dr < 0 ) then
!                 fraction = -1._dp
!               else if ( dr < (10 * f_ov * Hp) ) then
!                 fraction = exp((-2 * dr)/(f_ov * Hp))
!               else
!                 fraction = 0._dp
!               end if

!               s% adjust_mlt_gradT_fraction(k) = fraction
!         end do


!       end subroutine other_adjust_mlt_gradT_fraction_overshoot