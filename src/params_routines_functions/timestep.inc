subroutine determine_maximum_timestep(id, ierr)

    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    ! #### Determines maximum timestep in seconds ####
    !   Calculates thermal timescale in seconds.
    !   Calculates accretion timescale in seconds.
    !   Calculates sigmoid-shaped dt-profile, to "ease" model into burst
    ! ################################################

    integer  :: k
    real(dp) :: mass
    real(dp) :: mdot
    real(dp) :: luminosity
    real(dp) :: radius
    real(dp) :: timestep_factor
    real(dp) :: age_in_secs
    real(dp) :: duration
    real(dp) :: acc_timescale, therm_timescale, max_dt_in_years
    real(dp) :: maximum_timestep

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    mass = s% star_mass * Msun
    mdot = s% mstar_dot
    luminosity = s% L(1)
    radius = s% r(1)
    timestep_factor = s% x_ctrl(6)
    age_in_secs = s% star_age * secyer

    ! Thermal time scale
    therm_timescale = thermal_timescale(timestep_factor, mass, radius, luminosity)

    ! Accretion time scale
    acc_timescale = accretion_timescale(timestep_factor, mass, mdot)

    ! Maximum timestep increase post burst
    max_dt_in_years = set_max_timestep_after_burst(id, ierr, timestep_factor, age_in_secs)

    ! Determine which one dominates
    maximum_timestep = min(therm_timescale, &
                            acc_timescale, &
                            max_dt_in_years*secyer)

    ! Set maximum timestep in seconds
    s% max_timestep = maximum_timestep
    
    write(*,*) "max timestep in years:", s% max_timestep / secyer
    ! if (s% dt .gt. s% max_timestep) then
    !     write(*,*) 'timestep sizes:', s% dt/secyer, s% max_timestep/secyer
    !     s% dt = (s% dt) / 2
    !     retry_flag = .true.
    ! else
    !     retry_flag = .false.
    ! end if

end subroutine determine_maximum_timestep


function thermal_timescale(timestep_factor, mass, radius, luminosity) result(therm_timescale)
    implicit none
    ! mass, cgrav, rad and lum in cgs units
    real(dp) :: timestep_factor, mass, radius, luminosity
    real(dp) :: therm_timescale
    therm_timescale = timestep_factor * (mass)**2 * standard_cgrav / (2 * radius * luminosity)
end function thermal_timescale


function accretion_timescale(timestep_factor, mass, mdot) result(acc_timescale)
    implicit none
    ! mass and mdot in cgs units
    real(dp) :: timestep_factor, mass, mdot
    real(dp) :: acc_timescale
    acc_timescale = timestep_factor * mass / mdot 
end function accretion_timescale

function set_max_timestep_after_burst(id, ierr, timestep_factor, age_in_secs) result(max_dt_in_years)
    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    integer  :: k
    real(dp) :: timestep_factor
    real(dp) :: max_dt_in_years
    real(dp) :: age_in_secs
    

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (s% x_integer_ctrl(10) .ne. 1) return

    ! Set maximum timestep if bp% t_end_accretion is not calculated yet, i.e., burst hasen't ended
    write(*,*) 'end accretion:', bp% t_end_burst / secyer
    write(*,*) 1 / secyer

    if (bp% t_end_burst .lt. 0) then
        max_dt_in_years = 0.04
        write(*,*) "max dt in years:", max_dt_in_years, timestep_factor
        return
    end if

    ! Set maximum timestep if burst has ended
    if (age_in_secs .gt. bp% t_end_burst) then
        max_dt_in_years = 1d1 ! 10 years
    else
        max_dt_in_years = 0.04! two weeks
    end if
    write(*,*) "max dt in years:", max_dt_in_years, timestep_factor
    

    ! end_age = bp% target_age + bp% burst_rise_time + bp% burst_plateau_time + bp% burst_decay_time
    
    ! if ((age_in_secs .gt. end_age) .and. (s% mstar_dot .le. bp% mdot_ini)) then
    !     ! If we are beyond burst and mdot has returned to mdot_ini we can increase the timestep again
    !     max_dt_in_years = 1d1 ! 10 years
    ! else
    !     max_dt_in_years = 0.04! two weeks
    ! end if

    

end function set_max_timestep_after_burst

subroutine set_dX_div_X_after_seed_phase(id, ierr)
    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    integer  :: k

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (s% model_number .lt. 100) return

    if (s% x_integer_ctrl(10) .ne. 0) return

    s% dX_limit_species(1) = 'he3'
    s% dX_limit_species(2) = 'c12'
    s% dX_limit_species(3) = 'c13'
    s% dX_limit_species(4) = 'n14'
    s% dX_limit_species(5) = 'n15'
    s% dX_limit_species(6) = 'o16'
    s% dX_limit_species(7) = 'o17'
    s% dX_limit_species(8) = 'o18'
    s% dX_limit_species(9) = 'f19'
    s% dX_div_X_limit(1:9) = 1d-3
    s% dX_div_X_hard_limit(1:9) = -1
    s% dX_div_X_limit_min_X(1:9) = 1d-12
    s% dX_decreases_only(1:9) = .false.
    s% dX_div_X_at_high_T_limit(1:9) = 1d-3
    s% dX_div_X_at_high_T_hard_limit(1:9) = -1
    s% dX_div_X_at_high_T_limit_lgT_min(1:9) = 6

end subroutine set_dX_div_X_after_seed_phase
