subroutine mdot_single_burst(id, ierr)
    use star_def
    use const_def
    implicit none
    integer, intent(in) :: id
    integer, intent(out) :: ierr
    integer :: step, i

    real(dp) :: mdot_beg, burst_profile, decay_profile

   ! type(burst_parameters_single), intent(inout) :: bp
    type(star_info), pointer :: s

    ! write(*,*) 'mdot_single_burst'

    ierr = 0
    call star_ptr(id, s, ierr) ! retrieve star id

    ! call load_burst_parameters(id, ierr, bp)
    call key_points(id, ierr)
    ! write(*,*) 'mdot_ini:', bp% mdot_ini, bp% mdot_ini * secyer / Msun
    ! integrate sigmoid increase
    if (bp% sigmoid_mass .lt. 0) then
        bp% sigmoid_mass = integrate_sigmoid_mass(bp% sigmoid_mass, bp% mdot_ini) 
    end if

    if (bp% number_of_bursts .gt. 0) then
        if (bp% current_burst .eq. 1) then
            if (bp% target_age .lt. 0) then
                write(*,*) 'calculate target_age'
                bp% target_age = ((bp% burst_mass_start - bp% sigmoid_mass)/(bp% mdot_ini)) + & 
                                    (years_for_initial_sigmoid_increase * secyer - bp% burst_rise_time)
            end if
        else
            write(*,*) 'calculate target_age post burst'
            bp% target_age = ((bp% final_mass - bp% mass_end_burst)/(bp% mdot_ini)) + & 
                                        (bp% t_end_burst)
        end if
    else
        ! If we just want a constant accretion model with no bursts
        if (bp% target_age .lt. 0) then
            bp% target_age = ((bp% final_mass - bp% sigmoid_mass)/(bp% mdot_ini)) + 4000 * secyer
            write(*,*) 'end of accretion target age:', bp% target_age / secyer
        end if
    end if

    if (s% star_mass * Msun .lt. bp% final_mass) then
        if (bp% number_of_bursts .gt. 0) then
            if ((initial_increase .eqv. .true.) .and. (bp% no_more_bursts .eqv. .false.)) then
                mdot_beg = calculate_initial_increase(s% star_age * secyer, bp% mdot_ini)
                s% mstar_dot = MIN(mdot_beg, bp% mdot_ini)
                write(*,*) '#### 1 ####'
            else if ((initial_increase .eqv. .false.) .and. (bp% no_more_bursts .eqv. .false.)) then
                if (s% x_integer_ctrl(10) .ne. 0) then
                    burst_profile = calculate_burst_profile(id, ierr, bp) ! burst
                    s% mstar_dot = MAX(bp% mdot_ini, burst_profile)
                    write(*,*) '#### 2 ####'
                else
                    s% mstar_dot = bp% mdot_ini ! in saving pre-burst model or constant accretion scenario
                    write(*,*) '#### 3 ####'
                end if
            else if ((initial_increase .eqv. .false.) .and. (bp% no_more_bursts .eqv. .true.)) then
                decay_profile = calculate_decay_profile(id, ierr, bp)
                s% mstar_dot = MIN(bp% mdot_ini, decay_profile)
                write(*,*) '#### 4 ####'
            end if
        else
            bp% no_more_bursts = .true.
            if ((initial_increase .eqv. .true.) .and. (bp% no_more_bursts .eqv. .true.)) then
                mdot_beg = calculate_initial_increase(s% star_age * secyer, bp% mdot_ini)
                s% mstar_dot = MIN(mdot_beg, bp% mdot_ini)
                write(*,*) '#### 5 ####'
            else if ((initial_increase .eqv. .false.) .and. (bp% no_more_bursts .eqv. .true.)) then
                decay_profile = calculate_decay_profile(id, ierr, bp)
                s% mstar_dot = MIN(bp% mdot_ini, decay_profile)
                write(*,*) '#### 6 ####'
            end if
        end if
    else
        if (final_mass_reached_within_tolerance) then
            s% mstar_dot = 0.0d0
            write(*,*) '#### 7 ####'
        else
            write(*,*) 'still need to reach final mass within tolerance'
            s% mstar_dot = MIN(bp% mdot_ini, decay_profile)
        end if
    end if   

end subroutine mdot_single_burst


subroutine load_burst_parameters(id, ierr, bp)
    use star_def
    use const_def

    implicit none

    integer, intent(in) :: id
    integer, intent(out) :: ierr

    type(burst_parameters_single), intent(inout) :: bp
    type(star_info), pointer :: s
    
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    bp% burst_over          = .false.
    bp% no_more_bursts      = .false.
    bp% T                   = s% x_integer_ctrl(5)
    bp% final_mass          = s% x_ctrl(3) * Msun
    bp% mdot_ini            = s% x_ctrl(15) * Msun / secyer
    bp% sigmoid_mass        = -1
    bp% target_age          = -1
    bp% number_of_bursts    = s% x_integer_ctrl(1)
    bp% current_burst       = 1
    bp% t_end_burst         = -1
    bp% mass_end_burst      = -1
    bp% t_end_accretion     = -1

    bp% burst_rise_time     = s% x_ctrl(10) * secyer
    bp% burst_plateau_time  = s% x_ctrl(11) * secyer
    bp% burst_decay_time    = s% x_ctrl(12) * secyer
    bp% burst_mdot          = s% x_ctrl(13) * Msun / secyer
    bp% burst_mass_start    = s% x_ctrl(14) * Msun

end subroutine load_burst_parameters


function integrate_sigmoid_mass(sigmoid_mass, mdot_ini) result(sigmoid_m)
    use const_def, only: Msun, secyer
    real(dp) :: dt
    real(dp), allocatable :: ages(:), sigmoid_vals(:)
    integer :: nsteps, i
    real(dp) :: sigmoid_mass, sigmoid_m, mdot_init, mdot_ini 

    sigmoid_init = 1d-12
    years_for_initial_sigmoid_increase = 4000d0
    mdot_init = mdot_ini * secyer / Msun
    
    if (sigmoid_mass .lt. 0.0d0) then
        nsteps = 40000
        allocate(ages(nsteps), sigmoid_vals(nsteps))
        dt = (years_for_initial_sigmoid_increase - 0d0) / real(nsteps, dp)
        do i = 1, nsteps
            ages(i) = 0d0 + dt * real(i - 1, dp)
        end do

        do i = 1, nsteps
            sigmoid_vals(i) = sigmoid_init + (mdot_init - sigmoid_init) / &
                                 (1d0 + EXP( - (ages(i) - years_for_initial_sigmoid_increase/2d0)/(100d0) ))
        end do

        sigmoid_mass = 0.01d0
        do i = 1, nsteps-1 
            sigmoid_mass = sigmoid_mass + (sigmoid_vals(i) + sigmoid_vals(i+1)) * (ages(i+1) - ages(i)) / (2d0)
        end do
    end if
    sigmoid_m = sigmoid_mass * Msun
    write(*,*) '-----> Sigmoid mass:', sigmoid_m
end function integrate_sigmoid_mass


function calculate_burst_profile(id, ierr, bp) result(burst_profile)
    use star_def
    use const_def

    implicit none

    integer, intent(in) :: id
    integer, intent(out) :: ierr

    type(burst_parameters_single), intent(in) :: bp
    type(star_info), pointer :: s

    real(dp) :: burst_profile ! output
    real(dp) :: mdot_ini, mdot_burst
    real(dp) :: exp_rise, exp_decay
    real(dp) :: term_rise, term_decay, front_factor
    real(dp) :: rise_time, plateau_time, decay_time
    real(dp) :: age_in_secs 
    real(dp) :: elapsed_time

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    rise_time = bp% burst_rise_time
    plateau_time = bp% burst_plateau_time
    decay_time = bp% burst_decay_time
    mdot_burst = bp% burst_mdot
    mdot_ini = bp% mdot_ini

    age_in_secs = s% star_age * secyer
    elapsed_time = age_in_secs - bp% target_age

    front_factor = (mdot_burst-mdot_ini)
    exp_rise = EXP(-10d0 * (elapsed_time - rise_time/2d0)/rise_time)
    term_rise = 1d0/(1d0 + exp_rise)
    exp_decay = EXP(-10d0 * (elapsed_time - (rise_time + plateau_time + decay_time/2d0))/decay_time)
    term_decay = 1d0/(1d0 + exp_decay)

    burst_profile = mdot_ini + front_factor * &
                        ( &
                            term_rise - &
                            term_decay &
                        )      

end function calculate_burst_profile

function calculate_decay_profile(id, ierr, bp) result(decay_profile)
    use star_def
    use const_def
    ! use extras_check_model

    implicit none

    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s
    type(burst_parameters_single), intent(in) :: bp

    real(dp) :: target_age
    real(dp) :: decay_profile
    real(dp) :: mdot_ini
    real(dp) :: tau
    real(dp) :: elapsed_time
    real(dp) :: mini, maxi
    real(dp) :: mass_to_be_accreted, duration_of_decay

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    ! extras_check_model = keep_going

    mdot_ini = bp% mdot_ini
    mass_to_be_accreted = s% x_ctrl(4) * Msun
    if (target_age_for_decay .lt. 0) then
        target_age_for_decay = bp% target_age - (mass_to_be_accreted / bp% mdot_ini)
    end if
    ! 5e5 means we accrete the last 0.05 solar masses over 500.000 years.
    duration_of_decay = s% x_ctrl(16) ! standard value, we could do 10 * (target_age_for_decay - bp% target_age), but let's keep it simple for now
    tau = s% x_ctrl(5) * secyer ! characteristic decay scale to ensure we cover "mass_to_be_accreted" in "duration_of_decay" and cover mdot max and min
    ! If you want other decay profiles, use python script on to solve for tau that fits)
    elapsed_time = s% star_age * secyer - target_age_for_decay

    mini = 1d-9 * Msun / secyer ! min mdot
    maxi = mdot_ini ! max mdot
    
    ! if ((s% star_age * secyer .gt. target_age_for_decay + 10 * secyer) .and. &
    !         (.not. decay_at_end_of_accretion_can_begin)) then
    !     s% dt = s% dt * 0.9d0
    !     write(*,*) 'overshoot age for final decay -> retry'
    !     extras_check_model = retry
    ! else if ((s% star_age * secyer .lt. target_age_for_decay + 10 * secyer) &
    !         (s% star_age * secyer .gt. target_age_for_decay) &
    !         (.not. decay_at_end_of_accretion_can_begin)) then
    !     decay_at_end_of_accretion_can_begin = .true.
    ! end if
    
    if (decay_at_end_of_accretion_can_begin) then
        decay_profile = mini + (maxi - mini) * EXP(-(elapsed_time/tau))
    else
        decay_profile = maxi
    end if

    ! write(*,*) 'decay:', bp% target_age/secyer, mass_to_be_accreted / Msun, elapsed_time / secyer, decay_profile
end function calculate_decay_profile

function calculate_initial_increase(age_in_secs, mdot_ini) result(mdot_beg)
    real(dp) :: mdot_ini, age_in_secs, mdot_beg, init
    ! 1000 years sigmoid increase in accretion rate to "ease" seed into accretion
    init = sigmoid_init * Msun / secyer
    mdot_beg = init + (mdot_ini - init) / &
            ( &
                (1d0 + EXP(-(age_in_secs/secyer - years_for_initial_sigmoid_increase/2d0) / (100d0))) & 
                ! sigmoid decay centered at rise_time/2 with e-folding time 100 years
            )
end function calculate_initial_increase

subroutine key_points(id, ierr)
    use star_def
    use const_def

    implicit none

    integer, intent(in) :: id
    integer, intent(out) :: ierr
    type(star_info), pointer :: s

    real(dp) :: age_in_secs
    real(dp) :: target_age
    real(dp) :: burst_duration
    real(dp) :: number_of_bursts
    real(dp) :: burst_profile

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! ########### Saving key moments during the evolution ########### !

    if (bp% target_age .lt. 0) return ! failsafe so restarts are possible

    age_in_secs = s% star_age * secyer
    target_age = bp% target_age 
    burst_duration = bp% burst_rise_time + bp% burst_plateau_time + bp% burst_decay_time
    number_of_bursts = bp% number_of_bursts

    if (age_in_secs .le. 4000*secyer) then
        initial_increase = .true.
    else
        if (initial_increase) then
            write(*,*) '###### initial increase is over ######'
        end if
        initial_increase = .false.
    end if

    if (.not. seed_phase_over) then
        if (center_c12_ini .lt. 0) then
            center_c12_ini = s% center_c12
        end if
        if (s% center_c12 .lt. central_c12_decrease_to_mark_end_of_seed_phase * center_c12_ini) then
            seed_phase_over = .true.
            write(*,*) '###### Seed phase (marked as no c12 depletion) is over ######'
        end if
    end if

    ! burst_profile = calculate_burst_profile(id, ierr, bp)
    ! write(*,*) 'lol:', burst_profile * secyer / Msun, ((bp% mdot_ini * secyer / Msun + 0.001d-6) )
    if ((age_in_secs .gt. target_age + burst_duration)) then ! are we past burst duration?
        burst_profile = calculate_burst_profile(id, ierr, bp)
        if (burst_profile * secyer / Msun .le. (bp% mdot_ini * secyer / Msun + 0.0005d-6)) then ! have we returned to mdot init within tolerance?
            write(*,*) '###### mdot returned to mdot init after burst within tolerance ######'
            if (s% x_integer_ctrl(10) .gt. 0) then ! failsafe, not clean
                ! if x_integer_ctrl(10) = 2, current burst will become 2 and t_end_burst will be greater than 0.
                ! Then next "target age" will become end of accretion phase.
                bp% burst_over = .true.
                if (bp% t_end_burst .lt. 0) then
                    bp% t_end_burst = age_in_secs
                    bp% mass_end_burst = s% star_mass * Msun
                    bp% current_burst = bp% current_burst + 1
                    write(*,*) '###### A burst is completed ######'
                end if
            end if
        end if
    else
        if (s% x_integer_ctrl(10) .ne. 1) then
            bp% burst_over = .false.
        end if
    end if

    if ((bp% current_burst .gt. number_of_bursts) .and. (.not. bp% no_more_bursts)) then
        bp% no_more_bursts = .true.
        write(*,*) '###### Final burst over ######'
    end if

    if (s% x_integer_ctrl(10) .eq. 2) then
        bp% no_more_bursts = .true.
        write(*,*) '###### Constant accretion scenario --> no bursts ######'
    end if
    

end subroutine key_points