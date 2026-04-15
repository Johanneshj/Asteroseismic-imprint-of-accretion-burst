subroutine set_varcontrol(id, ierr) 

    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    real(dp) :: varcontrol_min, varcontrol_max
    real(dp) :: minimum_varcontrol_min, minimum_varcontrol_max
    real(dp) :: age_in_secs, elapsed_time
    real(dp) :: when_to_decrease, midpoint, steepness, term_1, term_2, front_factor, duration
    real(dp) :: varcontrol_profile, minimum_varcontrol_profile

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    varcontrol_min = 1d-6
    varcontrol_max = 1d-4

    minimum_varcontrol_min = varcontrol_min / 10d0
    minimum_varcontrol_max = varcontrol_max / 10d0

    age_in_secs = s% star_age * secyer
    elapsed_time = age_in_secs - bp% target_age
    
    ! ###########################################################
    ! These values give a reasonable profile for the varcontrol
    when_to_decrease = (5d0 * amount_of_mass_to_begin_decreasing_timestep_before_burst) / &
                            (bp% mdot_ini * secyer / Msun) * secyer
    midpoint = when_to_decrease / 2d0
    steepness = 0.07d0 * midpoint
    ! write(*,*) 'varcontrol stuff:', when_to_decrease, midpoint, steepness, elapsed_time, age_in_secs, bp% target_age
    ! ###########################################################

    duration = bp% burst_rise_time + &
                bp% burst_plateau_time + &
                bp% burst_decay_time

    front_factor = (varcontrol_max - varcontrol_min)
    term_1 = 1d0 / (1d0 + EXP((midpoint - elapsed_time) / steepness)) 
    term_2 = 1d0 / (1d0 + EXP(((midpoint + duration) + elapsed_time) / steepness))
    
    varcontrol_profile = varcontrol_min + front_factor * (term_1 + term_2)
    minimum_varcontrol_profile = varcontrol_profile / 10d0

    s% varcontrol_target = varcontrol_profile
    s% min_allowed_varcontrol_target = minimum_varcontrol_profile

end subroutine set_varcontrol

subroutine set_lgL_limit(id, ierr)
    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    real(dp) :: lgL_limit_min, lgL_limit_max
    real(dp) :: age_in_secs, elapsed_time
    real(dp) :: when_to_decrease, midpoint, steepness, term_1, term_2, front_factor, duration
    real(dp) :: lgL_limit_profile
    real(dp) :: M_minus_Mcc, lg_L_limit_vs_Mcc

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    if (.not. seed_phase_over) then
        s% delta_lgL_limit = 1d-2
        s% delta_lgL_hard_limit = -1
        return
    end if

    lgL_limit_min = 1d-4
    lgL_limit_max = 1d-2

    ! ########## Sigmoid profile ########## !
    age_in_secs = s% star_age * secyer
    elapsed_time = age_in_secs - bp% target_age
    ! ###########################################################
    ! These values give a reasonable profile for the delta_lgL_limit profile
    when_to_decrease = (5d0 * amount_of_mass_to_begin_decreasing_timestep_before_burst) / &
                            (bp% mdot_ini * secyer / Msun) * secyer
    midpoint = when_to_decrease / 2d0
    steepness = 0.07d0 * midpoint
    ! ###########################################################
    duration = bp% burst_rise_time + &
                bp% burst_plateau_time + &
                bp% burst_decay_time
    ! Sigmoid decrease/increase around burst
    front_factor = (lgL_limit_max - lgL_limit_min)
    term_1 = 1d0 / (1d0 + EXP((midpoint - elapsed_time) / steepness)) 
    term_2 = 1d0 / (1d0 + EXP(((midpoint + duration) + elapsed_time) / steepness))
    lgL_limit_profile = lgL_limit_min + front_factor * (term_1 + term_2)
    ! ##################################### !


    ! ! ########## Limit lgL based on convective core mass ########## !
    ! M_minus_Mcc = s% mass_conv_core / s% star_mass
    ! ! Sigmoid profile based on ratio between star mass and convective core mass
    ! ! If convective core grows (in mass) decrease delta_lgL_limit
    ! ! This works because, in general, the convective core grows during the pre-MS
    ! lg_L_limit_vs_Mcc = lgL_limit_min + (lgL_limit_max - lgL_limit_min) * &
    !                     1d0 / (1d0 + EXP(40d0*(M_minus_Mcc-0.3d0)))
    ! ! 0.3d0 means that the transition between max and min lgL_limit
    ! ! happens when Mcc is 20% of the total mass, with steepness 40 (smooth)
    ! ! ############################################################# !

    ! Failsafe; limit either based on convective core or burst
    s% delta_lgL_limit = lgL_limit_profile ! MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)
    s% delta_lgL_H_limit = lgL_limit_profile ! MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)
    s% delta_HR_limit = lgL_limit_profile ! MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)

    s% delta_lgL_hard_limit = 10d0 * MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)
    s% delta_lgL_H_hard_limit = 10d0 * MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)
    s% delta_HR_hard_limit = 10d0 * MIN(lgL_limit_profile, lg_L_limit_vs_Mcc)

end subroutine set_lgL_limit






