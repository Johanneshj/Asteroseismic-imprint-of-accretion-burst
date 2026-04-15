
subroutine save_profiles_during_burst(id, ierr, extras_check_model)
    use star_def
    use const_def
    
    integer, intent(in) :: id
    integer :: ierr

    real(dp) :: dt_save_in_years, dt_save_in_years_during_burst, dt_save_in_years_after_burst
    real(dp) :: dt_diff, dt_tolerance

    type (star_info), pointer :: s
    

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    extras_check_model = keep_going  

    if (s% x_integer_ctrl(10) .ne. 1) return

    if (old_age .lt. 0) then
        s% need_to_save_profiles_now = .true.
        old_age = s% star_age
    end if

    dt_save_in_years_during_burst = 10
    dt_save_in_years_after_burst = 50
    dt_tolerance = 0.01 

    if (.not. bp% burst_over) then ! save profiles during burst
        dt_diff = (s% star_age - old_age)
        if (dt_diff .gt. dt_save_in_years_during_burst) then
            if ((dt_diff - dt_save_in_years_during_burst) .gt. dt_tolerance) then
                s% dt = s% dt * 0.9d0
                write(*,*) 'overshoot age for profile save -> retry'
                extras_check_model = retry
            else
                s% need_to_save_profiles_now = .true.
                s% need_to_update_history_now = .true.
                s% save_profiles_model_priority = 4
                old_age = s% star_age
            end if
        end if
    else if (bp% burst_over) then ! save profiles after burst, that is, after mdot has returned to mdot_ini
        dt_diff = (s% star_age - old_age)
        if (dt_diff .gt. dt_save_in_years_after_burst) then
            if ((dt_diff - dt_save_in_years_after_burst) .gt. dt_tolerance) then
                s% dt = s% dt * 0.9d0
                write(*,*) 'overshoot age for profile save -> retry'
                extras_check_model = retry
            else
                s% need_to_save_profiles_now = .true.
                s% need_to_update_history_now = .true.
                s% save_profiles_model_priority = 5
                old_age = s% star_age
            end if
        end if
    end if

end subroutine save_profiles_during_burst