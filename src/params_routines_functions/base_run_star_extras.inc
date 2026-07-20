
         
         subroutine extras_controls(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            real(dp) :: H, He, r_h2, r_he3 ! definitions
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            ! this is the place to set any procedure pointers you want to change
            ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
   
            ! the extras functions in this file will not be called
            ! unless you set their function pointers as done below.
            ! otherwise we use a null_ version which does nothing (except warn).
   
            s% extras_startup => extras_startup
            s% extras_start_step => extras_start_step
            s% extras_check_model => extras_check_model
            s% extras_finish_step => extras_finish_step
            s% extras_after_evolve => extras_after_evolve
            s% how_many_extra_history_columns => how_many_extra_history_columns
            s% data_for_extra_history_columns => data_for_extra_history_columns
            s% how_many_extra_profile_columns => how_many_extra_profile_columns
            s% data_for_extra_profile_columns => data_for_extra_profile_columns  
            
            s% other_energy => energy_routine
            s% other_adjust_mdot => other_adjust_mdot
   
            s% how_many_extra_history_header_items => how_many_extra_history_header_items
            s% data_for_extra_history_header_items => data_for_extra_history_header_items
            s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
            s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

            ! He4 and He3 abundances
            He = 0.248d0 + (1.4d0 * s% initial_z)
            r_he3 = s% x_ctrl(2) * 1d-6 

            s% job% initial_he4 = He / (1d0 + r_he3)
            s% job% initial_he3 = s% job% initial_he4 * r_he3
            s% accretion_he4 = s% job% initial_he4 
            s% accretion_he3 = s% job% initial_he3
            
            ! Hydrogen and deuterium abundances
            H = 1d0 - He - s% initial_z
            r_h2 = s% x_ctrl(1) * 1d-6

            s% job% initial_h1 = H / (1d0 + r_h2)
            s% job% initial_h2 = s% job% initial_h1 * r_h2
            s% accretion_h1 = s% job% initial_h1
            s% accretion_h2 = s% job% initial_h2

            s% kap_rq% Zbase = s% initial_z
         
            s% how_many_other_mesh_fcns => how_many_mesh_fcns
            s% other_mesh_fcn_data => gradr_grada_mesh_fcn_data
            ! s% other_timestep_limit => limit_convective_core_changes
            ! s% other_timestep_limit => limit_gradr_changes
            ! s% other_adjust_mlt_gradT_fraction => adjust_mlt_gradT_fraction_conv_pen

         end subroutine extras_controls
         
         subroutine other_adjust_mdot(id, ierr)
            use star_def
            implicit none
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            integer :: step, i
            type(star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)

            s% mstar_dot = s% x_ctrl(15) * Msun / secyer

            call mdot_single_burst(id, ierr)

            ! write(*,*) 'initial_increase', initial_increase
            write(*,*) 'accretion rate:', s% mstar_dot * secyer / Msun

         end subroutine other_adjust_mdot
   
         subroutine extras_startup(id, restart, ierr)
            integer, intent(in) :: id
            logical, intent(in) :: restart
            integer, intent(out) :: ierr

            real(dp) :: r
            integer :: i

            type (star_info), pointer :: s
            real(kind=8) :: offset
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            if (s% x_logical_ctrl(2)) then
               above_1_5 = .false.
            else
               call load_burst_parameters(id, ierr, bp)
            end if

         end subroutine extras_startup

         integer function extras_start_step(id)
            integer, intent(in) :: id
            integer :: ierr

            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_start_step = 0 

         end function extras_start_step
   

         integer function extras_check_model(id)
            integer, intent(in) :: id
            integer :: ierr

            real(dp) :: log_Tsurf, log_Lsurf
            real(dp) :: dteff, dlogL, dmass
            real(dp) :: M_final
            real(dp) :: dt_save_in_years, dt_save_in_years_during_burst, dt_save_in_years_after_burst
            real(dp) :: dt_diff, dt_tolerance, mass_to_be_accreted
         
            integer :: f

            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_check_model = keep_going  

            ! %%%%%% initial model %%%%%%% !
            if (s% x_logical_ctrl(2)) then
               if ((s% r(1)/Rsun > s%x_ctrl(3)) .and. (.not. above_1_5)) then
                  above_1_5 = .true.
               end if

               if (above_1_5) then
                  write(*,*) s% r(1)/Rsun,s% x_ctrl(3)
                  if ((s% r(1)/Rsun < s% x_ctrl(3) + 0.0001) .and. (s% r(1)/Rsun >= s% x_ctrl(3))) then
                        extras_check_model = terminate    
                  else if (s% r(1)/Rsun < s% x_ctrl(3) + 0.0001)  then
                        extras_check_model = retry
                  end if
               end if
               if (extras_check_model == terminate) s% termination_code = t_extras_check_model
               return
            end if
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

            M_final = bp% final_mass
            
            call determine_maximum_timestep(id, ierr)

            if (s% x_logical_ctrl(20) .eqv. .true.) call control_mesh(id, ierr)


            write(*,*) s% max_timestep / secyer

            if ((s% x_integer_ctrl(10) .lt. 2) .and. (termination_age .lt. 0)) then ! set to 2 in no burst scenario, then model will terminate at 0.5 hydrogen abundance
               if (s% x_integer_ctrl(10) .eq. 0) then 
                  termination_age = bp% target_age / secyer - 10 ! save pre-burst model
               else if (s% x_integer_ctrl(10) .eq. 1) then
                  termination_age = bp% target_age / secyer + &
                     (bp% burst_rise_time + bp% burst_plateau_time + bp% burst_decay_time) / secyer & 
                        + 1000 ! save post-burst model
               end if
            end if
                  
            if ((s% star_age .gt. termination_age) .and. (.not. reached_target_age)) then
               reached_target_age = .true.
            end if

            if (s% x_integer_ctrl(10) .ne. 2) then
               if (reached_target_age) then
                  if ((s% star_age - termination_age) .gt. 0.01d0) then ! 
                     s% mstar_dot = bp% mdot_ini ! failsafe
                     s% dt = s% dt * 0.9d0
                     write(*,*) 'overshoot termination age -> retry'
                     extras_check_model = retry
                  else
                     write(*,*) 'within tolerance of termination age -> terminating'
                     extras_check_model = terminate
                  end if
               end if
            end if

            if (s% x_integer_ctrl(10) .eq. 1) then
               seed_phase_over = .true. ! failsafe, not clean
            end if

            write(*,*) 'target_age:', bp% target_age / secyer
            mass_to_be_accreted = s% x_ctrl(4) * Msun
            target_age_for_decay = bp% target_age - (mass_to_be_accreted / bp% mdot_ini)
            if (s% x_integer_ctrl(10) .eq. 2) then
               if ((s% star_age * secyer .gt. target_age_for_decay + 0.01 * secyer) .and. &
                        (.not. decay_at_end_of_accretion_can_begin)) then
                  s% dt = s% dt * 0.9d0
                  write(*,*) 'overshoot age for final decay -> retry'
                  extras_check_model = retry
               else if ((s% star_age * secyer .lt. target_age_for_decay + 0.01 * secyer) .and. &
                        (s% star_age * secyer .gt. target_age_for_decay) .and. &
                        (.not. decay_at_end_of_accretion_can_begin)) then
                  decay_at_end_of_accretion_can_begin = .true.
                  write(*,*) 'can decay begin?', decay_at_end_of_accretion_can_begin
               end if
            end if

            write(*,*) s% star_mass, bp% final_mass / Msun, bp% final_mass / Msun * (1 + 1d-9)
            write(*,*) final_mass_reached_within_tolerance


            ! #### CHECK IF FINAL MASS IS OVERSHOT #################
            ! We check if we overshoot final mass within tolerance.
            ! Important if we want to map the full decay.
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            ! Follwing are defined in logics.inc
            !  - decay_at_end_of_accretion_can_begin
            !  - final_mass_reached_within_tolerance
            ! Following are defined in burst_params.inc
            !  - bp% final_mass
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (decay_at_end_of_accretion_can_begin) then
               if (.not. final_mass_reached_within_tolerance) then
                  if (s% star_mass .gt. bp% final_mass / Msun) then
                     if (s% star_mass .gt. bp% final_mass / Msun * (1 + 1d-6)) then
                        write(*,*) 'overshoot final mass within tolerance -> retry'
                        s% dt = s% dt * 0.9d0
                        extras_check_model = retry
                     else if (s% star_mass .lt. bp% final_mass / Msun * (1 + 1d-6)) then
                        write(*,*) 'we are within tolerance'
                        final_mass_reached_within_tolerance = .true.
                     end if
                  end if
               end if
            end if
            ! ######################################################


            ! ##### Retries to save profiles during burst ##### !
            ! if (s% x_integer_ctrl(10) .eq. 1) then
            !    can_save_profile_during_burst = .false.
               
            !    dt_save_in_years_during_burst = 10
            !    dt_save_in_years_after_burst = 50
            !    dt_tolerance = 0.01 

            !    if ((old_age .lt. 0) .and. (s% model_number .eq. 100)) then
            !       can_save_profile_during_burst = .true.
            !       old_age = bp% target_age / secyer - 1
            !    end if

            !    write(*,*) 'old age:', old_age
               
            !    if (old_age .gt. 0) then
            !       if (.not. bp% burst_over) then ! save profiles during burst
            !          dt_diff = (s% star_age - old_age)
            !          if (dt_diff .gt. dt_save_in_years_during_burst) then
            !             if ((dt_diff - dt_save_in_years_during_burst) .gt. dt_tolerance) then
            !                s% dt = s% dt * 0.9d0
            !                write(*,*) 'overshoot age for profile save -> retry'
            !                extras_check_model = retry
            !             else
            !                ! s% need_to_save_profiles_now = .true.
            !                ! s% need_to_update_history_now = .true.
            !                ! s% save_profiles_model_priority = 4
            !                can_save_profile_during_burst = .true.
            !                old_age = s% star_age
            !             end if
            !          end if
            !       else if (bp% burst_over) then ! save profiles after burst, that is, after mdot has returned to mdot_ini
            !          dt_diff = (s% star_age - old_age)
            !          if (dt_diff .gt. dt_save_in_years_after_burst) then
            !             if ((dt_diff - dt_save_in_years_after_burst) .gt. dt_tolerance) then
            !                s% dt = s% dt * 0.9d0
            !                write(*,*) 'overshoot age for profile save -> retry'
            !                extras_check_model = retry
            !             else
            !                ! s% need_to_save_profiles_now = .true.
            !                ! s% need_to_update_history_now = .true.
            !                ! s% save_profiles_model_priority = 5
            !                can_save_profile_during_burst = .true.
            !                old_age = s% star_age
            !             end if
            !          end if
            !       end if
            !    end if

            ! end if

            ! if (s% x_integer_ctrl(10) .eq. 1) then
            !    if (MOD(s% model_number, 10) .eq. 0) then
            !       s% need_to_save_profiles_now = .true.
            !       s% need_to_update_history_now = .true.
            !       s% save_profiles_model_priority = 5
            !       write(*,*) '### SAVING PROFILE ###'
            !    end if
            ! end if

            ! ##### ##### ##### ##### ##### ##### ##### ##### !

         end function extras_check_model   
   
         integer function how_many_extra_history_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_columns = 3
         end function how_many_extra_history_columns
         
         
         subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            integer, intent(out) :: ierr
            real(dp) :: max_temp 
            real(dp) :: integral_N, Pi0
            real(dp), allocatable :: brunt_N(:)
            integer :: max_pos, k
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            
            if (s% x_logical_ctrl(2)) return ! initial model flag
   
            names(1) = 'delta_lgL_limit'
            vals(1) = s% delta_lgL_limit

            names(2) = 'varcontrol_target'
            vals(2) = s% varcontrol_target

            allocate(brunt_N(s% nz))

            names(3) = 'Pi0' ! Buoyancy travel time (buoyancy radius).

            brunt_N = sqrt(max(0._dp,s% brunt_N2))
            integral_N = 0.0_dp

            do k = 1, s%nz-1
               integral_N  = integral_N + brunt_N(k)*abs(s% rmid(k+1) - s% rmid(k)) / s% r(k)
            end do

            Pi0 = (2._dp*pow(pi,2))/integral_N
            vals(3) = Pi0
   
         end subroutine data_for_extra_history_columns
   
         
         integer function how_many_extra_profile_columns(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_columns = 2
         end function how_many_extra_profile_columns
         
         
         subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
            integer, intent(in) :: id, n, nz
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(nz,n)
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            integer :: k
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return

            if (s% x_logical_ctrl(2)) return ! initial model flag
            
            ! Add extra accretion heat to profile columns
            names(1) = 'Ladd'
            names(2) = "log_Ladd"

            do k = 1, s% nz
               vals(k,1) = s% extra_heat(k) %val
               vals(k,2) = LOG10( MAX(s% extra_heat(k) %val, 1d-10) )
            end do

         end subroutine data_for_extra_profile_columns
   
   
         integer function how_many_extra_history_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_history_header_items = 0
         end function how_many_extra_history_header_items


         subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_history_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return
         end subroutine data_for_extra_history_header_items


         integer function how_many_extra_profile_header_items(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            how_many_extra_profile_header_items = 1
         end function how_many_extra_profile_header_items
   

         subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
            integer, intent(in) :: id, n
            character (len=maxlen_profile_column_name) :: names(n)
            real(dp) :: vals(n)
            type(star_info), pointer :: s
            integer, intent(out) :: ierr
            ierr = 0
            call star_ptr(id,s,ierr)
            if(ierr/=0) return

            ! here is an example for adding an extra profile header item
            ! also set how_many_extra_profile_header_items
            ! names(1) = 'mixing_length_alpha'
            ! vals(1) = s% mixing_length_alpha
   
         end subroutine data_for_extra_profile_header_items
   
   
         ! returns either keep_going or terminate.
         ! note: cannot request retry; extras_check_model can do that.
         integer function extras_finish_step(id)
            integer, intent(in) :: id
            integer :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            extras_finish_step = keep_going
   
            ! to save a profile, 
               ! s% need_to_save_profiles_now = .true.
            ! to update the star log,
               ! s% need_to_update_history_now = .true.

            ! if (can_save_profile_during_burst) then
            !    write(*,*) '##### saving profile during burst #####'
            !    s% need_to_save_profiles_now = .true.
            !    s% need_to_update_history_now = .true.
            !    s% save_profiles_model_priority = 4
            ! end if

            ! ########## Leonard ##########
            ! Write logic for saving profiles
            ! #############################
   
            ! see extras_check_model for information about custom termination codes
            ! by default, indicate where (in the code) MESA terminated
            if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
         end function extras_finish_step
         
         
         subroutine extras_after_evolve(id, ierr)
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s
            ierr = 0
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
         end subroutine extras_after_evolve

         subroutine sort_dp_array(arr, n)
            implicit none
            integer, intent(in) :: n
            real(dp), intent(inout) :: arr(n)
            integer :: i, j
            real(dp) :: temp
         
            do i = 1, n - 1
               do j = i + 1, n
                  if (arr(j) < arr(i)) then
                     temp = arr(i)
                     arr(i) = arr(j)
                     arr(j) = temp
                  end if
               end do
            end do
         end subroutine

         subroutine energy_routine(id, ierr)
            use const_def
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            type (star_info), pointer :: s

            real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul, beta, Ladd
            real(dp) :: xposk, dif, xpos, thingie, thingie2, Mstar, mass_dif
            integer :: k, xposk_orig

            real(dp), parameter :: mdot_break = 3.5e-5
            real(dp), parameter :: alpha_high = 0.2
            real(dp), parameter :: alpha_low = 0.005
            real(dp), parameter :: delta_break = 2.0/3.0*1e-6
            real(dp) :: numerator, denominator

            ierr = 0
            call star_ptr(id, s, ierr)

            mass = s% mstar
            mdot = s% mstar_dot

            if (mdot == 0.) return

            radius = s% r(1)

            mcumul = 0.
            k = 0
            beta = s% x_ctrl(50)
            Ladd = beta * 5d-1 * standard_cgrav * s% star_mass * Msun * s% mstar_dot / s% r(1)           
            
            Mstar=(s% star_mass * Msun) ![g]
   
            xpos = s% x_ctrl(51)
         
            do k = 1, s% nz
               if (s% star_mdot > 0.0d0 .and. s% m(k) >= (1.0d0-xpos)*(s% star_mass)*Msun ) then
                     s% extra_heat(k) = 2.0d0 * Ladd / ( Mstar * xpos**2 ) * (s% m(k) / Mstar - (1.0d0 - xpos))
               end if
            end do
  
         end subroutine energy_routine