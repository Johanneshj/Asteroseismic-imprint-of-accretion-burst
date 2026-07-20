! Saving profiles
real(dp), save :: old_age = -1d99
real(dp), save :: age_post_burst = -1d99
real(dp), save :: age_pre_burst = -1d99
real(dp), save :: age_approach_zams = -1d99
real(dp), save :: termination_age = -1d99
real(dp), save :: old_mass = -1d99

! Parameters for initial sigmoid increase
! To be defined in params_routines_functions/single_burst.inc
! --> function integrate_sigmoid_mass(sigmoid_mass, mdot_ini) result(sigmoid_m)
real(dp) :: sigmoid_init
real(dp) :: years_for_initial_sigmoid_increase 

real(dp) :: age_init = -1d99

! Burst parameters
type :: burst_parameters_single
    logical  :: burst_over, no_more_bursts
    integer  :: T, current_burst, number_of_bursts
    real(dp) :: mdot_ini, final_mass, sigmoid_mass 
    real(dp) :: target_age, t_end_accretion
    real(dp) :: burst_rise_time
    real(dp) :: burst_plateau_time
    real(dp) :: burst_decay_time
    real(dp) :: burst_mdot
    real(dp) :: burst_mass_start
    real(dp) :: t_end_burst
    real(dp) :: mass_end_burst
end type burst_parameters_single
type(burst_parameters_single), save :: bp

! Parameters for timestep controls
real(dp) :: conv_core_top = -1
real(dp) :: amount_of_mass_to_begin_decreasing_timestep_before_burst = 1d-3 ! Msun
! real(dp), dimension(:), allocatable, save :: gradr_old(:), mass_old(:), gradr_interpolated(:)

real(dp) :: mass_end_of_seed_phase = 4d-1
real(dp) :: central_c12_decrease_to_mark_end_of_seed_phase = 0.99
real(dp) :: center_c12_ini = -1

real(dp) :: target_age_for_decay = -1d99