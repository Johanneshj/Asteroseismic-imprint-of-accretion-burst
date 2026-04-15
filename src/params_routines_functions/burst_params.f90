! ### Module containing burst parameters ###
module burst_parameters_single
    use const_def
    implicit none

    type :: burst_parameters_for_single_burst
        logical  :: burst_over, no_more_bursts
        integer  :: T, current_burst, number_of_bursts
        real(dp) :: mdot_ini, final_mass, sigmoid_mass, 
        real(dp) :: target_age, t_end_accretion
        real(dp) :: burst_rise_time
        real(dp) :: burst_plateau_time
        real(dp) :: burst_decay_time
        real(dp) :: burst_mdot
        real(dp) :: burst_mass_start
        real(dp) :: t_end_burst
        real(dp) :: mass_end_burst
    end type burst_parameters_for_single_burst

    type(burst_parameters_for_single_burst) :: bp

end module burst_parameters_single
