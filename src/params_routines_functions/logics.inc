logical :: initial_increase = .true.
logical :: final_mass_reached = .false.
logical :: seed_phase_over = .false.
logical :: sigmoid_mass_integrated = .false.
logical :: set_initial_composition = .false.
logical :: reached_target_age = .false.
logical :: decay_at_end_of_accretion_can_begin = .false.
logical :: can_save_profile_during_burst = .false.
logical :: final_mass_reached_within_tolerance = .false.
logical, save :: above_1_5
! logical :: retry_flag = .false.