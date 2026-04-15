function set_mesh_factor_towards_burst(id, ierr, min_factor, max_factor) result(mesh_factor)

    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    real(dp) :: min_factor, max_factor, mesh_factor
    real(dp) :: age_in_secs, elapsed_time
    real(dp) :: midpoint, steepness, term_1, term_2, front_factor, duration, when_to_increase
    real(dp) :: varcontrol_profile, minimum_varcontrol_profile

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! if (.not. seed_phase_over) return

    age_in_secs = s% star_age * secyer
    elapsed_time = age_in_secs - bp% target_age

    ! increase mesh along with decreasing varcontrol_target (see constrain_model_with_physics.inc)
    when_to_increase = (5d0 * amount_of_mass_to_begin_decreasing_timestep_before_burst) / &
                            (bp% mdot_ini * secyer / Msun) * secyer
    midpoint = when_to_increase / 2d0
    steepness = 0.07d0 * midpoint

    duration = bp% burst_rise_time + &
                bp% burst_plateau_time + &
                bp% burst_decay_time

    front_factor = min_factor-max_factor
    term_1 = 1d0 / (1d0 + EXP((midpoint - elapsed_time) / steepness))
    term_2 = 1d0 / (1d0 + EXP(((midpoint + duration) + elapsed_time) / steepness))

    mesh_factor = max_factor + front_factor * (term_1 + term_2) ! Increase and decrease mesh factor towards and after burst
    ! write(*,*) 'mesh factor:', mesh_factor

end function set_mesh_factor_towards_burst


subroutine how_many_mesh_fcns(id, n)
    integer, intent(in) :: id
    integer, intent(out) :: n
    n = 1
end subroutine how_many_mesh_fcns
      
      
subroutine gradr_grada_mesh_fcn_data( &
        id, nfcns, names, gval_is_xa_function, vals1, ierr)
    integer, intent(in) :: id
    integer, intent(in) :: nfcns
    character (len=*) :: names(:)
    logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
    real(dp), pointer :: vals1(:) ! =(nz, nfcns)
    integer, intent(out) :: ierr

    real(dp), pointer :: vals(:,:)
    type (star_info), pointer :: s
    integer :: nz, k, idx
    real(dp) :: weight, width, center
    real(dp), parameter :: maxval = 700d0 ! max value for tanh from crlibm

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    ! if (.not. seed_phase_over) return ! ignore seed phase for stability
    ! write(*,*) 'cc mass:', s% mass_conv_core
    weight = s% x_ctrl(7)
    width = s% x_ctrl(8)
    center = s% x_ctrl(9)

    names(1) = 'gradr_grada_function'
    gval_is_xa_function(1) = .false.

    ! idx = s% mixing_region_bottom(2) ! bottom of mixing region 2 should correspond to envelope
    ! write(*,*) 'idx:', idx
    nz = s% nz
    vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)

    do k=1,nz
        vals(k,1) = weight*tanh(min(maxval,(s% gradr(k)-s% grada(k)-center)/width))*width
    end do

end subroutine gradr_grada_mesh_fcn_data


subroutine control_mesh(id, ierr)

    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    if (s% x_integer_ctrl(10) .ne. 1) return

    ! Gradually increase mesh factors 
    s% xa_function_weight(1:2)     = set_mesh_factor(id, ierr, 60d0, 80d0) ! he3
    s% mesh_dlogX_dlogP_extra(1:2) = set_mesh_factor(id, ierr, 0.2d0, 0.3d0) ! c12
    s% T_function1_weight          = set_mesh_factor(id, ierr, 50d0, 75d0)
    s% P_function_weight           = set_mesh_factor(id, ierr, 20d0, 50d0)
    s% R_function_weight           = set_mesh_factor(id, ierr, 10d0, 20d0)

    ! This does not work right now!

end subroutine control_mesh

function set_mesh_factor(id, ierr, factor_ini, factor_end) result(sigmoid)
    use star_def
    use const_def

    implicit none
    integer, intent(in) :: id
    integer :: ierr
    integer :: step, i
    type(star_info), pointer :: s

    real(dp) :: years_for_increase, secs_for_increase, age_in_secs, target_age
    real(dp) :: factor, factor_end, factor_ini, sigmoid, exponential, steepness, midpoint, elapsed_time

    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    if (s% x_integer_ctrl(10) .ne. 1) return

    age_in_secs = s% star_age * secyer
    if (age_init .lt. 0) age_init = s% star_age * secyer

    years_for_increase = 10d0
    secs_for_increase = years_for_increase * secyer

    target_age = age_init + (secs_for_increase * secyer)
    elapsed_time = target_age - age_in_secs

    midpoint = secs_for_increase / 2d0
    steepness = midpoint * 0.1d0

    factor = (factor_end - factor_ini)
    exponential = EXP(-(midpoint - elapsed_time)/steepness)

    sigmoid = factor_ini + (factor_end - factor_ini) * 1d0 / (1d0 + exponential)
    sigmoid = sigmoid
    ! write(*,*) 'sigmoid:', sigmoid
    ! write(*,*) 'use other?', s% use_other_mesh_functions

end function set_mesh_factor