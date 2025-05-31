!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   PSyKAl lite code to compute the stencil adj_poly_adv_upd_lookup_kernel.
!!          PSyclone issue: #2932.
!> @details PSy-lite code required here as halo exchanges must be enforced
!!          for fields being used in calculations. The lookup table acts as
!!          a stencil operation.

MODULE invoke_adj_poly_adv_upd_lookup_mod
  USE constants_mod, ONLY: r_tran, r_def, i_def
  USE r_tran_field_mod, ONLY: r_tran_field_type, r_tran_field_proxy_type
  USE integer_field_mod, ONLY: integer_field_type, integer_field_proxy_type
  IMPLICIT NONE
  CONTAINS

  !> @brief PSy-lite code required here as halo exchanges must be enforced
  !>        for fields being used in calculations. The lookup table acts as
  !>        a stencil operation.
  SUBROUTINE invoke_adj_poly_adv_upd_lookup(advective, reconstruction, &
                                            lookup_poly_adv_upd_field, num_sets_poly_adv_upd_field, &
                                            wind, wind_dir, &
                                            nsets, nindices, stencil_extent)
    USE adj_poly_adv_upd_lookup_kernel_mod, ONLY: adj_poly_adv_upd_lookup_code
    USE mesh_mod, ONLY: mesh_type
    USE stencil_2D_dofmap_mod, ONLY: stencil_2D_dofmap_type, STENCIL_2D_CROSS
    INTEGER(KIND=i_def), intent(in) :: nsets, nindices
    TYPE(r_tran_field_type), intent(in) :: advective, reconstruction, wind, wind_dir
    TYPE(integer_field_type), intent(in) :: lookup_poly_adv_upd_field, num_sets_poly_adv_upd_field
    INTEGER(KIND=i_def), intent(in) :: stencil_extent
    INTEGER(KIND=i_def) df
    INTEGER(KIND=i_def) cell
    INTEGER(KIND=i_def) loop1_start, loop1_stop
    INTEGER(KIND=i_def) loop0_start, loop0_stop
    INTEGER(KIND=i_def) nlayers_advective
    INTEGER(KIND=i_def), pointer, dimension(:) :: num_sets_poly_adv_upd_field_data
    INTEGER(KIND=i_def), pointer, dimension(:) :: lookup_poly_adv_upd_field_data
    TYPE(integer_field_proxy_type) lookup_poly_adv_upd_field_proxy, num_sets_poly_adv_upd_field_proxy
    REAL(KIND=r_tran), pointer, dimension(:) :: wind_dir_data
    REAL(KIND=r_tran), pointer, dimension(:) :: wind_data
    REAL(KIND=r_tran), pointer, dimension(:) :: reconstruction_data
    REAL(KIND=r_tran), pointer, dimension(:) :: advective_data
    TYPE(r_tran_field_proxy_type) advective_proxy, reconstruction_proxy, wind_proxy, wind_dir_proxy
    INTEGER(KIND=i_def), pointer :: map_adspc1_reconstruction(:,:), &
&map_adspc2_lookup_poly_adv_upd_field(:,:), map_adspc3_num_sets_poly_adv_upd_field(:,:), &
&map_w2(:,:), map_wtheta(:,:)
    INTEGER(KIND=i_def) ndf_wtheta, undf_wtheta, ndf_adspc1_reconstruction, undf_adspc1_reconstruction, &
&ndf_adspc2_lookup_poly_adv_upd_field, undf_adspc2_lookup_poly_adv_upd_field, ndf_adspc3_num_sets_poly_adv_upd_field, &
&undf_adspc3_num_sets_poly_adv_upd_field, ndf_w2, undf_w2, ndf_aspc1_advective, undf_aspc1_advective
    INTEGER(KIND=i_def) max_halo_depth_mesh
    TYPE(mesh_type), pointer :: mesh
    INTEGER(KIND=i_def) wind_max_branch_length
    INTEGER(KIND=i_def), pointer :: wind_stencil_size(:,:)
    INTEGER(KIND=i_def), pointer :: wind_stencil_dofmap(:,:,:,:)
    TYPE(stencil_2D_dofmap_type), pointer :: wind_stencil_map

    nullify( num_sets_poly_adv_upd_field_data, lookup_poly_adv_upd_field_data, &
             wind_dir_data, wind_data, reconstruction_data, advective_data, &
             map_adspc1_reconstruction, map_adspc2_lookup_poly_adv_upd_field, &
             map_adspc3_num_sets_poly_adv_upd_field, map_w2, map_wtheta, &
             mesh, &
             wind_stencil_size, wind_stencil_dofmap, wind_stencil_map )

    !
    ! Initialise field and/or operator proxies
    !
    advective_proxy = advective%get_proxy()
    advective_data => advective_proxy%data
    reconstruction_proxy = reconstruction%get_proxy()
    reconstruction_data => reconstruction_proxy%data
    lookup_poly_adv_upd_field_proxy = lookup_poly_adv_upd_field%get_proxy()
    lookup_poly_adv_upd_field_data => lookup_poly_adv_upd_field_proxy%data
    num_sets_poly_adv_upd_field_proxy = num_sets_poly_adv_upd_field%get_proxy()
    num_sets_poly_adv_upd_field_data => num_sets_poly_adv_upd_field_proxy%data
    wind_proxy = wind%get_proxy()
    wind_data => wind_proxy%data
    wind_dir_proxy = wind_dir%get_proxy()
    wind_dir_data => wind_dir_proxy%data
    !
    ! Initialise number of layers
    !
    nlayers_advective = advective_proxy%vspace%get_nlayers()
    !
    ! Create a mesh object
    !
    mesh => advective_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Initialise stencil dofmaps
    !
    wind_stencil_map => wind_proxy%vspace%get_stencil_2D_dofmap(STENCIL_2D_CROSS,stencil_extent)
    wind_max_branch_length = stencil_extent + 1_i_def
    wind_stencil_dofmap => wind_stencil_map%get_whole_dofmap()
    wind_stencil_size => wind_stencil_map%get_stencil_sizes()
    !
    ! Look-up dofmaps for each function space
    !
    map_wtheta => advective_proxy%vspace%get_whole_dofmap()
    map_adspc1_reconstruction => reconstruction_proxy%vspace%get_whole_dofmap()
    map_adspc2_lookup_poly_adv_upd_field => lookup_poly_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_adspc3_num_sets_poly_adv_upd_field => num_sets_poly_adv_upd_field_proxy%vspace%get_whole_dofmap()
    map_w2 => wind_proxy%vspace%get_whole_dofmap()
    !
    ! Initialise number of DoFs for wtheta
    !
    ndf_wtheta = advective_proxy%vspace%get_ndf()
    undf_wtheta = advective_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc1_reconstruction
    !
    ndf_adspc1_reconstruction = reconstruction_proxy%vspace%get_ndf()
    undf_adspc1_reconstruction = reconstruction_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc2_lookup_poly_adv_upd_field
    !
    ndf_adspc2_lookup_poly_adv_upd_field = lookup_poly_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc2_lookup_poly_adv_upd_field = lookup_poly_adv_upd_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for adspc3_num_sets_poly_adv_upd_field
    !
    ndf_adspc3_num_sets_poly_adv_upd_field = num_sets_poly_adv_upd_field_proxy%vspace%get_ndf()
    undf_adspc3_num_sets_poly_adv_upd_field = num_sets_poly_adv_upd_field_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for w2
    !
    ndf_w2 = wind_proxy%vspace%get_ndf()
    undf_w2 = wind_proxy%vspace%get_undf()
    !
    ! Initialise number of DoFs for aspc1_advective
    !
    ndf_aspc1_advective = advective_proxy%vspace%get_ndf()
    undf_aspc1_advective = advective_proxy%vspace%get_undf()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = mesh%get_last_edge_cell()
    loop1_start = 1
    loop1_stop = advective_proxy%vspace%get_last_dof_halo(1)
    !
    ! Call kernels and communication routines
    !
    IF (wind_proxy%is_dirty(depth=stencil_extent)) THEN
      CALL wind_proxy%halo_exchange(depth=stencil_extent)
    END IF

    ! EDIT(JC): Call halo exchange for stencil-like accesses.
    if (advective_proxy%is_dirty(depth=stencil_extent)) then
      call advective_proxy%halo_exchange(depth=stencil_extent)
    end if
    if (wind_dir_proxy%is_dirty(depth=stencil_extent)) then
      call wind_dir_proxy%halo_exchange(depth=stencil_extent)
    end if

    !$omp parallel default(shared), private(cell)
    !$omp do schedule(static)
    DO cell = loop0_start, loop0_stop, 1
      CALL adj_poly_adv_upd_lookup_code(nlayers_advective, advective_data, reconstruction_data, lookup_poly_adv_upd_field_data, &
&num_sets_poly_adv_upd_field_data, wind_data, wind_stencil_size(:,cell), wind_max_branch_length, wind_stencil_dofmap(:,:,:,cell), &
&wind_dir_data, nsets, nindices, ndf_wtheta, undf_wtheta, map_wtheta(:,cell), ndf_adspc1_reconstruction, &
&undf_adspc1_reconstruction, map_adspc1_reconstruction(:,cell), ndf_adspc2_lookup_poly_adv_upd_field, &
&undf_adspc2_lookup_poly_adv_upd_field, map_adspc2_lookup_poly_adv_upd_field(:,cell), ndf_adspc3_num_sets_poly_adv_upd_field, &
&undf_adspc3_num_sets_poly_adv_upd_field, map_adspc3_num_sets_poly_adv_upd_field(:,cell), ndf_w2, undf_w2, map_w2(:,cell))
    END DO
    !$omp end do
    !$omp end parallel

    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    CALL advective_proxy%set_dirty()
    CALL reconstruction_proxy%set_dirty()
    !
    DO df = loop1_start, loop1_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      advective_data(df) = 0.0_r_def
    END DO
    !
    ! Set halos dirty/clean for fields modified in the above loop
    !
    CALL advective_proxy%set_dirty()
    CALL advective_proxy%set_clean(1)
    !
    !
  END SUBROUTINE invoke_adj_poly_adv_upd_lookup

END MODULE invoke_adj_poly_adv_upd_lookup_mod
