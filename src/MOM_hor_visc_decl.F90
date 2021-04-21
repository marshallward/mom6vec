subroutine horizontal_viscosity(u, v, h, diffu, diffv, MEKE, VarMix, G, GV, US, &
                                CS, OBC, BT, TD, ADp)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                       !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                       !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  type(MEKE_type),               pointer     :: MEKE   !< Pointer to a structure containing fields
                                                       !! related to Mesoscale Eddy Kinetic Energy.
  type(VarMix_CS),               pointer     :: VarMix !< Pointer to a structure with fields that
                                                       !! specify the spatially variable viscosities
  type(unit_scale_type),         intent(in)  :: US     !< A dimensional unit scaling type
  type(hor_visc_CS),             pointer     :: CS     !< Control structure returned by a previous
                                                       !! call to hor_visc_init.
  type(ocean_OBC_type), optional, pointer    :: OBC    !< Pointer to an open boundary condition type
  type(barotropic_CS),  optional, pointer    :: BT     !< Pointer to a structure containing
                                                       !! barotropic velocities.
  type(thickness_diffuse_CS), optional, pointer :: TD  !< Pointer to a structure containing
                                                       !! thickness diffusivities.
  type(accel_diag_ptrs), optional, pointer :: ADp      !< Acceleration diagnostic pointers

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    Del2u, &      ! The u-compontent of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_u, &        ! Thickness interpolated to u points [H ~> m or kg m-2].
    vort_xy_dy, & ! y-derivative of vertical vorticity (d/dy(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    div_xx_dx, &  ! x-derivative of horizontal divergence (d/dx(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    ubtav         ! zonal barotropic vel. ave. over baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G)) :: &
    Del2v, &      ! The v-compontent of the Laplacian of velocity [L-1 T-1 ~> m-1 s-1]
    h_v, &        ! Thickness interpolated to v points [H ~> m or kg m-2].
    vort_xy_dx, & ! x-derivative of vertical vorticity (d/dx(dv/dx - du/dy)) [L-1 T-1 ~> m-1 s-1]
    div_xx_dy, &  ! y-derivative of horizontal divergence (d/dy(du/dx + dv/dy)) [L-1 T-1 ~> m-1 s-1]
    vbtav         ! meridional barotropic vel. ave. over baroclinic time-step [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: &
    dudx_bt, dvdy_bt, & ! components in the barotropic horizontal tension [T-1 ~> s-1]
    div_xx, &     ! Estimate of horizontal divergence at h-points [T-1 ~> s-1]
    sh_xx, &      ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    sh_xx_bt, &   ! barotropic horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    str_xx,&      ! str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
    str_xx_GME,&  ! smoothed diagonal term in the stress tensor from GME [H L2 T-2 ~> m3 s-2 or kg s-2]
    bhstr_xx, &   ! A copy of str_xx that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    FrictWorkIntz, & ! depth integrated energy dissipated by lateral friction [R L2 T-3 ~> W m-2]
    grad_vort_mag_h, & ! Magnitude of vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_h_2d, & ! Magnitude of 2d vorticity gradient at h-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_h, & ! Laplacian of vorticity at h-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_h, &     ! Magnitude of divergence gradient at h-points [L-1 T-1 ~> m-1 s-1]
    dudx, dvdy, &    ! components in the horizontal tension [T-1 ~> s-1]
    grad_vel_mag_h, & ! Magnitude of the velocity gradient tensor squared at h-points [T-2 ~> s-2]
    grad_vel_mag_bt_h, & ! Magnitude of the barotropic velocity gradient tensor squared at h-points [T-2 ~> s-2]
    grad_d2vel_mag_h, & ! Magnitude of the Laplacian of the velocity vector, squared [L-2 T-2 ~> m-2 s-2]
    boundary_mask_h ! A mask that zeroes out cells with at least one land edge [nondim]

  real, dimension(SZIB_(G),SZJB_(G)) :: &
    dvdx, dudy, & ! components in the shearing strain [T-1 ~> s-1]
    dDel2vdx, dDel2udy, & ! Components in the biharmonic equivalent of the shearing strain [L-2 T-1 ~> m-2 s-1]
    dvdx_bt, dudy_bt,   & ! components in the barotropic shearing strain [T-1 ~> s-1]
    sh_xy,  &     ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
    sh_xy_bt, &   ! barotropic horizontal shearing strain (du/dy + dv/dx) inc. metric terms [T-1 ~> s-1]
    str_xy, &     ! str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2]
    str_xy_GME, & ! smoothed cross term in the stress tensor from GME [H L2 T-2 ~> m3 s-2 or kg s-2]
    bhstr_xy, &   ! A copy of str_xy that only contains the biharmonic contribution [H L2 T-2 ~> m3 s-2 or kg s-2]
    vort_xy, &    ! Vertical vorticity (dv/dx - du/dy) including metric terms [T-1 ~> s-1]
    Leith_Kh_q, & ! Leith Laplacian viscosity at q-points [L2 T-1 ~> m2 s-1]
    Leith_Ah_q, & ! Leith bi-harmonic viscosity at q-points [L4 T-1 ~> m4 s-1]
    grad_vort_mag_q, & ! Magnitude of vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vort_mag_q_2d, & ! Magnitude of 2d vorticity gradient at q-points [L-1 T-1 ~> m-1 s-1]
    Del2vort_q, & ! Laplacian of vorticity at q-points [L-2 T-1 ~> m-2 s-1]
    grad_div_mag_q, &  ! Magnitude of divergence gradient at q-points [L-1 T-1 ~> m-1 s-1]
    grad_vel_mag_q, &  ! Magnitude of the velocity gradient tensor squared at q-points [T-2 ~> s-2]
    hq, &         ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
                  ! This form guarantees that hq/hu < 4.
    grad_vel_mag_bt_q, &  ! Magnitude of the barotropic velocity gradient tensor squared at q-points [T-2 ~> s-2]
    boundary_mask_q ! A mask that zeroes out cells with at least one land edge [nondim]

  real, dimension(SZIB_(G),SZJB_(G),SZK_(GV)) :: &
    Ah_q, &      ! biharmonic viscosity at corner points [L4 T-1 ~> m4 s-1]
    Kh_q, &      ! Laplacian viscosity at corner points [L2 T-1 ~> m2 s-1]
    vort_xy_q, & ! vertical vorticity at corner points [T-1 ~> s-1]
    sh_xy_q,   & ! horizontal shearing strain at corner points [T-1 ~> s-1]
    GME_coeff_q, &  !< GME coeff. at q-points [L2 T-1 ~> m2 s-1]
    ShSt         ! A diagnostic array of shear stress [T-1 ~> s-1].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)+1) :: &
    KH_u_GME  !< interface height diffusivities in u-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)+1) :: &
    KH_v_GME  !< interface height diffusivities in v-columns [L2 T-1 ~> m2 s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    Ah_h, &          ! biharmonic viscosity at thickness points [L4 T-1 ~> m4 s-1]
    Kh_h, &          ! Laplacian viscosity at thickness points [L2 T-1 ~> m2 s-1]
    FrictWork, &     ! work done by MKE dissipation mechanisms [R L2 T-3 ~> W m-2]
    FrictWork_GME, & ! work done by GME [R L2 T-3 ~> W m-2]
    div_xx_h,      & ! horizontal divergence [T-1 ~> s-1]
    sh_xx_h,       & ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
    NoSt             ! A diagnostic array of normal stress [T-1 ~> s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: &
    grid_Re_Kh, &    ! Grid Reynolds number for Laplacian horizontal viscosity at h points [nondim]
    grid_Re_Ah, &    ! Grid Reynolds number for Biharmonic horizontal viscosity at h points [nondim]
    GME_coeff_h      ! GME coeff. at h-points [L2 T-1 ~> m2 s-1]
  real :: AhSm       ! Smagorinsky biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: AhLth      ! 2D Leith biharmonic viscosity [L4 T-1 ~> m4 s-1]
  real :: mod_Leith  ! nondimensional coefficient for divergence part of modified Leith
                     ! viscosity. Here set equal to nondimensional Laplacian Leith constant.
                     ! This is set equal to zero if modified Leith is not used.
  real :: Shear_mag_bc  ! Shear_mag value in backscatter [T-1 ~> s-1]
  real :: sh_xx_sq   ! Square of tension (sh_xx) [T-2 ~> s-2]
  real :: sh_xy_sq   ! Square of shearing strain (sh_xy) [T-2 ~> s-2]
  real :: h2uq, h2vq ! temporary variables [H2 ~> m2 or kg2 m-4].
  real :: hu, hv     ! Thicknesses interpolated by arithmetic means to corner
                     ! points; these are first interpolated to u or v velocity
                     ! points where masks are applied [H ~> m or kg m-2].
  real :: h_neglect  ! thickness so small it can be lost in roundoff and so neglected [H ~> m or kg m-2]
  real :: h_neglect3 ! h_neglect^3 [H3 ~> m3 or kg3 m-6]
  real :: h_min      ! Minimum h at the 4 neighboring velocity points [H ~> m]
  real :: Kh_scale  ! A factor between 0 and 1 by which the horizontal
                    ! Laplacian viscosity is rescaled [nondim]
  real :: RoScl     ! The scaling function for MEKE source term [nondim]
  real :: FatH      ! abs(f) at h-point for MEKE source term [T-1 ~> s-1]
  real :: local_strain ! Local variable for interpolating computed strain rates [T-1 ~> s-1].
  real :: meke_res_fn ! A copy of the resolution scaling factor if being applied to MEKE. Otherwise =1.
  real :: GME_coeff ! The GME (negative) viscosity coefficient [L2 T-1 ~> m2 s-1]
  real :: GME_coeff_limiter ! Maximum permitted value of the GME coefficient [L2 T-1 ~> m2 s-1]
  real :: FWfrac    ! Fraction of maximum theoretical energy transfer to use when scaling GME coefficient [nondim]
  real :: DY_dxBu   ! Ratio of meridional over zonal grid spacing at vertices [nondim]
  real :: DX_dyBu   ! Ratio of zonal over meridiononal grid spacing at vertices [nondim]
  real :: DY_dxCv   ! Ratio of meridional over zonal grid spacing at faces [nondim]
  real :: DX_dyCu   ! Ratio of zonal over meridional grid spacing at faces [nondim]
  real :: Sh_F_pow  ! The ratio of shear over the absolute value of f raised to some power and rescaled [nondim]
  real :: backscat_subround ! The ratio of f over Shear_mag that is so small that the backscatter
                    ! calculation gives the same value as if f were 0 [nondim].
  real :: H0_GME    ! Depth used to scale down GME coefficient in shallow areas [Z ~> m]
  real :: KE        ! Local kinetic energy [L2 T-2 ~> m2 s-2]
  real :: d_del2u   ! dy-weighted Laplacian(u) diff in x [L-2 T-1 ~> m-2 s-1]
  real :: d_del2v   ! dx-weighted Laplacian(v) diff in y [L-2 T-1 ~> m-2 s-1]
  real :: d_str     ! Stress tensor update [H L2 T-2 ~> m3 s-2 or kg s-2]
  real :: grad_vort ! Vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grad_vort_qg ! QG-based vorticity gradient magnitude [L-1 T-1 ~> m-1 s-1]
  real :: grid_Kh   ! Laplacian viscosity bound by grid [L2 T-1 ~> m2 s-1]
  real :: grid_Ah   ! Biharmonic viscosity bound by grid [L4 T-1 ~> m4 s-1]

  logical :: rescale_Kh, legacy_bound
  logical :: find_FrictWork
  logical :: apply_OBC = .false.
  logical :: use_MEKE_Ku
  logical :: use_MEKE_Au
  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k, n
  real :: inv_PI3, inv_PI2, inv_PI6

  ! Fields evaluated on active layers, used for constructing 3D stress fields
  ! NOTE: The position of these declarations can impact performance, due to the
  !   very large number of stack arrays in this function.  Move with caution!
  real, dimension(SZIB_(G),SZJB_(G)) :: &
    Ah, &           ! biharmonic viscosity (h or q) [L4 T-1 ~> m4 s-1]
    Kh, &           ! Laplacian  viscosity [L2 T-1 ~> m2 s-1]
    Shear_mag, &    ! magnitude of the shear [T-1 ~> s-1]
    vert_vort_mag, &  ! magnitude of the vertical vorticity gradient [L-1 T-1 ~> m-1 s-1]
    hrat_min, &     ! h_min divided by the thickness at the stress point (h or q) [nondim]
    visc_bound_rem  ! fraction of overall viscous bounds that remain to be applied [nondim]

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    hf_diffu_2d, &    ! Depth sum of hf_diffu, hf_diffv [L T-2 ~> m s-2]
    intz_diffu_2d     ! Depth-integral of diffu [L2 T-2 ~> m2 s-2]

  real, dimension(SZI_(G),SZJB_(G)) :: &
    hf_diffv_2d, &    ! Depth sum of hf_diffu, hf_diffv [L T-2 ~> m s-2]
    intz_diffv_2d     ! Depth-integral of diffv [L2 T-2 ~> m2 s-2]

