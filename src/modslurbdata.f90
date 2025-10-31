module modslurbdata
    use modprecision, only: field_r
    implicit none
    save
    public

    logical :: lslurb            ! On/off switch LSM


    !-- Derived type for the SLUrb model.
    type surf_slurb

        real(field_r), allocatable ::  dt_max(:,:)  !< time step limit for model physical processes (s)

        real(field_r), allocatable ::  dz_road(:,:,:)  !< road layer thickness (m)
        real(field_r), allocatable ::  dz_roof(:,:,:)  !< roof layer thickness (m)
        real(field_r), allocatable ::  dz_wall(:,:,:)  !< wall layer thickness (m)
        real(field_r), allocatable ::  dz_win(:,:,:)   !< window layer thickness (m)
        real(field_r), allocatable ::  zw_win(:,:,:)   !< cumulative window thickness (m)
        !
        !--    Tile-aggregated quantities.
        real(field_r), allocatable ::  albedo_urb(:,:)       !< effective urban albedo (.)
        real(field_r), allocatable ::  emiss_urb(:,:)        !< effective urban emissivity (.)
        real(field_r), allocatable ::  ol_urb(:,:)           !< urban Obukhov length (m)
        real(field_r), allocatable ::  qsws_urb(:,:)         !< total urban latent heat flux (W m^-2)
        real(field_r), allocatable ::  rad_lw_in_urb(:,:)   !< incoming longwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_lw_out_urb(:,:)  !< outgoing longwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_urb(:,:)   !< incoming shortwave radiation (W m^-2)
        real(field_r), allocatable ::  rad_sw_out_urb(:,:)  !< outgoing shortwave radiation (W m^-2)
        real(field_r), allocatable ::  ram_urb(:,:)         !< urban aerodynamic resistance for momentum (s m^-1)
        real(field_r), allocatable ::  rib_urb(:,:)         !< urban bulk-Richardson number (.)
        real(field_r), allocatable ::  shf_urb(:,:)         !< total urban sensible heat flux (W m^-2)
        real(field_r), allocatable ::  t_2m_urb(:,:)        !< urban 2-metre temperature (extrapolated) (K)
        real(field_r), allocatable ::  t_c_urb(:,:)         !< complete (area-weighted) urban surface temperature (K)
        real(field_r), allocatable ::  t_h_urb(:,:)         !< effective urban surface temperature (K)
        real(field_r), allocatable ::  t_rad_urb(:,:)       !< urban radiative surface temperature (K)
        real(field_r), allocatable ::  usws_urb(:,:)        !< urban momentum flux (u-component) (kg m^-1 s^-2)
        real(field_r), allocatable ::  vsws_urb(:,:)        !< urban momentum flux (v-component) (kg m^-1 s^-2)
        real(field_r), allocatable ::  thlskin(:,:)         !< urban skin liquid water potential temperature (K)
        real(field_r), allocatable ::  qtskin(:,:)          !< urban skin specific humidity TODOSELF (kg kg^-1)
        !
        !--    Model prognostic variables.
        real(field_r), pointer, contiguous ::  m_liq_road(:,:)    !< liquid water reservoir on roads (m^3 m^-2)
        real(field_r), pointer, contiguous ::  m_liq_road_p(:,:)  !< prog. liquid water reservoir on roads (m^3 m^-2)
        real(field_r), pointer, contiguous ::  m_liq_roof(:,:)    !< liquid water reservoir on roofs (m^3 m^-2)
        real(field_r), pointer, contiguous ::  m_liq_roof_p(:,:)  !< prog. liquid water reservoir on roofs (m^3 m^-2)
        real(field_r), pointer, contiguous ::  q_can(:,:)         !< canyon mixing ratio (kg kg^-1)
        real(field_r), pointer, contiguous ::  q_can_p(:,:)       !< prognostic canyon mixing ratio (kg kg^-1)
        real(field_r), pointer, contiguous ::  t_can(:,:)         !< canyon air temperature (K)
        real(field_r), pointer, contiguous ::  t_can_p(:,:)       !< prog. canyon temperature (K)

        real(field_r), pointer, contiguous ::  t_road(:,:,:)      !< road temperature (K)
        real(field_r), pointer, contiguous ::  t_road_p(:,:,:)    !< prog. road temperature (K)
        real(field_r), pointer, contiguous ::  t_roof(:,:,:)      !< roof temperature (K)
        real(field_r), pointer, contiguous ::  t_roof_p(:,:,:)    !< prog. roof temperature (K)
        real(field_r), pointer, contiguous ::  t_wall_a(:,:,:)    !< wall A temperature (K)
        real(field_r), pointer, contiguous ::  t_wall_a_p(:,:,:)  !< prog. wall A temperature (K)
        real(field_r), pointer, contiguous ::  t_wall_b(:,:,:)    !< wall B temperature (K)
        real(field_r), pointer, contiguous ::  t_wall_b_p(:,:,:)  !< prog. wall B temperature (K)
        real(field_r), pointer, contiguous ::  t_win_a(:,:,:)     !< window A temperature (K)
        real(field_r), pointer, contiguous ::  t_win_a_p(:,:,:)   !< prog. window A temperature (K)
        real(field_r), pointer, contiguous ::  t_win_b(:,:,:)     !< window B temperature (K)
        real(field_r), pointer, contiguous ::  t_win_b_p(:,:,:)   !< prog. window B temperature (K)
        !
        !--    Tendencies of the prognostic variables.
        real(field_r), allocatable ::  tm_liq_road(:,:)  !< road liquid water reservoir tendency (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tm_liq_roof(:,:)  !< roof liquid water reservoir tendency (m^3 m^-2 s^-1)
        real(field_r), allocatable ::  tq_can(:,:)       !< canyon mixing ratio tendency (kg kg^-1 s^-1)
        real(field_r), allocatable ::  tt_can(:,:)       !< canyon temperature tendency (K s^-1)

        real(field_r), allocatable ::  tt_road(:,:,:)    !< road temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_roof(:,:,:)    !< roof temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_wall_a(:,:,:)  !< wall A temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_wall_b(:,:,:)  !< wall B temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_win_a(:,:,:)   !< window A temperature tendency (K s^-1)
        real(field_r), allocatable ::  tt_win_b(:,:,:)   !< window B temperature tendency (K s^-1)
        !
        !--    Diagnostic surface thermodynamic variables.
        real(field_r), allocatable ::  pt_road(:,:)    !< road surface potential temperature (K)
        real(field_r), allocatable ::  pt_roof(:,:)    !< roof surface potential temperature (K)
        real(field_r), allocatable ::  pt_wall_a(:,:)  !< wall A surface potential temperature (K)
        real(field_r), allocatable ::  pt_wall_b(:,:)  !< wall B surface potential temperature (K)
        real(field_r), allocatable ::  pt_win_a(:,:)   !< window A surface potential temperature (K)
        real(field_r), allocatable ::  pt_win_b(:,:)   !< window A surface potential temperature (K)
        real(field_r), allocatable ::  q_road(:,:)     !< road surface mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  q_roof(:,:)     !< roof surface mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  qs_road(:,:)    !< road surface saturation mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  qs_roof(:,:)    !< roof surface saturation mixing ratio (kg kg^-1)
        real(field_r), allocatable ::  vpt_road(:,:)   !< road surface virtual potential temperature (K)
        real(field_r), allocatable ::  vpt_roof(:,:)   !< roof surface virtual potential temperature (K)
        !
        !--    Diagnostic internal sensible heat fluxes.
        real(field_r), allocatable ::  shf_can(:,:)       !< sensible heat flux between the street canyon and the atmosphere (W m^-2)
        real(field_r), allocatable ::  shf_external(:,:)  !< sensible heat flux external to the model (e.g. industry) (W m^-2)
        real(field_r), allocatable ::  shf_road(:,:)      !< road surface sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_roof(:,:)      !< roof surface sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_traffic(:,:)   !< traffic sensible heat flux (input-only) (W m^-2)
        real(field_r), allocatable ::  shf_wall_a(:,:)    !< wall A sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_wall_b(:,:)    !< wall B sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_win_a(:,:)     !< window A sensible heat flux (W m^-2)
        real(field_r), allocatable ::  shf_win_b(:,:)     !< window B sensible heat flux (W m^-2)
        !
        !--    Diagnostic internal latent heat fluxes.
        real(field_r), allocatable ::  qsws_can(:,:)       !< latent heat flux between the street canyon and the atmosphere (W m^-2)
        real(field_r), allocatable ::  qsws_external(:,:)  !< latent heat flux external to the model (e.g. industry) (W m^-2)
        real(field_r), allocatable ::  qsws_liq_road(:,:)  !< roof latent heat flux (liquid incl. precipitation) (W m^-2)
        real(field_r), allocatable ::  qsws_liq_roof(:,:)  !< roof latent heat flux (liquid incl. precipitation) (W m^-2)
        real(field_r), allocatable ::  qsws_road(:,:)      !< road latent heat flux (W m^-2)
        real(field_r), allocatable ::  qsws_roof(:,:)      !< roof latent heat flux (W m^-2)
        !
        !--    Liquid water coverages (storages).
        real(field_r), allocatable ::  c_liq_road(:,:)  !< liquid water coverage on road (.)
        real(field_r), allocatable ::  c_liq_roof(:,:)  !< liquid water coverage on roof (.)
        !
        !--    Diagnostic ground heat fluxes.
        real(field_r), allocatable ::  ghf_road(:,:)    !< road ground heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_roof(:,:)    !< roof indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_wall_a(:,:)  !< wall A indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_wall_b(:,:)  !< wall B indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_win_a(:,:)   !< window A indoor heat flux (W m^-2)
        real(field_r), allocatable ::  ghf_win_b(:,:)   !< window B indoor heat flux (W m^-2)
        !
        !--    Model internal radiation fluxes.
        real(field_r), allocatable ::  rad_lw_net_can(:,:)     !< net longwave radiative at canyon top (downwards) (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_road(:,:)    !< net longtwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_roof(:,:)    !< net longwave radiative flux on roof (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_urb(:,:)     !< urban aggegated net longwave radiative flux (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_wall_a(:,:)  !< net longwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_wall_b(:,:)  !< net longwave radiative flux wall B (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_win_a(:,:)   !< net longwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_lw_net_win_b(:,:)   !< net longwave radiative flux window B (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_road(:,:)     !< incoming shortwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_win_a(:,:)    !< incoming shortwave radiative flux on window A (W m^-2)
        real(field_r), allocatable ::  rad_sw_in_win_b(:,:)    !< incoming shortwave radiative flux on window B (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_road(:,:)    !< net shortwave radiative flux on road (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_roof(:,:)    !< net shortwave radiative flux on roof (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_urb(:,:)     !< urban aggegated net shortwave radiative flux (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_wall_a(:,:)  !< net shortwave radiative flux on wall A (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_wall_b(:,:)  !< net shortwave radiative flux on wall B (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_win_a(:,:)   !< net shortwave radiative flux on window A (W m^-2)
        real(field_r), allocatable ::  rad_sw_net_win_b(:,:)   !< net shortwave radiative flux on wall B (W m^-2)
        !
        !--    Surface layer model diagnostic variables.
        real(field_r), allocatable ::  ol_can(:,:)      !< canyon top Obukhov length (m)
        real(field_r), allocatable ::  ol_road(:,:)     !< road Obukhov length (m)
        real(field_r), allocatable ::  ol_roof(:,:)     !< roof Obukhov length (m)
        real(field_r), allocatable ::  pt_can(:,:)      !< street canyon virtual potential temperature (K)
        real(field_r), allocatable ::  rib_can(:,:)     !< canyon top bulk Richardson number (.)
        real(field_r), allocatable ::  rib_road(:,:)    !< road bulk Richardson number (.)
        real(field_r), allocatable ::  rib_roof(:,:)    !< roof bulk Richardson number (.)
        real(field_r), allocatable ::  us_can(:,:)      !< friction velocity for canyon resistance calculation (m s^-1)
        real(field_r), allocatable ::  uv_abs_can(:,:)  !< horizontal wind speed in street caynon at half-height (m s^-1)
        real(field_r), allocatable ::  uv_eff_can(:,:)  !< effective horizontal wind speed in street canyon at half-height (m s^-1)
        real(field_r), allocatable ::  vpt_can(:,:)     !< street canyon virtual potential temperature (K)
        !
        !--    Aerodynamic resistances for heat.
        real(field_r), allocatable ::  rah_can(:,:)     !< street canyon air aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_facade(:,:)  !< wall and window aerodynamic resistance for heat (combined) (s m^-1)
        real(field_r), allocatable ::  rah_road(:,:)    !< road aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_roof(:,:)    !< roof aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_wall_a(:,:)  !< wall A aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_wall_b(:,:)  !< wall B aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_win_a(:,:)   !< wall A aerodynamic resistance for heat (s m^-1)
        real(field_r), allocatable ::  rah_win_b(:,:)   !< wall B aerodynamic resistance for heat (s m^-1)
        !
        !--    Local friction velocities for roofs and roads.
        real(field_r), allocatable ::  us_road(:,:)  !< friction velocity for roads (m s^-1)
        real(field_r), allocatable ::  us_roof(:,:)  !< friction velocity for roofs (m s^-1)
        !
        !--    Diagnostic variables, defined at the first atmospheric grid level.
        real(field_r), allocatable ::  pt1(:,:)      !< potential temperature (K)
        real(field_r), allocatable ::  q1(:,:)       !< specific humidity (kg kg^-1)
        real(field_r), allocatable ::  us_urb(:,:)   !< friction velocity (m s^-1)
        real(field_r), allocatable ::  uv_abs1(:,:)  !< horizontal wind speed (m s^-1)
        real(field_r), allocatable ::  uv_eff1(:,:)  !< effective horizontal wind speed (m s^-1)
        real(field_r), allocatable ::  vpt1(:,:)     !< virtual potential temperature (K)
        !
        !--    Parameters for the whole urban tile.
        LOGICAL,  allocatable ::  anisotropic_canyon(:,:)  !< boolean flag to mark anisotropic canyon
        real(field_r), allocatable ::  f_bld(:,:)               !< fractional area occupied by buldings (plan area fraction) (.)
        real(field_r), allocatable ::  f_bld_frn(:,:)           !< frontal area fraction of buildings (.)
        real(field_r), allocatable ::  f_win(:,:)               !< window fraction (.)
        real(field_r), allocatable ::  h_bld(:,:)               !< building height (m)
        real(field_r), allocatable ::  hw_can(:,:)              !< canyon aspect ratio (.)
        real(field_r), allocatable ::  svf_road(:,:)            !< sky-view factor for road (.)
        real(field_r), allocatable ::  svf_wall(:,:)            !< sky-view-factor for walls (.)
        real(field_r), allocatable ::  theta_can(:,:)           !< canyon orientation / road direction in radians (.)
        real(field_r), allocatable ::  z0_urb(:,:)              !< aerodynamic roughness length of the urban surface (m)
        !
        !--    Material properties.
        real(field_r), allocatable ::  albedo_road(:,:)         !< albedo of the road (.)
        real(field_r), allocatable ::  albedo_roof(:,:)         !< albedo of the roof (.)
        real(field_r), allocatable ::  albedo_wall(:,:)         !< albedo of the wall (.)
        real(field_r), allocatable ::  albedo_wall_win(:,:)     !< weighted average of wall and window albedos for reflections (.)
        real(field_r), allocatable ::  albedo_win(:,:)          !< albedo of the window (.)
        real(field_r), allocatable ::  emiss_road(:,:)          !< emissivity of the road (.)
        real(field_r), allocatable ::  emiss_roof(:,:)          !< emissivity of the roof (.)
        real(field_r), allocatable ::  emiss_wall(:,:)          !< emissivity of the wall (.)
        real(field_r), allocatable ::  emiss_win(:,:)           !< emissivity of the window (.)
        real(field_r), allocatable ::  transmissivity_win(:,:)  !< transmissivity of the window layers (.)
        real(field_r), allocatable ::  z0_road(:,:)             !< aerodynamic roughness length for momentum for roads (m)
        real(field_r), allocatable ::  z0_roof(:,:)             !< aerodynamic roughness length for momentum of roofs (m)
        real(field_r), allocatable ::  z0_wall(:,:)             !< aerodynamic roughness length for walls and windows (m)
        real(field_r), allocatable ::  z0h_road(:,:)            !< aerodynamic roughness length for heat for roads (m)
        real(field_r), allocatable ::  z0h_roof(:,:)            !< aerodynamic roughness length for heat for roofs (m)

        real(field_r), allocatable ::  absorption_win(:,:,:)  !< fraction of absorbed shortwave radiation over glass sheet (.)
        real(field_r), allocatable ::  c_road(:,:,:)          !< total (specific c * layer depth) heat capacity of the road (J m^-2 K^-1)
        real(field_r), allocatable ::  c_roof(:,:,:)          !< total (specific c * layer depth) heat capacity of the roof (J m^-2 K^-1)
        real(field_r), allocatable ::  c_wall(:,:,:)          !< total (specific c * layer depth) heat capacity of the wall (J m^-2 K^-1)
        real(field_r), allocatable ::  c_win(:,:,:)           !< total (specific c * layer depth) heat heat capacity of the window (J m^-2 K^-1)
        real(field_r), allocatable ::  lambda_road(:,:,:)     !< thermal conductivity of the road (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_roof(:,:,:)     !< thermal conductivity of the roof (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_wall(:,:,:)     !< thermal conductivity of the wall (W m^-1 K^-1)
        real(field_r), allocatable ::  lambda_win(:,:,:)      !< effective thermal conductivity of the window (W m^-1 K^-1)
        !
        !--    Building indoor parameters.
        real(field_r), allocatable ::  t_indoor(:,:)  !< building indoor temperature (K)
        !
        !--    Soil parameters.
        real(field_r), allocatable ::  t_soil(:,:)  !< fixed soil top temperature (K)
        !
        !--    Pre-computed total layer conductivities.
        real(field_r), allocatable ::  conductivity_road(:,:,:)  !< total conductivity bewtween road layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_roof(:,:,:)  !< total conductivity between roof layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_wall(:,:,:)  !< total conductivity between wall layers (lambda_h / dz) (W m^-2 K^-1)
        real(field_r), allocatable ::  conductivity_win(:,:,:)   !< total conductivity between window layers (lambda_h / dz) (W m^-2 K^-1)
        !
        !--    Pre-computed variables and coefficients.
        real(field_r), allocatable ::  sw_ref_denom(:,:)      !< SW radiation reflection denominator (.)
        real(field_r), allocatable ::  uv_abs_can_coef(:,:)   !< coefficient for the canyon wind speed (.)
        real(field_r), allocatable ::  wall_hor_a_ratio(:,:)  !< wall-to-horizontal area ratio (unused) (.)
        real(field_r), allocatable ::  z_mo(:,:)              !< reference height for MOST for the atmosphere (m)
        real(field_r), allocatable ::  z_mo_can(:,:)          !< canyon reference height for MOST (canyon half-height) (m)

        real(field_r), allocatable ::  lw_road_coef(:,:,:)  !< LW radiation coefficients for roads (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_roof_coef(:,:,:)  !< LW radiation coefficients for roofs (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_wall_coef(:,:,:)  !< LW radiation coefficients for walls (W m^-2 K^-4 & . & W m^-2 K^-4 ...)
        real(field_r), allocatable ::  lw_win_coef(:,:,:)   !< LW radiation coefficients for walls (W m^-2 K^-4 & . & W m^-2 K^-4 ...)

    end type surf_slurb

    type(surf_slurb) :: slurb_tile
end module modslurbdata