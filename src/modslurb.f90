!> \file modslurb.f90
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2025 Delft University of Technology
!
module modslurb
    use netcdf
    use modprecision, only : field_r
    use ieee_arithmetic, only: ieee_is_nan
    use modslurbdata, only: slurb_tile
    use modglobal, only: g => grav, cp => cp, kappa => fkar, pi
    use modfields, only: rho_air_zw => rhobf, exnf

    implicit none
    save   
    
    real(field_r), allocatable ::  ln_z_z0_roof(:,:)   !< temporary array to store logarithm ZELFTODO (.)
    real(field_r), allocatable ::  ln_z_z0h_roof(:,:)  !< temporary array to store logarithm (.)
    real(field_r), allocatable ::  ln_z_z0_urb(:,:)    !< temporary array to store logarithm (.)
    real(field_r), allocatable ::  pt_surface(:,:)     !< temporary array to store weighted temperature (K)
    
    real(field_r), allocatable ::  ln_z_z0_road(:,:)   !< temporary array to store logarithm ZELFTODO (.)
    real(field_r), allocatable ::  ln_z_z0h_road(:,:)  !< temporary array to store logarithm (.)


    !
    !-- Target arrays for timelevel switching.
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  m_liq_road_1  !< target array for liquid water reservoir on roads (m^3 m^-2)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  m_liq_road_2  !< target array for liquid water reservoir on roads (m^3 m^-2)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  m_liq_roof_1  !< target array for liquid water reservoir on roofs (m^3 m^-2)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  m_liq_roof_2  !< target array for liquid water reservoir on roofs (m^3 m^-2)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  q_can_1       !< target array for canyon water mixing ratio (kg kg**-1)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  q_can_2       !< target array for canyon water mixing ratio (kg kg**-1)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_can_1  !< target array for canyon temperature used to change timelevels (K)
    REAL(field_r), DIMENSION(:,:), TARGET, ALLOCATABLE ::  t_can_2  !< target array for canyon temperature (K)

    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_road_1    !< target array for road temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_road_2    !< target array for road temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_roof_1    !< target array for roof temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_roof_2    !< target array for roof temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_wall_a_1  !< target array for wall A temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_wall_a_2  !< target array for wall A temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_wall_b_1  !< target array for wall B temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_wall_b_2  !< target array for wall B temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_win_a_1   !< target array for window A temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_win_a_2   !< target array for window A temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_win_b_1   !< target array for window B temperature (K)
    REAL(field_r), DIMENSION(:,:,:), TARGET, ALLOCATABLE ::  t_win_b_2   !< target array for window B temperature (K)
    !
    !-- Arrays for output temporal averaging.
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  albedo_urb_av         !< road liquid water coverage 
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  c_liq_road_av         !< road liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  c_liq_roof_av         !< roof liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  emiss_urb_av          !< road liquid water coverage
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_road_av           !< heat flux between the road bottom layer and soil
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_roof_av           !< heat flux between the roof bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_wall_a_av         !< heat flux between the wall a bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_wall_b_av         !< heat flux between the wall b bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_win_a_av          !< heat flux between the window a bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ghf_win_b_av          !< heat flux between the window b bottom layer and indoor air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  m_liq_road_av         !< liquid water reservoir on roads
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  m_liq_roof_av         !< liquid water reservoir on roofs
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_can_av             !< street canyon top obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_road_av            !< road obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_roof_av            !< roof obukhov length
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ol_urb_av             !< urban obukhov length for momentum flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_can_av             !< street canyon air potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_road_av            !< road surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_roof_av            !< roof surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_wall_a_av          !< wall a surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_wall_b_av          !< wall b surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_win_a_av           !< window a surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  pt_win_b_av           !< window b surface potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_can_av              !< street canyon water vapour mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_road_av             !< road surface mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  q_roof_av             !< roof surface mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qs_road_av            !< road surface saturation mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qs_roof_av            !< roof surface saturation mixing ratio
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_can_av           !< latent heat flux between the street canyon and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_external_av      !< latent heat flux external to the model (e.g. from industry)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_road_av          !< latent heat flux between the road and the street canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_roof_av          !< latent heat flux between the roof and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  qsws_urb_av           !< urban aggregated latent heat flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_road_av    !< road surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_roof_av    !< roof surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_urb_av     !< urban aggergated net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_wall_a_av  !< wall a surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_wall_b_av  !< wall b surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_win_a_av   !< window a surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_net_win_b_av   !< window b surface net longwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_road_av    !< road surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_roof_av    !< roof surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_urb_av     !< aggegated urban surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_wall_a_av  !< wall a surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_wall_b_av  !< wall b surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_win_a_av   !< window a surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_net_win_b_av   !< window b surface net shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_tr_win_a_av    !< window a surface transmitted shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_tr_win_b_av    !< window b surface transmitted shortwave radiative flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_can_av            !< street canyon aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_road_av           !< road aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_roof_av           !< roof aerodynamic resistance for heat
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_wall_a_av         !< wall A aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_wall_b_av         !< wall B aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_win_a_av          !< window A aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_win_b_av          !< window B aerodynamic resistance for heat (DOE-2)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rah_facade_av         !< wall and window aerodynamic resistance for heat (combined)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  ram_urb_av            !< urban aerodynamic resistance for momentum
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_can_av            !< street canyon top bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_road_av           !< road bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  rib_roof_av           !< roof bulk richardson number
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_can_av            !< sensible heat flux between the street canyon and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_external_av       !< sensible heat flux external to the model (e.g. from industry)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_road_av           !< sensible heat flux between the road and the street canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_roof_av           !< sensible heat flux between the roof and the atmosphere
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_traffic_av        !< sensible heat flux from traffic to the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_urb_av            !< urban aggregated sensible heat flux
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_wall_a_av         !< sensible heat flux between the wall a and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_wall_b_av         !< sensible heat flux between the wall b and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_win_a_av          !< sensible heat flux between the window a and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  shf_win_b_av          !< sensible heat flux between the window b and the canyon air
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_2m_urb_av           !< extrapolated 2-metre urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_c_urb_av            !< complete urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_can_av              !< street canyon air temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_h_urb_av            !< effective urban surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_rad_urb_av          !< effective urban surface radiative temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_road_av        !< road surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_roof_av        !< roof surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_wall_a_av      !< wall a surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_wall_b_av      !< wall b surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_win_a_av       !< window a surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  t_surf_win_b_av       !< window b surface temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_can_av             !< friction velocity for street canyons
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_road_av            !< friction velocity for roads
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_roof_av            !< friction velocity for roofs
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  us_urb_av             !< urban friction velocity
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  usws_urb_av           !< urban momentum flux (u-component)
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  uv_abs_can_av         !< street canyon wind speed
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  uv_eff_can_av         !< street canyon effective wind speed
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_can_av            !< street canyon air virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_road_av           !< road surface virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vpt_roof_av           !< roof surface virtual potential temperature
    ! REAL(field_r), DIMENSION(:,:), ALLOCATABLE ::  vsws_urb_av           !< urban momentum flux (v-component
    !
    !-- Arrays for output of unmodified LSM fluxes (2D).
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  shf_lsm_av   !< sensible heat flux from lsm surfaces
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  qsws_lsm_av  !< latent heat flux from lsm surfaces

    !
    !-- Arrays for output temporal averaging (2D).
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_road_av    !< road temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_roof_av    !< roof temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_wall_a_av  !< wall a temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_wall_b_av  !< wall b temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_win_a_av   !< window a temperature (all layers)
    ! REAL(field_r), DIMENSION(:,:,:), ALLOCATABLE ::  t_win_b_av   !< window b temperature (all layers)

    REAL(field_r) ::  dt_slurb = HUGE( 1.0_field_r )  !< maximum allowed timestep of SLUrb

    !
    !-- Model constants.
    REAL(field_r) ::  drho_l_lv  !< (rho_l * l_v)**-1
    REAL(field_r) ::  rho_lv     !< rho_surface * l_v

    !
    !-- Parameter defaults.
    REAL(field_r), PARAMETER ::  m_liq_max_road = 1.0E-3_field_r  !< maximum capacity of the liquid water reservoir on roads (i,j) (m^3 m^-2)
    REAL(field_r), PARAMETER ::  m_liq_max_roof = 1.0E-3_field_r  !< maximum capacity of the liquid water reservoir on roofs (i,j) (m^3 m^-2)
    REAL(field_r), PARAMETER ::  rah_max   = 1.0E6_field_r        !< maximum aerodynamic resistance for scalars (s m^-1)
    REAL(field_r), PARAMETER ::  rah_min   = 1.0_field_r          !< minimum aerodynamic resistance for scalars (s m^-1)
    REAL(field_r), PARAMETER ::  ram_min   = 1.0_field_r          !< minimum aerodynamic resistance for momentum (s m^-1) (TODOSELF)
    REAL(field_r), PARAMETER ::  urb_thres = 1.0E-2_field_r       !< minimum urban fraction to consider (1%) (.)
    REAL(field_r), PARAMETER ::  us_min    = 1.0E-8_field_r       !< minimum friction velocity (m s^-1)
    REAL(field_r), PARAMETER ::  zeta_min  = 1.0E-3_field_r       !< minimum stability parameter absolute value (neutral limit) (.)
    !
    !-- slurb_parameters namelist defaults.
    CHARACTER(LEN=20) ::  aero_roughness_heat = 'kanda'                !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  facade_resistance_parametrization = 'doe-2'  !< SLURrb namelist parameter
    CHARACTER(LEN=20) ::  street_canyon_wspeed_factor = 'surfex'       !< SLURrb namelist parameter

    integer ::  building_type = 2     !< SLURrb namelist parameter
    integer ::  n_layers_roads = 4    !< SLURrb namelist parameter
    integer ::  n_layers_roofs = 4    !< SLURrb namelist parameter
    integer ::  n_layers_walls = 4    !< SLURrb namelist parameter
    integer ::  n_layers_windows = 4  !< SLURrb namelist parameter
    integer ::  pavement_type = 2     !< SLURrb namelist parameter

    LOGICAL ::  anisotropic_street_canyons = .FALSE.  !< SLURrb namelist parameter
    LOGICAL ::  moist_physics = .TRUE.                !< SLURrb namelist parameter

    REAL(field_r) ::  building_frontal_area_fraction = -9999.0_field_r  !< SLURrb namelist parameter (.)
    REAL(field_r) ::  building_height = -9999.0_field_r                 !< SLURrb namelist parameter (m)
    REAL(field_r) ::  building_indoor_temperature =  -9999.0_field_r    !< SLURrb namelist parameter (K)
    REAL(field_r) ::  building_plan_area_fraction = -9999.0_field_r     !< SLURrb namelist parameter (.)
    REAl(field_r) ::  deep_soil_temperature = -9999.0_field_r           !< SLURrb namelist parameter (K)
    REAL(field_r) ::  qsws_external = 0.0_field_r                       !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  shf_external = 0.0_field_r                        !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  shf_traffic = 0.0_field_r                         !< SLURrb namelist parameter (W m^-2 s^-1)
    REAL(field_r) ::  street_canyon_aspect_ratio = -9999.0_field_r      !< SLURrb namelist parameter (.)
    REAL(field_r) ::  street_canyon_orientation = -9999.0_field_r       !< SLURrb namelist parameter (.)
    REAL(field_r) ::  urban_fraction = -9999.0_field_r                  !< SLURrb namelist parameter (.)
    REAL(field_r) ::  urban_roughness_length = -9999.0_field_r          !< SLURrb namelist parameter (m)
    REAL(field_r) ::  window_fraction = -9999.0_field_r                 !< SLURrb namelist parameter (.)


    REAL(field_r), PARAMETER ::  ol_max   = 1.0E6_field_r   !< allowed absolute maximum value Obukhov length (m)
    REAL(field_r), PARAMETER ::  ol_min   = 1.0E-6_field_r  !< allowed absolute minimum value Obukhov length (m)
    REAL(field_r), PARAMETER ::  ol_tol   = 1.0E-4_field_r  !< convergence limit for Obukhov length, relative tolerance (m)
    REAL(field_r), PARAMETER ::  rib_max  = 1.0E1_field_r   !< maximum bulk Richardson number (absolute value) (.)

    integer :: ibc_pt_b = 0 ! indicates dirichlet bc
    !SELFTODO NU 0 OM OL BETER TE BEREKENEN

    ! CHARACTER (LEN=20)   ::  bc_pt_b = 'dirichlet'                        !< namelist parameter
    ! CHARACTER (LEN=20)   ::  bc_pt_t = 'initial_gradient'                 !< namelist parameter


    !-- Internal logical switches for character-based namelist settings.
    !TODO ADD CHECKS
    LOGICAL ::  facade_rah_doe       = .TRUE.  !< facade resistance parameterization using DOE-2
    LOGICAL ::  facade_rah_kray      = .FALSE.  !< facade resistance parameterization using Krayenhoff&Voogt (2007)
    LOGICAL ::  facade_rah_rowley    = .FALSE.  !< facade resistance parameterization using Rowley (1932)
    LOGICAL ::  roughness_kanda      = .FALSE.  !< roughness parameterization of horizontal surfaces using Kanda et al. (2007)
    LOGICAL ::  uv_can_factor_kray   = .TRUE.  !< street canyon wind speed factor following Krayenhoff&Voogt (2007)
    LOGICAL ::  uv_can_factor_masson = .FALSE.  !< street canyon wind speed factor following Masson (2000)
    LOGICAL ::  uv_can_factor_surfex = .FALSE.  !< street canyon wind speed factor following the SURFEX model

    !-- Default subsurface layer configuration.
    INTEGER ::  nzt_wall  !< top of the wall model (outer surface)
    INTEGER ::  nzb_wall  !< bottom of the wall model (inside surface)
    INTEGER ::  nzt_win   !< top of the window model (outer surface)
    INTEGER ::  nzb_win   !< bottom of the window model (inside surface)
    INTEGER ::  nzt_roof  !< top of the roof model (outer surface)
    INTEGER ::  nzb_roof  !< bottom the roof model (inside surface)
    INTEGER ::  nzt_road  !< top of the road model
    INTEGER ::  nzb_road  !< bottom of the road model

    real(field_r), allocatable :: fraction_slurb(:,:) !< (.)

    !-- Default surface description.
    REAL(field_r), DIMENSION(0:45,1:6) ::  building_pars_slurb  !< building default parameters derived from USM
    REAL(field_r), DIMENSION(0:14,1:5) ::  pavement_pars_slurb  !< pavement default parameters derived from LSM

    real(field_r) :: output_fill_value = -99999.0_field_r
    logical :: data_output_raw = .false.
    logical :: spinup = .false.
    logical :: calc_t_2m = .true.
    logical :: calc_t_c = .true.
    logical :: calc_t_h = .true.

    REAL(field_r) ::  tsc(10) = (/ 1.0_field_r, 1.0_field_r, 0.0_field_r, 0.0_field_r, &    !< array used for controlling time-integration at different substeps
                 0.0_field_r, 0.0_field_r, 0.0_field_r, 0.0_field_r, 0.0_field_r, 0.0_field_r /)

    real(field_r) :: namemiss_road = -9999
    real(field_r) :: namemiss_wall = -9999
    real(field_r) :: namemiss_win = -9999
    real(field_r) :: namemiss_roof = -9999
    real(field_r) :: namalbedo_road = -9999
    real(field_r) :: namalbedo_wall = -9999
    real(field_r) :: namalbedo_win = -9999
    real(field_r) :: namalbedo_roof = -9999
    real(field_r) :: namlambda_road = -9999
    real(field_r) :: namlambda_wall = -9999
    real(field_r) :: namlambda_win = -9999
    real(field_r) :: namlambda_roof = -9999
    

contains


subroutine slurb_read_namelist
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror
    use modmpi,      only : myid, comm3d, mpierr, D_MPI_BCAST

    
    implicit none

    integer :: ierr

    real(field_r), allocatable :: canyon_orientation_tmp(:,:)

    ! Namelist definition
    namelist /NAMSLURB/ &
        urban_fraction, urban_roughness_length, building_plan_area_fraction, building_frontal_area_fraction, building_height, window_fraction,&
        street_canyon_aspect_ratio, building_type, pavement_type, anisotropic_street_canyons, street_canyon_orientation, deep_soil_temperature,building_indoor_temperature,shf_external,qsws_external


    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMSLURB, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMSLURB')
        write(6, NAMSLURB)
        close(ifnamopt)
    end if


    ! some possible defaults
    ! urban_fraction = 0.5
    ! urban_roughness_length = 0.5
    ! building_plan_area_fraction = 0.3
    ! building_frontal_area_fraction = 0.11
    ! building_height = 20.0
    ! window_fraction = 0.2
    ! street_canyon_aspect_ratio = 0.5
    ! building_type = 2
    ! pavement_type = 2
    ! anisotropic_street_canyons = .true.
    ! street_canyon_orientation = 0
    ! deep_soil_temperature = 288.0

    ! Broadcast namelist values to all MPI tasks
    call D_MPI_BCAST(urban_fraction,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(urban_roughness_length,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_plan_area_fraction,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_frontal_area_fraction,             1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_height,       1, 0, comm3d, mpierr)
    call D_MPI_BCAST(window_fraction,   1, 0, comm3d, mpierr)
    call D_MPI_BCAST(street_canyon_aspect_ratio,       1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_type, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(pavement_type, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(anisotropic_street_canyons, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(street_canyon_orientation, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(deep_soil_temperature, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(building_indoor_temperature, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(shf_external, 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(qsws_external, 1, 0, comm3d, mpierr)
end subroutine slurb_read_namelist

subroutine slurb_read_RADnamelist
    use modglobal,   only : ifnamopt, fname_options, checknamelisterror
    use modmpi,      only : myid, comm3d, mpierr, D_MPI_BCAST

    
    implicit none

    integer :: ierr

    real(field_r), allocatable :: canyon_orientation_tmp(:,:)

    ! Namelist definition
    namelist /NAMSLURBRAD/ &
        namemiss_road, namemiss_wall, namemiss_win, namemiss_roof, namalbedo_road, namalbedo_wall, namalbedo_win, namalbedo_roof, namlambda_road,namlambda_wall,namlambda_win,namlambda_roof


    ! Read namelist
    if (myid == 0) then
        open(ifnamopt, file=fname_options, status='old', iostat=ierr)
        read(ifnamopt, NAMSLURBRAD, iostat=ierr)
        call checknamelisterror(ierr, ifnamopt, 'NAMSLURBRAD')
        write(6, NAMSLURBRAD)
        close(ifnamopt)
    end if


    ! Broadcast namelist values to all MPI tasks
    call D_MPI_BCAST(namemiss_road,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namemiss_wall,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namemiss_win,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namemiss_roof,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namalbedo_road,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namalbedo_wall,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namalbedo_win,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namalbedo_roof,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namlambda_road,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namlambda_wall,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namlambda_win,  1, 0, comm3d, mpierr)
    call D_MPI_BCAST(namlambda_roof,  1, 0, comm3d, mpierr)
end subroutine slurb_read_RADnamelist
subroutine slurb_bulk_allocations
    use modglobal, only: i2, j2
    use, intrinsic :: IEEE_ARITHMETIC

    allocate(fraction_slurb(i2,j2))

    allocate(ln_z_z0_roof(i2,j2))
    allocate(ln_z_z0h_roof(i2,j2))
    allocate(ln_z_z0_urb(i2,j2))
    allocate(pt_surface(i2,j2))
    allocate(ln_z_z0_road(i2,j2))
    allocate(ln_z_z0h_road(i2,j2))
    !-- Bulk allocation
    ALLOCATE( slurb_tile%dz_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%dz_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%dz_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%dz_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%zw_win(nzt_win:nzb_win,i2,j2) )

    ALLOCATE( slurb_tile%t_c_urb(i2,j2) )
    ALLOCATE( slurb_tile%t_rad_urb(i2,j2) )
    ALLOCATE( slurb_tile%t_h_urb(i2,j2) )
    ALLOCATE( slurb_tile%t_2m_urb(i2,j2) )
    ALLOCATE( slurb_tile%shf_urb(i2,j2) )
    ALLOCATE( slurb_tile%qsws_urb(i2,j2) )
    ALLOCATE( slurb_tile%ol_urb(i2,j2) )
    ALLOCATE( slurb_tile%rib_urb(i2,j2) )
    ALLOCATE( slurb_tile%ram_urb(i2,j2) )
    ALLOCATE( slurb_tile%usws_urb(i2,j2) )
    ALLOCATE( slurb_tile%vsws_urb(i2,j2) )
    ALLOCATE( slurb_tile%thlskin(i2,j2) )
    ALLOCATE( slurb_tile%qtskin(i2,j2) )

    ALLOCATE( slurb_tile%albedo_urb(i2,j2) )
    ALLOCATE( slurb_tile%emiss_urb(i2,j2) )

    ALLOCATE( slurb_tile%t_indoor(i2,j2) )
    ALLOCATE( slurb_tile%t_soil(i2,j2) )

    ALLOCATE( slurb_tile%tt_can(i2,j2) )
    ALLOCATE( slurb_tile%tt_wall_a(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%tt_wall_b(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%tt_win_a(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%tt_win_b(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%tt_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%tt_road(nzt_road:nzb_road,i2,j2) )

    ALLOCATE( slurb_tile%pt_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%pt_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%pt_win_a(i2,j2) )
    ALLOCATE( slurb_tile%pt_win_b(i2,j2) )
    ALLOCATE( slurb_tile%pt_roof(i2,j2) )
    ALLOCATE( slurb_tile%pt_road(i2,j2) )

    ALLOCATE( slurb_tile%shf_can(i2,j2) )
    ALLOCATE( slurb_tile%shf_roof(i2,j2) )
    ALLOCATE( slurb_tile%shf_road(i2,j2) )
    ALLOCATE( slurb_tile%shf_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%shf_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%shf_win_a(i2,j2) )
    ALLOCATE( slurb_tile%shf_win_b(i2,j2) )

    ALLOCATE( slurb_tile%shf_external(i2,j2) )
    ALLOCATE( slurb_tile%shf_traffic(i2,j2) )

    ALLOCATE( slurb_tile%ghf_road(i2,j2) )
    ALLOCATE( slurb_tile%ghf_roof(i2,j2) )
    ALLOCATE( slurb_tile%ghf_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%ghf_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%ghf_win_a(i2,j2) )
    ALLOCATE( slurb_tile%ghf_win_b(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_in_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_out_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_out_urb(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_urb(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_urb(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_can(i2,j2) )

    ALLOCATE( slurb_tile%rad_lw_net_roof(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_roof(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_road(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_wall_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_wall_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_win_a(i2,j2) )
    ALLOCATE( slurb_tile%rad_lw_net_win_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_net_win_b(i2,j2) )
    ALLOCATE( slurb_tile%rad_sw_in_win_b(i2,j2) )

    ALLOCATE( slurb_tile%pt_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_eff_can(i2,j2) )
    ALLOCATE( slurb_tile%us_can(i2,j2) )
    ALLOCATE( slurb_tile%rib_can(i2,j2) )
    ALLOCATE( slurb_tile%ol_can(i2,j2) )

    ALLOCATE( slurb_tile%rib_roof(i2,j2) )
    ALLOCATE( slurb_tile%ol_roof(i2,j2) )
    ALLOCATE( slurb_tile%rib_road(i2,j2) )
    ALLOCATE( slurb_tile%ol_road(i2,j2) )

    ALLOCATE( slurb_tile%us_roof(i2,j2) )
    ALLOCATE( slurb_tile%us_road(i2,j2) )

    ALLOCATE( slurb_tile%hw_can(i2,j2) )
    ALLOCATE( slurb_tile%anisotropic_canyon(i2,j2) )
    ALLOCATE( slurb_tile%theta_can(i2,j2) )
    ALLOCATE( slurb_tile%h_bld(i2,j2) )
    ALLOCATE( slurb_tile%f_bld(i2,j2) )
    ALLOCATE( slurb_tile%f_bld_frn(i2,j2) )
    ALLOCATE( slurb_tile%f_win(i2,j2) )
    ALLOCATE( slurb_tile%svf_road(i2,j2) )
    ALLOCATE( slurb_tile%svf_wall(i2,j2) )
    ALLOCATE( slurb_tile%z0_urb(i2,j2) )

    ALLOCATE( slurb_tile%rah_roof(i2,j2) )
    ALLOCATE( slurb_tile%rah_road(i2,j2) )
    ALLOCATE( slurb_tile%rah_can(i2,j2) )

    IF ( facade_rah_doe )  THEN
       ALLOCATE( slurb_tile%rah_wall_a(i2,j2) )
       ALLOCATE( slurb_tile%rah_wall_b(i2,j2) )
       ALLOCATE( slurb_tile%rah_win_a(i2,j2) )
       ALLOCATE( slurb_tile%rah_win_b(i2,j2) )
    ELSE
       ALLOCATE( slurb_tile%rah_facade(i2,j2) )
    ENDIF

    ALLOCATE( slurb_tile%lambda_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%c_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%albedo_roof(i2,j2) )
    ALLOCATE( slurb_tile%emiss_roof(i2,j2) )
    ALLOCATE( slurb_tile%z0_roof(i2,j2) )
    ALLOCATE( slurb_tile%z0h_roof(i2,j2) )
    ALLOCATE( slurb_tile%lambda_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%c_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%albedo_wall(i2,j2) )
    ALLOCATE( slurb_tile%emiss_wall(i2,j2) )
    ALLOCATE( slurb_tile%z0_wall(i2,j2) )
    ALLOCATE( slurb_tile%lambda_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%c_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%albedo_wall_win(i2,j2) )
    ALLOCATE( slurb_tile%albedo_win(i2,j2) )
    ALLOCATE( slurb_tile%emiss_win(i2,j2) )
    ALLOCATE( slurb_tile%transmissivity_win(i2,j2) )
    ALLOCATE( slurb_tile%absorption_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%lambda_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%c_road(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( slurb_tile%albedo_road(i2,j2) )
    ALLOCATE( slurb_tile%emiss_road(i2,j2) )
    ALLOCATE( slurb_tile%z0_road(i2,j2) )
    ALLOCATE( slurb_tile%z0h_road(i2,j2) )

    ALLOCATE( slurb_tile%conductivity_roof(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_wall(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_win(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( slurb_tile%conductivity_road(nzt_road:nzb_road,i2,j2) )

    ALLOCATE( slurb_tile%z_mo(i2,j2) )
    ALLOCATE( slurb_tile%z_mo_can(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs_can_coef(i2,j2) )
    ALLOCATE( slurb_tile%wall_hor_a_ratio(i2,j2) )


    ALLOCATE( slurb_tile%lw_roof_coef(1:2,i2,j2) )
    ALLOCATE( slurb_tile%lw_road_coef(1:4,i2,j2) )
    ALLOCATE( slurb_tile%lw_wall_coef(1:6,i2,j2) )
    ALLOCATE( slurb_tile%lw_win_coef(1:6,i2,j2) )
    ALLOCATE( slurb_tile%sw_ref_denom(i2,j2) )

    ALLOCATE( slurb_tile%us_urb(i2,j2) )
    ALLOCATE( slurb_tile%uv_eff1(i2,j2) )
    ALLOCATE( slurb_tile%uv_abs1(i2,j2) )
    ALLOCATE( slurb_tile%pt1(i2,j2) )

    ALLOCATE( t_can_1(i2,j2) )
    ALLOCATE( t_can_2(i2,j2) )
    ALLOCATE( t_wall_a_1(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( t_wall_a_2(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( t_wall_b_1(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( t_wall_b_2(nzt_wall:nzb_wall,i2,j2) )
    ALLOCATE( t_win_a_1(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( t_win_a_2(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( t_win_b_1(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( t_win_b_2(nzt_win:nzb_win,i2,j2) )
    ALLOCATE( t_roof_1(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( t_roof_2(nzt_roof:nzb_roof,i2,j2) )
    ALLOCATE( t_road_1(nzt_road:nzb_road,i2,j2) )
    ALLOCATE( t_road_2(nzt_road:nzb_road,i2,j2) )

    IF ( moist_physics )  THEN
       ALLOCATE( slurb_tile%tq_can(i2,j2))
       ALLOCATE( slurb_tile%tm_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%tm_liq_road(i2,j2) )

       ALLOCATE( slurb_tile%vpt_roof(i2,j2) )
       ALLOCATE( slurb_tile%vpt_road(i2,j2) )

       ALLOCATE( slurb_tile%q_roof(i2,j2) )
       ALLOCATE( slurb_tile%q_road(i2,j2) )
       ALLOCATE( slurb_tile%qs_roof(i2,j2) )
       ALLOCATE( slurb_tile%qs_road(i2,j2) )

       ALLOCATE( slurb_tile%qsws_can(i2,j2) )
       ALLOCATE( slurb_tile%qsws_roof(i2,j2) )
       ALLOCATE( slurb_tile%qsws_road(i2,j2) )
       ALLOCATE( slurb_tile%qsws_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%qsws_liq_road(i2,j2) )

       ALLOCATE( slurb_tile%c_liq_roof(i2,j2) )
       ALLOCATE( slurb_tile%c_liq_road(i2,j2) )

       ALLOCATE( slurb_tile%vpt_can(i2,j2) )

       ALLOCATE( slurb_tile%q1(i2,j2) )
       ALLOCATE( slurb_tile%vpt1(i2,j2) )

       ALLOCATE( slurb_tile%qsws_external(i2,j2) )

       ALLOCATE( q_can_1(i2,j2) )
       ALLOCATE( q_can_2(i2,j2) )
       ALLOCATE( m_liq_roof_1(i2,j2) )
       ALLOCATE( m_liq_roof_2(i2,j2) )
       ALLOCATE( m_liq_road_1(i2,j2) )
       ALLOCATE( m_liq_road_2(i2,j2) )

    ENDIF

    ALLOCATE( slurb_tile%dt_max(i2,j2) )




    fraction_slurb(:,:) = ieee_value(fraction_slurb,ieee_signaling_nan)

    ln_z_z0_roof(:,:) = ieee_value(ln_z_z0_roof,ieee_signaling_nan)
    ln_z_z0h_roof(:,:) = ieee_value(ln_z_z0h_roof,ieee_signaling_nan)
    ln_z_z0_urb(:,:) = ieee_value(ln_z_z0_urb,ieee_signaling_nan)
    pt_surface(:,:) = ieee_value(pt_surface,ieee_signaling_nan)
    ln_z_z0_road(:,:) = ieee_value(ln_z_z0_road,ieee_signaling_nan)
    ln_z_z0h_road(:,:) = ieee_value(ln_z_z0h_road,ieee_signaling_nan)
    !-- Bulk allocation
    slurb_tile%dz_roof(:,:,:) = ieee_value(slurb_tile%dz_roof,ieee_signaling_nan)
    slurb_tile%dz_wall(:,:,:) = ieee_value(slurb_tile%dz_wall,ieee_signaling_nan)
    slurb_tile%dz_road(:,:,:) = ieee_value(slurb_tile%dz_road,ieee_signaling_nan)
    slurb_tile%dz_win(:,:,:) = ieee_value(slurb_tile%dz_win,ieee_signaling_nan)
    slurb_tile%zw_win(:,:,:) = ieee_value(slurb_tile%zw_win,ieee_signaling_nan)

    slurb_tile%t_c_urb(:,:) = ieee_value(slurb_tile%t_c_urb,ieee_signaling_nan)
    slurb_tile%t_rad_urb(:,:) = ieee_value(slurb_tile%t_rad_urb,ieee_signaling_nan)
    slurb_tile%t_h_urb(:,:) = ieee_value(slurb_tile%t_h_urb,ieee_signaling_nan)
    slurb_tile%t_2m_urb(:,:) = ieee_value(slurb_tile%t_2m_urb,ieee_signaling_nan)
    slurb_tile%shf_urb(:,:) = ieee_value(slurb_tile%shf_urb,ieee_signaling_nan)
    slurb_tile%qsws_urb(:,:) = ieee_value(slurb_tile%qsws_urb,ieee_signaling_nan)
    slurb_tile%ol_urb(:,:) = ieee_value(slurb_tile%ol_urb,ieee_signaling_nan)
    slurb_tile%rib_urb(:,:) = ieee_value(slurb_tile%rib_urb,ieee_signaling_nan)
    slurb_tile%ram_urb(:,:) = ieee_value(slurb_tile%ram_urb,ieee_signaling_nan)
    slurb_tile%usws_urb(:,:) = ieee_value(slurb_tile%usws_urb,ieee_signaling_nan)
    slurb_tile%vsws_urb(:,:) = ieee_value(slurb_tile%vsws_urb,ieee_signaling_nan)
    slurb_tile%thlskin(:,:) = ieee_value(slurb_tile%thlskin,ieee_signaling_nan)
    slurb_tile%qtskin(:,:) = ieee_value(slurb_tile%qtskin,ieee_signaling_nan)

    slurb_tile%albedo_urb(:,:) = ieee_value(slurb_tile%albedo_urb,ieee_signaling_nan)
    slurb_tile%emiss_urb(:,:) = ieee_value(slurb_tile%emiss_urb,ieee_signaling_nan)

    slurb_tile%t_indoor(:,:) = ieee_value(slurb_tile%t_indoor,ieee_signaling_nan)
    slurb_tile%t_soil(:,:) = ieee_value(slurb_tile%t_soil,ieee_signaling_nan)

    slurb_tile%tt_can(:,:) = ieee_value(slurb_tile%tt_can,ieee_signaling_nan)
    slurb_tile%tt_wall_a(:,:,:) = ieee_value(slurb_tile%tt_wall_a,ieee_signaling_nan)
    slurb_tile%tt_wall_b(:,:,:) = ieee_value(slurb_tile%tt_wall_b,ieee_signaling_nan)
    slurb_tile%tt_win_a(:,:,:) = ieee_value(slurb_tile%tt_win_a,ieee_signaling_nan)
    slurb_tile%tt_win_b(:,:,:) = ieee_value(slurb_tile%tt_win_b,ieee_signaling_nan)
    slurb_tile%tt_roof(:,:,:) = ieee_value(slurb_tile%tt_roof,ieee_signaling_nan)
    slurb_tile%tt_road(:,:,:) = ieee_value(slurb_tile%tt_road,ieee_signaling_nan)

    slurb_tile%pt_wall_a(:,:) = ieee_value(slurb_tile%pt_wall_a,ieee_signaling_nan)
    slurb_tile%pt_wall_b(:,:) = ieee_value(slurb_tile%pt_wall_b,ieee_signaling_nan)
    slurb_tile%pt_win_a(:,:) = ieee_value(slurb_tile%pt_win_a,ieee_signaling_nan)
    slurb_tile%pt_win_b(:,:) = ieee_value(slurb_tile%pt_win_b,ieee_signaling_nan)
    slurb_tile%pt_roof(:,:) = ieee_value(slurb_tile%pt_roof,ieee_signaling_nan)
    slurb_tile%pt_road(:,:) = ieee_value(slurb_tile%pt_road,ieee_signaling_nan)

    slurb_tile%shf_can(:,:) = ieee_value(slurb_tile%shf_can,ieee_signaling_nan)
    slurb_tile%shf_roof(:,:) = ieee_value(slurb_tile%shf_roof,ieee_signaling_nan)
    slurb_tile%shf_road(:,:) = ieee_value(slurb_tile%shf_road,ieee_signaling_nan)
    slurb_tile%shf_wall_a(:,:) = ieee_value(slurb_tile%shf_wall_a,ieee_signaling_nan)
    slurb_tile%shf_wall_b(:,:) = ieee_value(slurb_tile%shf_wall_b,ieee_signaling_nan)
    slurb_tile%shf_win_a(:,:) = ieee_value(slurb_tile%shf_win_a,ieee_signaling_nan)
    slurb_tile%shf_win_b(:,:) = ieee_value(slurb_tile%shf_win_b,ieee_signaling_nan)

    slurb_tile%shf_external(:,:) = ieee_value(slurb_tile%shf_external,ieee_signaling_nan)
    slurb_tile%shf_traffic(:,:) = ieee_value(slurb_tile%shf_traffic,ieee_signaling_nan)

    slurb_tile%ghf_road(:,:) = ieee_value(slurb_tile%ghf_road,ieee_signaling_nan)
    slurb_tile%ghf_roof(:,:) = ieee_value(slurb_tile%ghf_roof,ieee_signaling_nan)
    slurb_tile%ghf_wall_a(:,:) = ieee_value(slurb_tile%ghf_wall_a,ieee_signaling_nan)
    slurb_tile%ghf_wall_b(:,:) = ieee_value(slurb_tile%ghf_wall_b,ieee_signaling_nan)
    slurb_tile%ghf_win_a(:,:) = ieee_value(slurb_tile%ghf_win_a,ieee_signaling_nan)
    slurb_tile%ghf_win_b(:,:) = ieee_value(slurb_tile%ghf_win_b,ieee_signaling_nan)

    slurb_tile%rad_lw_in_urb(:,:) = ieee_value(slurb_tile%rad_lw_in_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_in_urb(:,:) = ieee_value(slurb_tile%rad_sw_in_urb,ieee_signaling_nan)
    slurb_tile%rad_lw_out_urb(:,:) = ieee_value(slurb_tile%rad_lw_out_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_out_urb(:,:) = ieee_value(slurb_tile%rad_sw_out_urb,ieee_signaling_nan)

    slurb_tile%rad_lw_net_urb(:,:) = ieee_value(slurb_tile%rad_lw_net_urb,ieee_signaling_nan)
    slurb_tile%rad_sw_net_urb(:,:) = ieee_value(slurb_tile%rad_sw_net_urb,ieee_signaling_nan)

    slurb_tile%rad_lw_net_can(:,:) = ieee_value(slurb_tile%rad_lw_net_can,ieee_signaling_nan)

    slurb_tile%rad_lw_net_roof(:,:) = ieee_value(slurb_tile%rad_lw_net_roof,ieee_signaling_nan)
    slurb_tile%rad_sw_net_roof(:,:) = ieee_value(slurb_tile%rad_sw_net_roof,ieee_signaling_nan)
    slurb_tile%rad_lw_net_road(:,:) = ieee_value(slurb_tile%rad_lw_net_road,ieee_signaling_nan)
    slurb_tile%rad_sw_net_road(:,:) = ieee_value(slurb_tile%rad_sw_net_road,ieee_signaling_nan)
    slurb_tile%rad_sw_in_road(:,:) = ieee_value(slurb_tile%rad_sw_in_road,ieee_signaling_nan)
    slurb_tile%rad_lw_net_wall_a(:,:) = ieee_value(slurb_tile%rad_lw_net_wall_a,ieee_signaling_nan)
    slurb_tile%rad_sw_net_wall_a(:,:) = ieee_value(slurb_tile%rad_sw_net_wall_a,ieee_signaling_nan)
    slurb_tile%rad_lw_net_wall_b(:,:) = ieee_value(slurb_tile%rad_lw_net_wall_b,ieee_signaling_nan)
    slurb_tile%rad_sw_net_wall_b(:,:) = ieee_value(slurb_tile%rad_sw_net_wall_b,ieee_signaling_nan)
    slurb_tile%rad_lw_net_win_a(:,:) = ieee_value(slurb_tile%rad_lw_net_win_a,ieee_signaling_nan)
    slurb_tile%rad_sw_net_win_a(:,:) = ieee_value(slurb_tile%rad_sw_net_win_a,ieee_signaling_nan)
    slurb_tile%rad_sw_in_win_a(:,:) = ieee_value(slurb_tile%rad_sw_in_win_a,ieee_signaling_nan)
    slurb_tile%rad_lw_net_win_b(:,:) = ieee_value(slurb_tile%rad_lw_net_win_b,ieee_signaling_nan)
    slurb_tile%rad_sw_net_win_b(:,:) = ieee_value(slurb_tile%rad_sw_net_win_b,ieee_signaling_nan)
    slurb_tile%rad_sw_in_win_b(:,:) = ieee_value(slurb_tile%rad_sw_in_win_b,ieee_signaling_nan)

    slurb_tile%pt_can(:,:) = ieee_value(slurb_tile%pt_can,ieee_signaling_nan)
    slurb_tile%uv_abs_can(:,:) = ieee_value(slurb_tile%uv_abs_can,ieee_signaling_nan)
    slurb_tile%uv_eff_can(:,:) = ieee_value(slurb_tile%uv_eff_can,ieee_signaling_nan)
    slurb_tile%us_can(:,:) = ieee_value(slurb_tile%us_can,ieee_signaling_nan)
    slurb_tile%rib_can(:,:) = ieee_value(slurb_tile%rib_can,ieee_signaling_nan)
    slurb_tile%ol_can(:,:) = ieee_value(slurb_tile%ol_can,ieee_signaling_nan)

    slurb_tile%rib_roof(:,:) = ieee_value(slurb_tile%rib_roof,ieee_signaling_nan)
    slurb_tile%ol_roof(:,:) = ieee_value(slurb_tile%ol_roof,ieee_signaling_nan)
    slurb_tile%rib_road(:,:) = ieee_value(slurb_tile%rib_road,ieee_signaling_nan)
    slurb_tile%ol_road(:,:) = ieee_value(slurb_tile%ol_road,ieee_signaling_nan)

    slurb_tile%us_roof(:,:) = ieee_value(slurb_tile%us_roof,ieee_signaling_nan)
    slurb_tile%us_road(:,:) = ieee_value(slurb_tile%us_road,ieee_signaling_nan)

    slurb_tile%hw_can(:,:) = ieee_value(slurb_tile%hw_can,ieee_signaling_nan)
    ! slurb_tile%anisotropic_canyon(:,:) = ieee_value(slurb_tile%anisotropic_canyon,ieee_signaling_nan)
    slurb_tile%theta_can(:,:) = ieee_value(slurb_tile%theta_can,ieee_signaling_nan)
    slurb_tile%h_bld(:,:) = ieee_value(slurb_tile%h_bld,ieee_signaling_nan)
    slurb_tile%f_bld(:,:) = ieee_value(slurb_tile%f_bld,ieee_signaling_nan)
    slurb_tile%f_bld_frn(:,:) = ieee_value(slurb_tile%f_bld_frn,ieee_signaling_nan)
    slurb_tile%f_win(:,:) = ieee_value(slurb_tile%f_win,ieee_signaling_nan)
    slurb_tile%svf_road(:,:) = ieee_value(slurb_tile%svf_road,ieee_signaling_nan)
    slurb_tile%svf_wall(:,:) = ieee_value(slurb_tile%svf_wall,ieee_signaling_nan)
    slurb_tile%z0_urb(:,:) = ieee_value(slurb_tile%z0_urb,ieee_signaling_nan)

    slurb_tile%rah_roof(:,:) = ieee_value(slurb_tile%rah_roof,ieee_signaling_nan)
    slurb_tile%rah_road(:,:) = ieee_value(slurb_tile%rah_road,ieee_signaling_nan)
    slurb_tile%rah_can(:,:) = ieee_value(slurb_tile%rah_can,ieee_signaling_nan)

    IF ( facade_rah_doe )  THEN
       slurb_tile%rah_wall_a(:,:) = ieee_value(slurb_tile%rah_wall_a,ieee_signaling_nan)
       slurb_tile%rah_wall_b(:,:) = ieee_value(slurb_tile%rah_wall_b,ieee_signaling_nan)
       slurb_tile%rah_win_a(:,:) = ieee_value(slurb_tile%rah_win_a,ieee_signaling_nan)
       slurb_tile%rah_win_b(:,:) = ieee_value(slurb_tile%rah_win_b,ieee_signaling_nan)
    ELSE
       slurb_tile%rah_facade(:,:) = ieee_value(slurb_tile%rah_facade,ieee_signaling_nan)
    ENDIF

    slurb_tile%lambda_roof(:,:,:) = ieee_value(slurb_tile%lambda_roof,ieee_signaling_nan)
    slurb_tile%c_roof(:,:,:) = ieee_value(slurb_tile%c_roof,ieee_signaling_nan)
    slurb_tile%albedo_roof(:,:) = ieee_value(slurb_tile%albedo_roof,ieee_signaling_nan)
    slurb_tile%emiss_roof(:,:) = ieee_value(slurb_tile%emiss_roof,ieee_signaling_nan)
    slurb_tile%z0_roof(:,:) = ieee_value(slurb_tile%z0_roof,ieee_signaling_nan)
    slurb_tile%z0h_roof(:,:) = ieee_value(slurb_tile%z0h_roof,ieee_signaling_nan)
    slurb_tile%lambda_wall(:,:,:) = ieee_value(slurb_tile%lambda_wall,ieee_signaling_nan)
    slurb_tile%c_wall(:,:,:) = ieee_value(slurb_tile%c_wall,ieee_signaling_nan)
    slurb_tile%albedo_wall(:,:) = ieee_value(slurb_tile%albedo_wall,ieee_signaling_nan)
    slurb_tile%emiss_wall(:,:) = ieee_value(slurb_tile%emiss_wall,ieee_signaling_nan)
    slurb_tile%z0_wall(:,:) = ieee_value(slurb_tile%z0_wall,ieee_signaling_nan)
    slurb_tile%lambda_win(:,:,:) = ieee_value(slurb_tile%lambda_win,ieee_signaling_nan)
    slurb_tile%c_win(:,:,:) = ieee_value(slurb_tile%c_win,ieee_signaling_nan)
    slurb_tile%albedo_wall_win(:,:) = ieee_value(slurb_tile%albedo_wall_win,ieee_signaling_nan)
    slurb_tile%albedo_win(:,:) = ieee_value(slurb_tile%albedo_win,ieee_signaling_nan)
    slurb_tile%emiss_win(:,:) = ieee_value(slurb_tile%emiss_win,ieee_signaling_nan)
    slurb_tile%transmissivity_win(:,:) = ieee_value(slurb_tile%transmissivity_win,ieee_signaling_nan)
    slurb_tile%absorption_win(:,:,:) = ieee_value(slurb_tile%absorption_win,ieee_signaling_nan)
    slurb_tile%lambda_road(:,:,:) = ieee_value(slurb_tile%lambda_road,ieee_signaling_nan)
    slurb_tile%c_road(:,:,:) = ieee_value(slurb_tile%c_road,ieee_signaling_nan)
    slurb_tile%albedo_road(:,:) = ieee_value(slurb_tile%albedo_road,ieee_signaling_nan)
    slurb_tile%emiss_road(:,:) = ieee_value(slurb_tile%emiss_road,ieee_signaling_nan)
    slurb_tile%z0_road(:,:) = ieee_value(slurb_tile%z0_road,ieee_signaling_nan)
    slurb_tile%z0h_road(:,:) = ieee_value(slurb_tile%z0h_road,ieee_signaling_nan)

    slurb_tile%conductivity_roof(:,:,:) = ieee_value(slurb_tile%conductivity_roof,ieee_signaling_nan)
    slurb_tile%conductivity_wall(:,:,:) = ieee_value(slurb_tile%conductivity_wall,ieee_signaling_nan)
    slurb_tile%conductivity_win(:,:,:) = ieee_value(slurb_tile%conductivity_win,ieee_signaling_nan)
    slurb_tile%conductivity_road(:,:,:) = ieee_value(slurb_tile%conductivity_road,ieee_signaling_nan)

    slurb_tile%z_mo(:,:) = ieee_value(slurb_tile%z_mo,ieee_signaling_nan)
    slurb_tile%z_mo_can(:,:) = ieee_value(slurb_tile%z_mo_can,ieee_signaling_nan)
    slurb_tile%uv_abs_can_coef(:,:) = ieee_value(slurb_tile%uv_abs_can_coef,ieee_signaling_nan)
    slurb_tile%wall_hor_a_ratio(:,:) = ieee_value(slurb_tile%wall_hor_a_ratio,ieee_signaling_nan)


    slurb_tile%lw_roof_coef(:,:,:) = ieee_value(slurb_tile%lw_roof_coef,ieee_signaling_nan)
    slurb_tile%lw_road_coef(:,:,:) = ieee_value(slurb_tile%lw_road_coef,ieee_signaling_nan)
    slurb_tile%lw_wall_coef(:,:,:) = ieee_value(slurb_tile%lw_wall_coef,ieee_signaling_nan)
    slurb_tile%lw_win_coef(:,:,:) = ieee_value(slurb_tile%lw_win_coef,ieee_signaling_nan)
    slurb_tile%sw_ref_denom(:,:) = ieee_value(slurb_tile%sw_ref_denom,ieee_signaling_nan)

    slurb_tile%us_urb(:,:) = ieee_value(slurb_tile%us_urb,ieee_signaling_nan)
    slurb_tile%uv_eff1(:,:) = ieee_value(slurb_tile%uv_eff1,ieee_signaling_nan)
    slurb_tile%uv_abs1(:,:) = ieee_value(slurb_tile%uv_abs1,ieee_signaling_nan)
    slurb_tile%pt1(:,:) = ieee_value(slurb_tile%pt1,ieee_signaling_nan)

    t_can_1(:,:) = ieee_value(t_can_1,ieee_signaling_nan)
    t_can_2(:,:) = ieee_value(t_can_2,ieee_signaling_nan)
    t_wall_a_1(:,:,:) = ieee_value(t_wall_a_1,ieee_signaling_nan)
    t_wall_a_2(:,:,:) = ieee_value(t_wall_a_2,ieee_signaling_nan)
    t_wall_b_1(:,:,:) = ieee_value(t_wall_b_1,ieee_signaling_nan)
    t_wall_b_2(:,:,:) = ieee_value(t_wall_b_2,ieee_signaling_nan)
    t_win_a_1(:,:,:) = ieee_value(t_win_a_1,ieee_signaling_nan)
    t_win_a_2(:,:,:) = ieee_value(t_win_a_2,ieee_signaling_nan)
    t_win_b_1(:,:,:) = ieee_value(t_win_b_1,ieee_signaling_nan)
    t_win_b_2(:,:,:) = ieee_value(t_win_b_2,ieee_signaling_nan)
    t_roof_1(:,:,:) = ieee_value(t_roof_1,ieee_signaling_nan)
    t_roof_2(:,:,:) = ieee_value(t_roof_2,ieee_signaling_nan)
    t_road_1(:,:,:) = ieee_value(t_road_1,ieee_signaling_nan)
    t_road_2(:,:,:) = ieee_value(t_road_2,ieee_signaling_nan)

    IF ( moist_physics )  THEN
       slurb_tile%tq_can(:,:) = ieee_value(slurb_tile%tq_can,ieee_signaling_nan)
       slurb_tile%tm_liq_roof(:,:) = ieee_value(slurb_tile%tm_liq_roof,ieee_signaling_nan)
       slurb_tile%tm_liq_road(:,:) = ieee_value(slurb_tile%tm_liq_road,ieee_signaling_nan)

       slurb_tile%vpt_roof(:,:) = ieee_value(slurb_tile%vpt_roof,ieee_signaling_nan)
       slurb_tile%vpt_road(:,:) = ieee_value(slurb_tile%vpt_road,ieee_signaling_nan)

       slurb_tile%q_roof(:,:) = ieee_value(slurb_tile%q_roof,ieee_signaling_nan)
       slurb_tile%q_road(:,:) = ieee_value(slurb_tile%q_road,ieee_signaling_nan)
       slurb_tile%qs_roof(:,:) = ieee_value(slurb_tile%qs_roof,ieee_signaling_nan)
       slurb_tile%qs_road(:,:) = ieee_value(slurb_tile%qs_road,ieee_signaling_nan)

       slurb_tile%qsws_can(:,:) = ieee_value(slurb_tile%qsws_can,ieee_signaling_nan)
       slurb_tile%qsws_roof(:,:) = ieee_value(slurb_tile%qsws_roof,ieee_signaling_nan)
       slurb_tile%qsws_road(:,:) = ieee_value(slurb_tile%qsws_road,ieee_signaling_nan)
       slurb_tile%qsws_liq_roof(:,:) = ieee_value(slurb_tile%qsws_liq_roof,ieee_signaling_nan)
       slurb_tile%qsws_liq_road(:,:) = ieee_value(slurb_tile%qsws_liq_road,ieee_signaling_nan)

       slurb_tile%c_liq_roof(:,:) = ieee_value(slurb_tile%c_liq_roof,ieee_signaling_nan)
       slurb_tile%c_liq_road(:,:) = ieee_value(slurb_tile%c_liq_road,ieee_signaling_nan)

       slurb_tile%vpt_can(:,:) = ieee_value(slurb_tile%vpt_can,ieee_signaling_nan)

       slurb_tile%q1(:,:) = ieee_value(slurb_tile%q1,ieee_signaling_nan)
       slurb_tile%vpt1(:,:) = ieee_value(slurb_tile%vpt1,ieee_signaling_nan)

       slurb_tile%qsws_external(:,:) = ieee_value(slurb_tile%qsws_external,ieee_signaling_nan)

       q_can_1(:,:) = ieee_value(q_can_1,ieee_signaling_nan)
       q_can_2(:,:) = ieee_value(q_can_2,ieee_signaling_nan)
       m_liq_roof_1(:,:) = ieee_value(m_liq_roof_1,ieee_signaling_nan)
       m_liq_roof_2(:,:) = ieee_value(m_liq_roof_2,ieee_signaling_nan)
       m_liq_road_1(:,:) = ieee_value(m_liq_road_1,ieee_signaling_nan)
       m_liq_road_2(:,:) = ieee_value(m_liq_road_2,ieee_signaling_nan)

    ENDIF

    slurb_tile%dt_max(:,:) = ieee_value(slurb_tile%dt_max,ieee_signaling_nan)
end subroutine slurb_bulk_allocations
subroutine initslurb
    use modglobal,   only : i1, j1, i2, j2, zf
    use modlsmdata, only : ilu, nlu, tile
    implicit none

    integer :: i, j

    real(field_r), allocatable :: canyon_orientation_tmp(:,:)


    call slurb_read_namelist



    !-- Initialize bounds for subsurface layers.
    nzt_road = 1
    nzb_road = n_layers_roads
    nzt_roof = 1
    nzb_roof = n_layers_roofs
    nzt_wall = 1
    nzb_wall = n_layers_walls
    nzt_win  = 1
    nzb_win  = n_layers_windows

    call slurb_bulk_allocations

    ! do j=2,j1
    !   do i=2,i1
    !    fraction_slurb(i,j) = 1
    !   enddo
    ! enddo
    do ilu=1,nlu
        if (tile(ilu)%lushort == "slb") then
            do j=2,j1
                do i=2,i1
                    fraction_slurb(i,j) = tile(ilu)%frac(i,j) ! (.)
                enddo
            enddo
        endif
    end do




    !
    !-- Process variables related to urban form.
    !
    !-- Internally, f_bld refers to building plan area fraction of the urban surface. However, it is
    !-- more common to report the building plan area fraction as a fraction of the total surface, e.g.
    !-- in the case of LCZs. We want the user input correspond to the latter, thus the scaling.
    slurb_tile%f_bld(:,:) = -9999.0_field_r
    ! CALL get_grid_variable_1d_real( 'building_plan_area_fraction', slurb_tile%f_bld,                     &
    !                                 building_plan_area_fraction )
    ! CALL check_grid_variable_1d_real( 'building_plan_area_fraction', slurb_tile%f_bld,                   &
    !                                   TINY( 1.0_field_r ), 0.99_field_r)
    slurb_tile%f_bld(:,:) = 0.5_field_r !TODOSELF deze is nu 0.5 voor nu
    do j=2,j1
      do i=2,i1
        ! IF ( slurb_tile%f_bld(i,j) > fr_urb(j,i) ) THEN
        !     WRITE( message_string, * ) 'building_plan_area_fraction = ', slurb_tile%f_bld(i,j),              &
        !                                 ' is higher than urban_fraction = ', fr_urb(j,i),             &
        !                                 ' for grid cell (j,i) = ', j, i, '.'
        !     CALL message( 'slurb_init', 'SLU1020', 2, 2, 0, 6, 0 )
        ! ENDIF
        if ( fraction_slurb(i,j) /= 0) then
            slurb_tile%f_bld(i, j) = slurb_tile%f_bld(i, j) / fraction_slurb(i, j)
        endif
        enddo
    enddo

    slurb_tile%f_bld_frn(:,:) = -9999.0_field_r
    slurb_tile%f_bld_frn(:,:) = 0.2_field_r !TODOSELF even vastgezet
    ! CALL get_grid_variable_1d_real( 'building_frontal_area_fraction', slurb_tile%f_bld_frn,              &
    !                                 building_frontal_area_fraction )
    ! CALL check_grid_variable_1d_real( 'building_frontal_area_fraction', slurb_tile%f_bld_frn,            &
    !                                   0.0_field_r, HUGE( 1.0_field_r ) )

    slurb_tile%h_bld(:,:) = -9999.0_field_r
    slurb_tile%h_bld(:,:) = building_height !TODOSELF even vastgezet (m)
    ! CALL get_grid_variable_1d_real( 'building_height', slurb_tile%h_bld, building_height )
    ! CALL check_grid_variable_1d_real( 'building_height', slurb_tile%h_bld, 0.0_field_r, 1000.0_field_r )

    !
    !-- Urban surface and street canyon MOST heights.
    do j=2,j1
      do i=2,i1
       slurb_tile%z_mo(i,j) = 0.5_field_r * (zf(2) - zf(1)) ! (m)
    !    slurb_tile%z_mo(i,j) = 0.5_field_r *  dzw(topo_top_ind(j,i,0)+1)
       slurb_tile%z_mo_can(i,j) = 0.5_field_r * building_height ! (m)
    !    slurb_tile%z_mo_can(i,j) = 0.5_field_r * slurb_tile%h_bld(i,j)
      enddo
    enddo
    !
    !-- Process canyon direction information if anisotropic street canyons are enabled.
    IF ( anisotropic_street_canyons )  THEN
       ALLOCATE( canyon_orientation_tmp(i2,j2) )
       canyon_orientation_tmp(:,:) = -9999.0_field_r
        IF ( street_canyon_orientation /= -9999.0_field_r  )  THEN
            canyon_orientation_tmp(:,:) = street_canyon_orientation
        ENDIF
    !    CALL get_grid_variable_1d_real( 'street_canyon_orientation', canyon_orientation_tmp,        &
    !                                    street_canyon_orientation )
    !
    !--    In order to make it possible to have a mix of isotropic and anisotropic tiles, use
    !--    anisotropic canyons only in the case a canyon orientation has been given either for all tiles
    !--    in the namelist or per-patch basis in input file.
        do j=2,j1
            do i=2,i1
                IF ( canyon_orientation_tmp(i,j) /= -9999.0_field_r )  THEN
                    slurb_tile%anisotropic_canyon(i,j) = .TRUE.
                    slurb_tile%theta_can(i,j) = canyon_orientation_tmp(i,j) * ( pi / 180.0_field_r )
                ELSE
                    slurb_tile%anisotropic_canyon(i,j) = .FALSE.
                    slurb_tile%theta_can(i,j) = -9999.0_field_r
                ENDIF
            enddo
        enddo
       DEALLOCATE( canyon_orientation_tmp )
    ELSE
        do j=2,j1
            do i=2,i1
                slurb_tile%anisotropic_canyon(i,j) = .FALSE.
                slurb_tile%theta_can(i,j) = -9999.0_field_r
            enddo
        enddo
    ENDIF

    ! CALL get_grid_variable_1d_real( 'street_canyon_aspect_ratio', slurb_tile%hw_can,                     &
    !                                 street_canyon_aspect_ratio )
    ! CALL check_grid_variable_1d_real( 'street_canyon_aspect_ratio', slurb_tile%hw_can,                   &
    !                                   TINY( 1.0_field_r ), HUGE( 0.0_field_r ) )

    ! CALL get_grid_variable_1d_real( 'z0_urb', slurb_tile%z0_urb, urban_roughness_length )
    ! CALL check_grid_variable_1d_real( 'z0_urb', slurb_tile%z0_urb, TINY( 1.0_field_r ), MINVAL( slurb_tile%z_mo ) )

    do j=2,j1
        do i=2,i1
            slurb_tile%z0_urb(i,j) = urban_roughness_length ! (m)
        enddo
    enddo
    ! do j=2,j1
    !     do i=2,i1
    !         slurb_tile%dt_max(i,j) = HUGE( 1.0_field_r )
    !     enddo
    ! enddo
    slurb_tile%dt_max(:,:) = HUGE( 1.0_field_r ) ! (s)

    call slurb_swap_timelevel(0)

    call process_surface_parameters

    call precompute_latent_variables

    call init_slurb_variables

    ! call radiation
    
    



end subroutine initslurb

subroutine exitslurb
    implicit none

    DEALLOCATE(fraction_slurb)

    DEALLOCATE(ln_z_z0_roof)
    DEALLOCATE(ln_z_z0h_roof)
    DEALLOCATE(ln_z_z0_urb)
    DEALLOCATE(pt_surface)
    DEALLOCATE(ln_z_z0_road)
    DEALLOCATE(ln_z_z0h_road)
    !-- Bulk allocation
    DEALLOCATE( slurb_tile%dz_roof)
    DEALLOCATE( slurb_tile%dz_wall)
    DEALLOCATE( slurb_tile%dz_road)
    DEALLOCATE( slurb_tile%dz_win)
    DEALLOCATE( slurb_tile%zw_win)

    DEALLOCATE( slurb_tile%t_c_urb)
    DEALLOCATE( slurb_tile%t_rad_urb)
    DEALLOCATE( slurb_tile%t_h_urb)
    DEALLOCATE( slurb_tile%t_2m_urb)
    DEALLOCATE( slurb_tile%shf_urb)
    DEALLOCATE( slurb_tile%qsws_urb)
    DEALLOCATE( slurb_tile%ol_urb)
    DEALLOCATE( slurb_tile%rib_urb)
    DEALLOCATE( slurb_tile%ram_urb)
    DEALLOCATE( slurb_tile%usws_urb)
    DEALLOCATE( slurb_tile%vsws_urb)
    DEALLOCATE( slurb_tile%thlskin)
    DEALLOCATE( slurb_tile%qtskin)

    DEALLOCATE( slurb_tile%albedo_urb)
    DEALLOCATE( slurb_tile%emiss_urb)

    DEALLOCATE( slurb_tile%t_indoor)
    DEALLOCATE( slurb_tile%t_soil)

    DEALLOCATE( slurb_tile%tt_can)
    DEALLOCATE( slurb_tile%tt_wall_a)
    DEALLOCATE( slurb_tile%tt_wall_b)
    DEALLOCATE( slurb_tile%tt_win_a)
    DEALLOCATE( slurb_tile%tt_win_b)
    DEALLOCATE( slurb_tile%tt_roof)
    DEALLOCATE( slurb_tile%tt_road)

    DEALLOCATE( slurb_tile%pt_wall_a)
    DEALLOCATE( slurb_tile%pt_wall_b)
    DEALLOCATE( slurb_tile%pt_win_a)
    DEALLOCATE( slurb_tile%pt_win_b)
    DEALLOCATE( slurb_tile%pt_roof)
    DEALLOCATE( slurb_tile%pt_road)

    DEALLOCATE( slurb_tile%shf_can)
    DEALLOCATE( slurb_tile%shf_roof)
    DEALLOCATE( slurb_tile%shf_road)
    DEALLOCATE( slurb_tile%shf_wall_a)
    DEALLOCATE( slurb_tile%shf_wall_b)
    DEALLOCATE( slurb_tile%shf_win_a)
    DEALLOCATE( slurb_tile%shf_win_b)

    DEALLOCATE( slurb_tile%shf_external)
    DEALLOCATE( slurb_tile%shf_traffic)

    DEALLOCATE( slurb_tile%ghf_road)
    DEALLOCATE( slurb_tile%ghf_roof)
    DEALLOCATE( slurb_tile%ghf_wall_a)
    DEALLOCATE( slurb_tile%ghf_wall_b)
    DEALLOCATE( slurb_tile%ghf_win_a)
    DEALLOCATE( slurb_tile%ghf_win_b)

    DEALLOCATE( slurb_tile%rad_lw_in_urb)
    DEALLOCATE( slurb_tile%rad_sw_in_urb)
    DEALLOCATE( slurb_tile%rad_lw_out_urb)
    DEALLOCATE( slurb_tile%rad_sw_out_urb)

    DEALLOCATE( slurb_tile%rad_lw_net_urb)
    DEALLOCATE( slurb_tile%rad_sw_net_urb)

    DEALLOCATE( slurb_tile%rad_lw_net_can)

    DEALLOCATE( slurb_tile%rad_lw_net_roof)
    DEALLOCATE( slurb_tile%rad_sw_net_roof)
    DEALLOCATE( slurb_tile%rad_lw_net_road)
    DEALLOCATE( slurb_tile%rad_sw_net_road)
    DEALLOCATE( slurb_tile%rad_sw_in_road)
    DEALLOCATE( slurb_tile%rad_lw_net_wall_a)
    DEALLOCATE( slurb_tile%rad_sw_net_wall_a)
    DEALLOCATE( slurb_tile%rad_lw_net_wall_b)
    DEALLOCATE( slurb_tile%rad_sw_net_wall_b)
    DEALLOCATE( slurb_tile%rad_lw_net_win_a)
    DEALLOCATE( slurb_tile%rad_sw_net_win_a)
    DEALLOCATE( slurb_tile%rad_sw_in_win_a)
    DEALLOCATE( slurb_tile%rad_lw_net_win_b)
    DEALLOCATE( slurb_tile%rad_sw_net_win_b)
    DEALLOCATE( slurb_tile%rad_sw_in_win_b)

    DEALLOCATE( slurb_tile%pt_can)
    DEALLOCATE( slurb_tile%uv_abs_can)
    DEALLOCATE( slurb_tile%uv_eff_can)
    DEALLOCATE( slurb_tile%us_can)
    DEALLOCATE( slurb_tile%rib_can)
    DEALLOCATE( slurb_tile%ol_can)

    DEALLOCATE( slurb_tile%rib_roof)
    DEALLOCATE( slurb_tile%ol_roof)
    DEALLOCATE( slurb_tile%rib_road)
    DEALLOCATE( slurb_tile%ol_road)

    DEALLOCATE( slurb_tile%us_roof)
    DEALLOCATE( slurb_tile%us_road)

    DEALLOCATE( slurb_tile%hw_can)
    DEALLOCATE( slurb_tile%anisotropic_canyon)
    DEALLOCATE( slurb_tile%theta_can)
    DEALLOCATE( slurb_tile%h_bld)
    DEALLOCATE( slurb_tile%f_bld)
    DEALLOCATE( slurb_tile%f_bld_frn)
    DEALLOCATE( slurb_tile%f_win)
    DEALLOCATE( slurb_tile%svf_road)
    DEALLOCATE( slurb_tile%svf_wall)
    DEALLOCATE( slurb_tile%z0_urb)

    DEALLOCATE( slurb_tile%rah_roof)
    DEALLOCATE( slurb_tile%rah_road)
    DEALLOCATE( slurb_tile%rah_can)

    IF ( facade_rah_doe )  THEN
       DEALLOCATE( slurb_tile%rah_wall_a)
       DEALLOCATE( slurb_tile%rah_wall_b)
       DEALLOCATE( slurb_tile%rah_win_a)
       DEALLOCATE( slurb_tile%rah_win_b)
    ELSE
       DEALLOCATE( slurb_tile%rah_facade)
    ENDIF

    DEALLOCATE( slurb_tile%lambda_roof)
    DEALLOCATE( slurb_tile%c_roof)
    DEALLOCATE( slurb_tile%albedo_roof)
    DEALLOCATE( slurb_tile%emiss_roof)
    DEALLOCATE( slurb_tile%z0_roof)
    DEALLOCATE( slurb_tile%z0h_roof)
    DEALLOCATE( slurb_tile%lambda_wall)
    DEALLOCATE( slurb_tile%c_wall)
    DEALLOCATE( slurb_tile%albedo_wall)
    DEALLOCATE( slurb_tile%emiss_wall)
    DEALLOCATE( slurb_tile%z0_wall)
    DEALLOCATE( slurb_tile%lambda_win)
    DEALLOCATE( slurb_tile%c_win)
    DEALLOCATE( slurb_tile%albedo_wall_win)
    DEALLOCATE( slurb_tile%albedo_win)
    DEALLOCATE( slurb_tile%emiss_win)
    DEALLOCATE( slurb_tile%transmissivity_win)
    DEALLOCATE( slurb_tile%absorption_win)
    DEALLOCATE( slurb_tile%lambda_road)
    DEALLOCATE( slurb_tile%c_road)
    DEALLOCATE( slurb_tile%albedo_road)
    DEALLOCATE( slurb_tile%emiss_road)
    DEALLOCATE( slurb_tile%z0_road)
    DEALLOCATE( slurb_tile%z0h_road)

    DEALLOCATE( slurb_tile%conductivity_roof)
    DEALLOCATE( slurb_tile%conductivity_wall)
    DEALLOCATE( slurb_tile%conductivity_win)
    DEALLOCATE( slurb_tile%conductivity_road)

    DEALLOCATE( slurb_tile%z_mo)
    DEALLOCATE( slurb_tile%z_mo_can)
    DEALLOCATE( slurb_tile%uv_abs_can_coef)
    DEALLOCATE( slurb_tile%wall_hor_a_ratio)


    DEALLOCATE( slurb_tile%lw_roof_coef)
    DEALLOCATE( slurb_tile%lw_road_coef)
    DEALLOCATE( slurb_tile%lw_wall_coef)
    DEALLOCATE( slurb_tile%lw_win_coef)
    DEALLOCATE( slurb_tile%sw_ref_denom)

    DEALLOCATE( slurb_tile%us_urb)
    DEALLOCATE( slurb_tile%uv_eff1)
    DEALLOCATE( slurb_tile%uv_abs1)
    DEALLOCATE( slurb_tile%pt1)

    DEALLOCATE( t_can_1)
    DEALLOCATE( t_can_2)
    DEALLOCATE( t_wall_a_1)
    DEALLOCATE( t_wall_a_2)
    DEALLOCATE( t_wall_b_1)
    DEALLOCATE( t_wall_b_2)
    DEALLOCATE( t_win_a_1)
    DEALLOCATE( t_win_a_2)
    DEALLOCATE( t_win_b_1)
    DEALLOCATE( t_win_b_2)
    DEALLOCATE( t_roof_1)
    DEALLOCATE( t_roof_2)
    DEALLOCATE( t_road_1)
    DEALLOCATE( t_road_2)

    IF ( moist_physics )  THEN
       DEALLOCATE( slurb_tile%tq_can)
       DEALLOCATE( slurb_tile%tm_liq_roof)
       DEALLOCATE( slurb_tile%tm_liq_road)

       DEALLOCATE( slurb_tile%vpt_roof)
       DEALLOCATE( slurb_tile%vpt_road)

       DEALLOCATE( slurb_tile%q_roof)
       DEALLOCATE( slurb_tile%q_road)
       DEALLOCATE( slurb_tile%qs_roof)
       DEALLOCATE( slurb_tile%qs_road)

       DEALLOCATE( slurb_tile%qsws_can)
       DEALLOCATE( slurb_tile%qsws_roof)
       DEALLOCATE( slurb_tile%qsws_road)
       DEALLOCATE( slurb_tile%qsws_liq_roof)
       DEALLOCATE( slurb_tile%qsws_liq_road)

       DEALLOCATE( slurb_tile%c_liq_roof)
       DEALLOCATE( slurb_tile%c_liq_road)

       DEALLOCATE( slurb_tile%vpt_can)

       DEALLOCATE( slurb_tile%q1)
       DEALLOCATE( slurb_tile%vpt1)

       DEALLOCATE( slurb_tile%qsws_external)

       DEALLOCATE( q_can_1)
       DEALLOCATE( q_can_2)
       DEALLOCATE( m_liq_roof_1)
       DEALLOCATE( m_liq_roof_2)
       DEALLOCATE( m_liq_road_1)
       DEALLOCATE( m_liq_road_2)
    ENDIF

end subroutine exitslurb

subroutine do_slurb
    use modglobal, only: ntrun

    call slurb_update_external_vars

    call calc_urban_resistances

    call calc_canyon_resistances

    call slurb_energy_balance_model

    CALL slurb_canyon_model

    CALL slurb_urban_aggregation_model

    ! CALL slurb_atmospheric_model_coupler

    call slurb_swap_timelevel(mod(ntrun, 2))
end subroutine do_slurb


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates model external variables, e.g. the variables defined at the first atmospheric level based
!> on atmospheric simulation state as well as the temporally dynamic SLUrb input variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_update_external_vars

    use modglobal, only : cp, rlv, cu, cv, i1, j1, ep
    use modfields, only : ql0, u0, v0, qt0, exnf, thl0, rhobf
    implicit none
    INTEGER ::  i      !< loop index
    INTEGER::  j      !< loop index
    INTEGER ::  k_atm  !< k index of the first atmospheric level
    INTEGER ::  m      !< loop index of surface tiles
    INTEGER ::  t      !< current timestep index
    INTEGER ::  tm     !< previous timestep index

    REAL(field_r) ::  fac_dt  !< factor for linear interpolation between timesteps
    REAL(field_r) ::  vtws    !< buoyancy flux (m K s^-1)
    REAL(field_r) ::  ws      !< free convection velocity scale

    real(field_r) :: du, dv

    real :: rho_cp !< cp * rho (J m^-3 K^-1)
    rho_cp = cp * rho_air_zw(1) !TODOSELF


    do j=2,j1
      do i=2,i1
        k_atm = 1
    !        i = slurb_tile%i(m)
    !        j = slurb_tile%j(m)
    !        k_atm = topo_top_ind(j,i,0) + 1
    ! !
    ! !--    Calculate the pt, vpt and q for atmosphere depending on what modules are enabled.
    !        IF ( bulk_cloud_model )  THEN
    !           slurb_tile%pt1(m) = pt(k_atm,j,i) + lv_d_cp * (1 / exnf(k_atm)) * ql(k_atm,j,i)
    !           slurb_tile%q1(m) = q(k_atm,j,i) - ql(k_atm,j,i)
    !           slurb_tile%vpt1(m) = slurb_tile%pt1(m) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(m) )
    !        ELSEIF ( cloud_droplets )  THEN
    !           slurb_tile%pt1(m) = pt(k_atm,j,i) + lv_d_cp * (1 / exnf(k_atm)) * ql(k_atm,j,i)
    !           slurb_tile%q1(m) = q(k_atm,j,i)
    !           slurb_tile%vpt1(m) = slurb_tile%pt1(m) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(m) )
    !        ELSE
    !           slurb_tile%pt1(m) = pt(k_atm,j,i)
    !           IF ( moist_physics )  THEN
    !              slurb_tile%q1(m) = q(k_atm,j,i)
    !              slurb_tile%vpt1(m) = slurb_tile%pt1(m) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(m) )
    !           ENDIF
    !        ENDIF

    !        slurb_tile%uv_abs1(m) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
    !                                ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 )

            ! K = K + J/kg /(J/kg K^-1) * (kg/kg)
            slurb_tile%pt1(i,j)  = thl0(i, j, k_atm) + (rlv/(cp * exnf(k_atm)))  * ql0(i,j,k_atm)
            slurb_tile%q1(i,j)   = qt0(i, j, k_atm) - ql0(i, j, k_atm) !TODOSELF BUG
            slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )

            du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
            dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
            slurb_tile%uv_abs1(i,j) = sqrt(du**2 + dv**2)
            ! slurb_tile%uv_abs1(i,j) = max(0.1, sqrt(du**2 + dv**2)) DALES VERSION




            !--    Calculate surface-parallel absolute velocity uv_eff1 at cell center using
            !--    free convection scale (w_star, for unstable cases).
            ! m K s^-1 = (J^-1 kg K)(kg^-1 m^3) J s^-1 m^-2 + (J/kg)/(J/kg/K)*W/m^2
            vtws = (1/(rho_cp)) * slurb_tile%shf_urb(i,j) + (rlv / cp) * slurb_tile%qsws_urb(i,j) ! (m K s^-1)
            !
            !--    No scaling for stable cases:
            vtws = MERGE( vtws, 0.0_field_r, vtws > 0.0_field_r )
            ! m/s = (m s^-2 K^-1 m ?)^(1/3)
            ! m/s = (m^2 s^-2 K^-1 (m K s^-1))^(1/3)
            ! m/s = (m^3 s^-3)^(1/3)
            ws = ( g / slurb_tile%pt1(i,j) * slurb_tile%z_mo(i,j) * vtws )**( 1.0_field_r / 3.0_field_r )  ! (m s^-1)

            !    slurb_tile%uv_eff1(m) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
            !                            ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 + ws**2 )
            slurb_tile%uv_eff1(i, j) = sqrt(du**2 + dv**2 + ws**2)
        enddo
    enddo

    ! IF ( slurb_dynamic%ntime > 0 )  CALL update_dynamic_inputs

end subroutine slurb_update_external_vars
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the heat and momentum fluxes between the atmosphere and the urban surface.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_urban_resistances
    use modglobal, only : i1, j1
    implicit none
    integer i, j


    !
    !-- Compute friction velocity and aerodynamic resistance for momentum for the whole urban surface.
    !-- As SLUrb doesn't explicitly compute the momentum flux for each individual surface, and as
    !-- pressure drag needs to be included in the total urban drag, us_urb and rah_urb are computed
    !-- using given roughness length for whole urban fabric (z0_urb, user input). To compute the MOST
    !-- stability corrections, we follow the SURFEX implementation where weighted pt/vpt from canyons
    !-- and roofs is used to represent the pt/vpt at roof level. For urban heat fluxes, aggregated
    !-- values from roofs and canyons are directly used, so rah_urb is not needed.
        IF ( moist_physics )  THEN

        do j=2,j1
            do i=2,i1
                pt_surface(i,j) = slurb_tile%f_bld(i,j)              * slurb_tile%vpt_roof(i,j) +                          &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%vpt_can(i,j)
                CALL calc_rib( slurb_tile%vpt1(i,j), pt_surface(i,j), slurb_tile%rib_urb(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo
       

    ELSE

        do j=2,j1
            do i=2,i1
                pt_surface(i,j) = slurb_tile%f_bld(i,j)              * slurb_tile%pt_roof(i,j) +                           &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%pt_can(i,j)
                CALL calc_rib( slurb_tile%pt1(i,j), pt_surface(i,j), slurb_tile%rib_urb(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo

    ENDIF

    do j=2,j1
        do i=2,i1
            CALL calc_ol( ln_z_z0_urb(i,j), ln_z_z0_urb(i,j), slurb_tile%ol_urb(i,j), slurb_tile%rib_urb(i,j), slurb_tile%z0_urb(i,j),       &
                  slurb_tile%z0_urb(i,j), slurb_tile%z_mo(i,j) )
        enddo
    enddo



    do j=2,j1
        do i=2,i1
       slurb_tile%us_urb(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                  &
                        ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                   &
                          psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_urb(i,j) ) +                                 &
                          psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_urb(i,j) ) )
        enddo
    enddo

    !
    !-- Ensure physical friction velocity (might be needed due to instabilities in e.g. initialization)
    do j=2,j1
        do i=2,i1
       IF ( slurb_tile%us_urb(i,j) <= us_min ) slurb_tile%us_urb(i,j) = us_min

       slurb_tile%ram_urb(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_urb(i,j) ) *                                     &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                  &
                           psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_urb(i,j) ) +                                &
                           psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_urb(i,j) ) )

       IF ( slurb_tile%ram_urb(i,j) < ram_min )  slurb_tile%ram_urb(i,j) = ram_min
        enddo
    enddo
    !
    !-- For street canyons, effective mixing between canyon half-height and roof height is assumed, thus
    !-- z_mo is used as as reference height when considering atmosphere-street canyon air mixing.
    !-- This is equivalent to mixing of canyon air between the roof top level and the first atm grid
    !-- level. For canyons, roughness length for the whole urban fabric (z0_urb) is used instead of
    !-- local z0m/z0h. This is based on an assumption that turbulence can effectively mix the two air
    !-- masses (canyon air and the atmospheric air). The same assumption is used in TEB/SURFEX.
    !-- Using e.g. the Kanda et al. (2007) parametrization or any other surface parametrization for
    !-- canyon z0h would yield unrealistically low mixing.
    !
    !-- Update z0h for roofs following Kanda et al. (2007) parametrization if enabled.
    IF ( roughness_kanda )  THEN
        do j=2,j1
            do i=2,i1
                slurb_tile%z0h_roof(i,j) = slurb_tile%z0_roof(i,j) * 7.4_field_r *                                            &
                                    EXP( -1.29_field_r * SQRT( SQRT( slurb_tile%z0_roof(i,j) * slurb_tile%us_roof(i,j) /       &
                                                                1.461E-5_field_r) ) )
                ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) )
            enddo
        enddo
    ENDIF

    IF ( moist_physics )  THEN
        do j=2,j1
            do i=2, i1
                CALL calc_rib( slurb_tile%vpt1(i,j), slurb_tile%vpt_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
                CALL calc_rib( slurb_tile%vpt1(i,j), slurb_tile%vpt_can(i,j),  slurb_tile%rib_can(i,j),  slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
        enddo
    enddo
    ELSE
        do j=2,j1
            do i=2, i1
                CALL calc_rib( slurb_tile%pt1(i,j), slurb_tile%pt_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
                CALL calc_rib( slurb_tile%pt1(i,j), slurb_tile%pt_can(i,j),  slurb_tile%rib_can(i,j),  slurb_tile%uv_eff1(i,j), slurb_tile%z_mo(i,j) )
            enddo
        enddo
    ENDIF
    do j=2,j1
        do i=2, i1
            ! write(*,*), "now_ol"
            CALL calc_ol( ln_z_z0_roof(i,j), ln_z_z0h_roof(i,j), slurb_tile%ol_roof(i,j), slurb_tile%rib_roof(i,j), slurb_tile%z0_roof(i,j), &
                        slurb_tile%z0h_roof(i,j), slurb_tile%z_mo(i,j) )
            ! write(*,*), i,j
            CALL calc_ol( ln_z_z0_urb(i,j), ln_z_z0_urb(i,j), slurb_tile%ol_can(i,j), slurb_tile%rib_can(i,j), slurb_tile%z0_urb(i,j),       &
                        slurb_tile%z0_urb(i,j), slurb_tile%z_mo(i,j) )
            ! write(*,*), "alldone"
        enddo
    enddo
    ! calc_obuk_dirichlet( &
    !                 tile%obuk(i,j), du_tot(i,j), tile%db(i,j), real(zf(1), 8), tile%z0m(i,j), tile%z0h(i,j))

    !
    !-- Compute the local friction velocity for roof and canyon.
    do j=2,j1
      do i=2,i1
       slurb_tile%us_roof(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                 &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_roof(i,j) ) -                                 &
                           psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_roof(i,j) ) +                               &
                           psi_m( slurb_tile%z0_roof(i,j) / slurb_tile%ol_roof(i,j) ) )

    !
    !--    For canyons, use urban roughness length (assume the air mixes efficiently
    !--    between the canyon air and atmosphere).
       slurb_tile%us_can(i,j) = kappa * slurb_tile%uv_eff1(i,j) /                                                  &
                        ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                   &
                          psi_m( slurb_tile%z_mo(i,j) / slurb_tile%ol_can(i,j) ) +                                 &
                          psi_m( slurb_tile%z0_urb(i,j) / slurb_tile%ol_can(i,j) ) )

    !
    !--    Ensure physical friction velocity.
       IF ( slurb_tile%us_roof(i,j) <= us_min )  slurb_tile%us_roof(i,j) = us_min
      enddo
    enddo

    !
    !-- Compute the aerodynamic resistances for heat.
    do j=2,j1
      do i=2,i1
       slurb_tile%rah_roof(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_roof(i,j) ) *                                   &
                          ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) ) -                               &
                            psi_h( slurb_tile%z_mo(i,j) / slurb_tile%ol_roof(i,j) ) +                              &
                            psi_h( slurb_tile%z0h_roof(i,j) / slurb_tile%ol_roof(i,j) ) )

       slurb_tile%rah_can(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_can(i,j) ) *                                     &
                         ( LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j) ) -                                  &
                           psi_h( slurb_tile%z_mo(i,j) / slurb_tile%ol_can(i,j) ) +                                &
                           psi_h( slurb_tile%z0_urb(i,j) / slurb_tile%ol_can(i,j) ) )

       IF ( slurb_tile%rah_roof(i,j) < rah_min )  slurb_tile%rah_roof(i,j) = rah_min
       IF ( slurb_tile%rah_roof(i,j) > rah_max )  slurb_tile%rah_roof(i,j) = rah_max
    !
    !--    Use ram_min for canyon air as turbulence is able to mix the air.
       IF ( slurb_tile%rah_can(i,j) < ram_min )  slurb_tile%rah_can(i,j) = ram_min
      enddo
    enddo

 END SUBROUTINE calc_urban_resistances



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Model for the surface resistances within the street canyon.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_canyon_resistances
    use modglobal, only : i1, j1
    implicit none

    integer ::  i       !< loop index x-direction
    integer ::  j       !< loop index y-direction
    integer ::  k_topo  !< k-index of topography
    ! integer ::  m       !< running index of surface tiles


    !
    !-- Calculate logarithms of ratio z/z0.
    !>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
    !>        and stored in surf_slurb, like it is done for the other surface types, too.
    do j=2,j1
      do i=2,i1
       ln_z_z0_road(i,j)  = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0_road(i,j)  )
       ln_z_z0h_road(i,j) = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0h_road(i,j) )
       ln_z_z0_roof(i,j)  = LOG( slurb_tile%z_mo(i,j)     / slurb_tile%z0_roof(i,j)  )
       ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j)     / slurb_tile%z0h_roof(i,j) )
      enddo
    enddo

    !
    !-- Update z0h for roads following Kanda et al. (2007) parametrization if necessary.
    IF ( roughness_kanda )  THEN
    do j=2,j1
      do i=2,i1
            slurb_tile%z0h_road(i,j) = slurb_tile%z0_road(i,j) * 7.4_field_r *                                            &
                                EXP( -1.29_field_r * SQRT( SQRT( slurb_tile%z0_road(i,j) * slurb_tile%us_road(i,j) /       &
                                                            1.461E-5_field_r) ) )
            ln_z_z0h_road(i,j) = LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0h_road(i,j) )
             enddo
        enddo
    ENDIF

    !
    !-- Compute the new Obukhov length for road.
    IF ( moist_physics )  THEN
        do j=2,j1
            do i=2,i1
                CALL calc_rib( slurb_tile%vpt_can(i,j), slurb_tile%vpt_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%uv_eff_can(i,j),        &
                             slurb_tile%z_mo_can(i,j) )
            enddo
        enddo
    ELSE
        do j=2,j1
            do i=2,i1
                CALL calc_rib( slurb_tile%pt_can(i,j), slurb_tile%pt_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%uv_eff_can(i,j),          &
                               slurb_tile%z_mo_can(i,j) )
            enddo
        enddo
    ENDIF

    do j=2,j1
      do i=2,i1
        CALL calc_ol( ln_z_z0_road(i,j), ln_z_z0h_road(i,j), slurb_tile%ol_road(i,j), slurb_tile%rib_road(i,j), slurb_tile%z0_road(i,j), &
                    slurb_tile%z0h_road(i,j), slurb_tile%z_mo_can(i,j) )
      enddo
    enddo

    !
    !-- Compute the local friction velocity for roads.
    do j=2,j1
      do i=2,i1
        slurb_tile%us_road(i,j) = kappa * slurb_tile%uv_eff_can(i,j) /                                              &
                            ( LOG( slurb_tile%z_mo_can(i,j) / slurb_tile%z0_road(i,j) ) -                             &
                            psi_m( slurb_tile%z_mo_can(i,j) / slurb_tile%ol_road(i,j) ) +                           &
                            psi_m( slurb_tile%z0_road(i,j) / slurb_tile%ol_road(i,j) ) )
        enddo
    enddo

    !
    !-- The resistance between the street canyon air and facades (walls and windows).
    IF ( facade_rah_doe )  THEN

        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0)  ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_wall_a(i,j) = rah_doe2( k_topo, slurb_tile%t_can(i,j), slurb_tile%t_wall_a(nzt_wall,i,j),         &
                                                slurb_tile%uv_eff_can(i,j), .TRUE. )
                IF ( slurb_tile%rah_wall_a(i,j) < rah_min )  slurb_tile%rah_wall_a(i,j) = rah_min
                IF ( slurb_tile%rah_wall_a(i,j) > rah_max )  slurb_tile%rah_wall_a(i,j) = rah_max
                IF ( slurb_tile%f_win(i,j) /= 0.0_field_r )  THEN
                    slurb_tile%rah_win_a(i,j) = rah_doe2( k_topo, slurb_tile%t_can(i,j), slurb_tile%t_win_a(nzt_win,i,j),         &
                                                slurb_tile%uv_eff_can(i,j), .FALSE. )
                    IF ( slurb_tile%rah_win_a(i,j) < rah_min )  slurb_tile%rah_win_a(i,j) = rah_min
                    IF ( slurb_tile%rah_win_a(i,j) > rah_max )  slurb_tile%rah_win_a(i,j) = rah_max
                ENDIF

                IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
                    slurb_tile%rah_wall_b(i,j) = rah_doe2( k_topo, slurb_tile%t_can(i,j), slurb_tile%t_wall_b(nzt_wall,i,j),      &
                                                    slurb_tile%uv_eff_can(i,j), .TRUE. )
                    IF ( slurb_tile%rah_wall_b(i,j) < rah_min )  slurb_tile%rah_wall_b(i,j) = rah_min
                    IF ( slurb_tile%rah_wall_b(i,j) > rah_max )  slurb_tile%rah_wall_b(i,j) = rah_max
                    IF ( slurb_tile%f_win(i,j) /= 0.0_field_r )  THEN
                        slurb_tile%rah_win_b(i,j) = rah_doe2( k_topo, slurb_tile%t_can(i,j), slurb_tile%t_win_b(nzt_win,i,j),      &
                                                    slurb_tile%uv_eff_can(i,j), .FALSE. )
                        IF ( slurb_tile%rah_win_b(i,j) < rah_min )  slurb_tile%rah_win_b(i,j) = rah_min
                        IF ( slurb_tile%rah_win_b(i,j) > rah_max )  slurb_tile%rah_win_b(i,j) = rah_max
                    ENDIF
                ENDIF
             enddo
        enddo

    ELSEIF ( facade_rah_kray )  THEN

        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0) ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_facade(i,j) = rah_kray( k_topo, slurb_tile%z0_wall(i,j), slurb_tile%uv_eff_can(i,j) )
                IF ( slurb_tile%rah_facade(i,j) < rah_min )  slurb_tile%rah_facade(i,j) = rah_min
                IF ( slurb_tile%rah_facade(i,j) > rah_max )  slurb_tile%rah_facade(i,j) = rah_max
            enddo
        enddo

    ELSEIF ( facade_rah_rowley )  THEN
    !
    !--    Rowley et al. (1930) , Cole and Sturrock (1977)  Mills (1993).
        do j=2,j1
            do i=2,i1
                ! k_topo = topo_top_ind(j,i,0) ! ZELFTODO
                k_topo = 1
                slurb_tile%rah_facade(i,j) = cp * rho_air_zw(k_topo) / ( 11.8_field_r + 4.2_field_r * slurb_tile%uv_eff_can(i,j) )
                IF ( slurb_tile%rah_facade(i,j) < rah_min )  slurb_tile%rah_facade(i,j) = rah_min
                IF ( slurb_tile%rah_facade(i,j) > rah_max )  slurb_tile%rah_facade(i,j) = rah_max
            enddo
        enddo
    ENDIF

    do j=2,j1
        do i=2,i1
            slurb_tile%rah_road(i,j) = 1.0_field_r / ( kappa * slurb_tile%us_can(i,j) ) *                                    &
                                ( ln_z_z0h_road(i,j) -                                                     &
                                    psi_h( slurb_tile%z_mo_can(i,j) / slurb_tile%ol_road(i,j) ) +                          &
                                    psi_h( slurb_tile%z0h_road(i,j) / slurb_tile%ol_road(i,j) ) )

            IF ( slurb_tile%rah_road(i,j) < rah_min )  slurb_tile%rah_road(i,j) = rah_min
            IF ( slurb_tile%rah_road(i,j) > rah_max )  slurb_tile%rah_road(i,j) = rah_max
        enddo
    enddo

 END SUBROUTINE calc_canyon_resistances



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_ol(ln_z_z0, ln_z_z0h, ol, rib, z0, z0h, z_mo )

    IMPLICIT NONE

    real, intent(in)    ::  ln_z_z0   !< logarithm (z/z0)
    real, intent(in)    ::  ln_z_z0h  !< logarithm (z/z0h)
    real, intent(inout) ::  ol        !< Obukhov length
    real, intent(in)    ::  rib       !< Richardson flux number
    real, intent(in)    ::  z0        !< rougness length for momentum
    real, intent(in)    ::  z0h       !< rougness length for scalar quantities
    real, intent(in)    ::  z_mo      !< constant flux layer height

    integer ::  iter  !< Newton iteration step

    ! LOGICAL ::  convergence_reached  !< convergence switch for vectorization

    real ::  f        !< function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    real ::  f_d_ol   !< derivative of f
    real ::  ol_l     !< lower bound of L for Newton iteration
    real ::  ol_m     !< previous value of L for Newton iteration
    real ::  ol_prev  !< previous time step value of L
    real ::  ol_u     !< upper bound of L for Newton iteration

    ! real ::  ol_prev_vec  !< temporary array required for vectorization

    !
    !-- Calculate the Obukhov length using Newton iteration.

    !
    !--       Store current value in case the Newton iteration fails.
            ol_prev = ol
    !
    !--       Flip the sign of the initial Obukhov length if the stability has changed from stable to
    !--       unstable or vice versa and set it to a moderate value. A moderate value is also chosen,
    !--       if the Obukhov length from the last time step reached the maximum threshold value.
            IF ( rib * ol < 0.0_field_r  .OR.  ABS( ol ) == ol_max )  THEN
                IF ( rib > 0.0_field_r )  ol =  100.0_field_r
                IF ( rib < 0.0_field_r )  ol = -100.0_field_r
            ENDIF
    !
    !--       Iteration to find Obukhov length.
            iter = 0
            DO
                iter = iter + 1
    !
    !--          In case of divergence, use the value of the previous time step.
                IF ( iter > 1000 )  THEN
                ol = ol_prev
                EXIT
                ENDIF

    !
    !--          Calculate step size for central difference.
                ol_m = ol
                ol_l = ol_m - 0.001_field_r * ol_m
                ol_u = ol_m + 0.001_field_r * ol_m

                IF ( ibc_pt_b /= 1 )  THEN
    !
    !--             Calculate f = Ri - [...]/[...]^2 = 0.
                f = rib - ( z_mo / ol_m ) * ( ln_z_z0h - psi_h( z_mo / ol_m )          &
                                                                + psi_h( z0h  / ol_m ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_m )          &
                                                                + psi_m( z0   / ol_m ) )**2
    !
    !--             Calculate df/dL.
                f_d_ol = ( - ( z_mo / ol_u ) * ( ln_z_z0h - psi_h( z_mo / ol_u )          &
                                                                + psi_h( z0h  / ol_u ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_u )          &
                                                                + psi_m( z0   / ol_u ) )**2     &
                            + ( z_mo / ol_l ) * ( ln_z_z0h - psi_h( z_mo / ol_l )          &
                                                                + psi_h( z0h  / ol_l ) )        &
                                                / ( ln_z_z0  - psi_m( z_mo / ol_l )          &
                                                                + psi_m( z0   / ol_l ) )**2     &
                            ) / ( ol_u - ol_l )
                ELSE
    !
    !--             Calculate f = Ri - 1 /[...]^3 = 0.
                f = rib - ( z_mo / ol_m ) /                                                  &
                                ( ln_z_z0 - psi_m( z_mo / ol_m ) + psi_m( z0 / ol_m ) )**3

    !
    !--             Calculate df/dL.
                f_d_ol = ( - ( z_mo / ol_u ) / ( ln_z_z0 - psi_m( z_mo / ol_u )           &
                                                                + psi_m( z0   / ol_u ) )**3      &
                            + ( z_mo / ol_l ) / ( ln_z_z0 - psi_m( z_mo / ol_l )           &
                                                                + psi_m( z0   / ol_l ) )**3      &
                            ) / ( ol_u - ol_l )
                ENDIF
    !
    !--          Calculate new L.
                ol = ol_m - f / f_d_ol
    !
    !--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
    !--          ensure convergence. If the sign is not the same, the above calculated Obukhov length
    !--          obviously overshooted to the opposite side, so the next iteration should start with
    !--          a smaller value.
                IF ( ol * ol_m < 0.0_field_r )  ol = ol_m * 0.5_field_r
    !
    !--          In the deep neutral zone, set L to the maximum allowed value.
                IF ( ABS( ol ) > ol_max )  THEN
                ol = SIGN( ol_max, ol )
                EXIT
                ENDIF
    !
    !--          Assure that Obukhov length does not become zero.
                IF ( ABS( ol ) < ol_min )  THEN
                ol = SIGN( ol_min, ol )
                EXIT
                ENDIF
    !
    !--          Check for convergence.
                IF ( ABS( ( ol - ol_m ) /  ol ) < ol_tol )  EXIT

            ENDDO

 END SUBROUTINE calc_ol



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the bulk Richardson number for given surface (z0) temperature.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rib( pt1, pt_surface, rib, uvw_abs, z_mo )
    use modmpi, only: myid
    implicit none

    real, intent(in)  ::  pt1          !< potential temperature at first grid level
    real, intent(in)  ::  pt_surface   !< skin-surface potential temperature
    real, intent(out) ::  rib          !< Richardson flux number
    real, intent(in)  ::  uvw_abs      !< absolute surface-parallel velocity on grid center
    real, intent(in)  ::  z_mo         !< constant flux layer height

    ! write(*,*) g, z_mo, pt1, pt_surface, uvw_abs, pt1
    !-- Evaluate bulk Richardson number.
    rib = g * z_mo * ( pt1 - pt_surface ) / ( uvw_abs**2 * pt1 + 1.0E-20_field_r )
    ! rib = 5
    !TODOSELF TODOSELF !!!! TURN OFF
    !
    !-- For the SLUrb model, limit to |rib| < |rib_max| to dampen possible instabilities during
    !-- initialization.
    IF ( ABS( rib ) > rib_max )  rib = SIGN( rib_max, rib )

 END SUBROUTINE calc_rib



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for momentum.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_m( zeta )

    IMPLICIT NONE

    REAL(field_r), INTENT(IN) ::  zeta   !< Stability parameter z/L

    REAL(field_r) ::  psi_m  !< Integrated similarity function result
    REAL(field_r) ::  x      !< dummy variable

    REAL(field_r), PARAMETER ::  a = 1.0_field_r            !< constant
    REAL(field_r), PARAMETER ::  b = 0.66666666666_field_r  !< constant
    REAL(field_r), PARAMETER ::  c = 5.0_field_r            !< constant
    REAL(field_r), PARAMETER ::  d = 0.35_field_r           !< constant
    REAL(field_r), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(field_r), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_field_r )  THEN
       x = SQRT( SQRT( 1.0_field_r  - 16.0_field_r * zeta ) )
       psi_m = pi * 0.5_field_r - 2.0_field_r * ATAN( x ) + LOG( ( 1.0_field_r + x )**2                           &
               * ( 1.0_field_r + x**2 ) * 0.125_field_r )
    ELSE

       psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta - bc_d_d
!
!--    Old version for stable conditions (only valid for z/L < 0.5) psi_m = - 5.0_field_r * zeta

    ENDIF

 END FUNCTION psi_m


!--------------------------------------------------------------------------------------------------!
! Description:
!------------
!> Integrated stability function for heat and moisture.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION psi_h( zeta )

    IMPLICIT NONE

    REAL(field_r), INTENT(IN) ::  zeta   !< stability parameter z/L

    REAL(field_r) ::  psi_h  !< integrated similarity function result
    REAL(field_r) ::  x      !< dummy variable

    REAL(field_r), PARAMETER ::  a = 1.0_field_r            !< constant
    REAL(field_r), PARAMETER ::  b = 0.66666666666_field_r  !< constant
    REAL(field_r), PARAMETER ::  c = 5.0_field_r            !< constant
    REAL(field_r), PARAMETER ::  d = 0.35_field_r           !< constant
    REAL(field_r), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(field_r), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_field_r )  THEN
       x = SQRT( 1.0_field_r  - 16.0_field_r * zeta )
       psi_h = 2.0_field_r * LOG( (1.0_field_r + x ) / 2.0_field_r )
    ELSE
       psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_field_r                                 &
               + 0.66666666666_field_r * a * zeta )**1.5_field_r - bc_d_d + 1.0_field_r
!
!--    Old version for stable conditions (only valid for z/L < 0.5)
!--    psi_h = - 5.0_field_r * zeta
    ENDIF

 END FUNCTION psi_h


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates stability function for momentum
!>
!> @author Hauke Wurps
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION phi_m( zeta )

    IMPLICIT NONE

    REAL(field_r), INTENT(IN) ::  zeta   !< stability parameter z/L

    REAL(field_r) ::  phi_m  !< value of the function

    REAL(field_r), PARAMETER ::  a = 16.0_field_r  !< constant
    REAL(field_r), PARAMETER ::  c = 5.0_field_r   !< constant

    IF ( zeta < 0.0_field_r )  THEN
       phi_m = 1.0_field_r / SQRT( SQRT( 1.0_field_r - a * zeta ) )
    ELSE
       phi_m = 1.0_field_r + c * zeta
    ENDIF

 END FUNCTION phi_m



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute aerodynamic resistance for heat for vertical surfaces following DOE-2 parametrization,
!> which takes natural convection into account. Average of leeward and windward sides.
!> Source: EnegyPlus 23.2.0 Engineering Reference p.68.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION rah_doe2( k_topo, t_air, t_surf, u_eff, rough )

    LOGICAL, INTENT(IN) ::  rough  !< flag for rough surface, true for walls, false for windows

    integer, INTENT(IN) ::  k_topo  !< k-index of topography

    REAL(field_r), INTENT(IN) ::  t_air   !< temperature of adjacent air
    REAL(field_r), INTENT(IN) ::  t_surf  !< surface temperature
    REAL(field_r), INTENT(IN) ::  u_eff   !< effective wind speed

    REAL(field_r), PARAMETER ::  r_f = 1.52_field_r  !< surface roughness multiplier

    REAL(field_r) ::  chtcn       !< convective heat transfer coefficient for natural convection
    REAL(field_r) ::  chtcs       !< convective heat transfer coefficient for smooth surface
    REAL(field_r) ::  chtcs_lee   !< convective heat transfer coefficient for smooth surface (leeward)
    REAL(field_r) ::  chtcs_wind  !< convective heat transfer coefficient for smooth surface (windward)
    REAL(field_r) ::  rah_doe2    !< resulting resistance


    chtcn = 1.31_field_r * ABS( t_air - t_surf )**0.33333_field_r

    chtcs_lee  = SQRT( chtcn**2 + ( 2.86_field_r * u_eff**0.617_field_r )**2 )
    chtcs_wind = SQRT( chtcn**2 + ( 2.38_field_r * u_eff**0.89_field_r  )**2 )

    chtcs = 0.5 * ( chtcs_lee + chtcs_wind )

    IF ( rough )  THEN
       rah_doe2 = cp * rho_air_zw(k_topo) / ( chtcn + r_f * ( chtcs - chtcn ) )
    ELSE
       rah_doe2 = cp * rho_air_zw(k_topo) / chtcs
    ENDIF

 END FUNCTION rah_doe2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute aerodynamic resistance for heat for vertical surfaces following
!> Krayenhoff & Voogt (2007).
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION rah_kray( k_topo, z0, u_eff )

    integer, INTENT(IN) ::  k_topo  !< k-index of topography

    REAL(field_r), INTENT(IN) ::  u_eff  !< effective wind speed
    REAL(field_r), INTENT(IN) ::  z0     !< roughness length for momentum

    REAL(field_r) ::  kray_coeff  !< denominator for the parametrization
    REAL(field_r) ::  rah_kray    !< resulting resistance


!
!-- Compute denominator first, ensuring it is a positive number.
    kray_coeff = MAX( z0 * 1000.0_field_r * ( 11.8_field_r + 4.2_field_r * u_eff ) - 4.0_field_r, 1.0E-3_field_r )

    rah_kray = cp * rho_air_zw(k_topo) / kray_coeff

 END FUNCTION



!--------------------------------------------------------------------------------------------------!
!   MODLULE PREDEFINED PARAMETERS
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Default parameters for the building types. These are based on the urban surface mod.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_default_pars

!
!-- Residential, < 1950.
    building_pars_slurb(:,1) = (/                                                                  &
       0.18_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.90_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_field_r,        &   !< parameter 43  - [-] transmissivity
       0.12_field_r,        &   !< parameter 44  - [-] albedo
       0.91_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, 1950 - 2000.
    building_pars_slurb(:,2) = (/                                                                  &
       0.25_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_field_r,        &   !< parameter 14  - [-] albedo
       0.95_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_field_r,        &   !< parameter 43  - [-] transmissivity
       0.15_field_r,        &   !< parameter 44  - [-] albedo
       0.87_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Residential, > 2000.
    building_pars_slurb(:,3) = (/                                                                  &
       0.29_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.92_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_field_r,        &   !< parameter 43  - [-] transmissivity
       0.18_field_r,        &   !< parameter 44  - [-] albedo
       0.80_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, < 1950.
    building_pars_slurb(:,4) = (/                                                                  &
       0.18_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.02_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.51200E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.70965E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.120_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.90_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.18_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.18_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       1.5120E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.5120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.810_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.810_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.45_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.45_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.45_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.45_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.70_field_r,        &   !< parameter 43  - [-] transmissivity
       0.12_field_r,        &   !< parameter 44  - [-] albedo
       0.91_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, 1950 - 2000.
    building_pars_slurb(:,5) = (/                                                                  &
       0.25_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.15_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.20_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       1.70000E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.07920E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       2.11200E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.160_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.046_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       2.100_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.10_field_r,        &   !< parameter 14  - [-] albedo
       0.95_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.06_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.24_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       2.1120E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.046_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       2.100_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.30_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.18_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.18_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.18_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.18_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.65_field_r,        &   !< parameter 43  - [-] transmissivity
       0.15_field_r,        &   !< parameter 44  - [-] albedo
       0.87_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Office, > 2000.
    building_pars_slurb(:,6) = (/                                                                  &
       0.29_field_r,        &   !< parameter 0   - [-] window fraction
       0.02_field_r,        &   !< parameter 1   - [m] 1st roof layer thickness (outside)
       0.04_field_r,        &   !< parameter 2   - [m] 2nd roof layer thickness
       0.30_field_r,        &   !< parameter 3   - [m] 3rd roof layer thickness
       0.02_field_r,        &   !< parameter 4   - [m] 4th roof layer thickness (inside)
       3.75360E6_field_r,   &   !< parameter 5   - [J/(m3*K)] specific heat capacity 1st roof layer (outside)
       0.70965E6_field_r,   &   !< parameter 6   - [J/(m3*K)] specific heat capacity 2nd roof layer
       0.07920E6_field_r,   &   !< parameter 7   - [J/(m3*K)] specific heat capacity 3rd roof layer
       1.52600E6_field_r,   &   !< parameter 8   - [J/(m3*K)] specific heat capacity 4th roof layer (inside)
       0.520_field_r,       &   !< parameter 9   - [W/(m*K)] thermal conductivity 1st roof layer (outside)
       0.120_field_r,       &   !< parameter 10  - [W/(m*K)] thermal conductivity 2nd roof layer
       0.035_field_r,       &   !< parameter 11  - [W/(m*K)] thermal conductivity 3rd roof layer
       0.700_field_r,       &   !< parameter 12  - [W/(m*K)] thermal conductivity 4th roof layer (inside)
       0.15_field_r,        &   !< parameter 13  - [m] z0 roughness length for momentum
       0.17_field_r,        &   !< parameter 14  - [-] albedo
       0.92_field_r,        &   !< parameter 15  - [-] emissivity
       0.02_field_r,        &   !< parameter 16  - [m] 1st wall layer thickness (outside)
       0.20_field_r,        &   !< parameter 17  - [m] 2nd wall layer thickness
       0.36_field_r,        &   !< parameter 18  - [m] 3rd wall layer thickness
       0.02_field_r,        &   !< parameter 19  - [m] 4th wall layer thickness
       1.5200E6_field_r,    &   !< parameter 20  - [J/(m3*K)] specific heat capacity 1st wall layer (outside)
       0.0792E6_field_r,    &   !< parameter 21  - [J/(m3*K)] specific heat capacity 2nd wall layer
       1.3400E6_field_r,    &   !< parameter 22  - [J/(m3*K)] specific heat capacity 3rd wall layer
       1.5260E6_field_r,    &   !< parameter 23  - [J/(m3*K)] specific heat capacity 4th wall layer (inside)
       0.930_field_r,       &   !< parameter 24  - [W/(m*K)] thermal conductivity 1st wall layer (outside)
       0.035_field_r,       &   !< parameter 25  - [W/(m*K)] thermal conductivity 2nd wall layer
       0.680_field_r,       &   !< parameter 26  - [W/(m*K)] thermal conductivity 3rd wall layer
       0.700_field_r,       &   !< parameter 27  - [W/(m*K)] thermal conductivity 4th wall layer (inside)
       0.001_field_r,       &   !< parameter 28  - [m] z0 roughness length for momentum
       0.37_field_r,        &   !< parameter 29  - [-] albedo
       0.93_field_r,        &   !< parameter 30  - [-] emissivity
       0.02_field_r,        &   !< parameter 31  - [m] 1st window layer thickness (glass sheet + air total) (outside)
       0.02_field_r,        &   !< parameter 32  - [m] 2rd window layer thickness
       0.02_field_r,        &   !< parameter 33  - [m] 3rd window layer thickness
       0.02_field_r,        &   !< parameter 34  - [m] 4th window layer thickness (inside)
       1.736E6_field_r,     &   !< parameter 35  - [J/(m3*K)] specific heat capacity 1st window layer (outside)
       1.736E6_field_r,     &   !< parameter 36  - [J/(m3*K)] specific heat capacity 2nd window layer
       1.736E6_field_r,     &   !< parameter 37  - [J/(m3*K)] specific heat capacity 3rd window layer
       1.736E6_field_r,     &   !< parameter 38  - [J/(m3*K)] specific heat capacity 4th window layer (inside)
       0.11_field_r,        &   !< parameter 39  - [W/(m*K)] thermal conductivity 1st window layer (outside)
       0.11_field_r,        &   !< parameter 40  - [W/(m*K)] thermal conductivity 2nd window layer
       0.11_field_r,        &   !< parameter 41  - [W/(m*K)] thermal conductivity 3rd window layer
       0.11_field_r,        &   !< parameter 42  - [W/(m*K)] thermal conductivity 4th window layer (inside)
       0.57_field_r,        &   !< parameter 43  - [-] transmissivity
       0.18_field_r,        &   !< parameter 44  - [-] albedo
       0.80_field_r         &   !< parameter 45  - [-] emissivity
    /)

!
!-- Asphalt concrete mix (I-II), stone aggregate(III), gravel and soil(IV), PALM-LSM default.
    pavement_pars_slurb(:,1) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.00E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.00E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.00_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.00_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.93_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Asphalt concrete (I-II), stone aggregate (III), gravel and soil (IV), Masson et al. (2002).
    pavement_pars_slurb(:,2) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       1.74E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       1.74E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       0.82_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       0.82_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.10_field_r,      &   !< parameter 13  - [-] albedo
       0.95_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Concrete (Portland concrete, I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002) and Yaghoobian et al. (2009).
    pavement_pars_slurb(:,3) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.11E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.11E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       1.51_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       1.51_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.30_field_r,      &   !< parameter 13  - [-] albedo
       0.90_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Sett (I-II), stone aggregate (III), gravel and soil (IV),Masson et al. (2002), Oke (1987)
!-- and Mandanici et al. (2016).
    pavement_pars_slurb(:,4) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.95_field_r       &   !< parameter 14  - [-] emissivity
    /)

!
!-- Pavement stones (I-II), stone aggregate (III), gravel and soil (IV),
!-- Masson et al. (2002), Oke (1987) and Göttsche & Hulley (2012).
    pavement_pars_slurb(:,5) = (/                                                                  &
       0.01_field_r,      &   !< parameter 0   - [m] 1st pavement layer thickness (top)
       0.04_field_r,      &   !< parameter 1   - [m] 2nd pavement layer thickness
       0.20_field_r,      &   !< parameter 2   - [m] 3rd pavement layer thickness
       1.00_field_r,      &   !< parameter 3   - [m] 4th pavement layer thickness (bottom)
       2.25E6_field_r,    &   !< parameter 4   - [J/(m3*K)] heat capacity 1st pavement layer (top)
       2.25E6_field_r,    &   !< parameter 5   - [J/(m3*K)] heat capacity 2nd pavement layer
       2.00E6_field_r,    &   !< parameter 6   - [J/(m3*K)] heat capacity 3rd pavement layer
       1.40E6_field_r,    &   !< parameter 7   - [J/(m3*K)] heat capacity 4th pavement layer (bottom)
       2.19_field_r,      &   !< parameter 8   - [W/(m*K)] thermal conductivity 1st pavement layer (top)
       2.19_field_r,      &   !< parameter 9   - [W/(m*K)] thermal conductivity 2nd pavement layer
       2.10_field_r,      &   !< parameter 10  - [W/(m*K)] thermal conductivity 3rd pavement layer
       0.40_field_r,      &   !< parameter 11  - [W/(m*K)] thermal conductivity 4th pavement layer (bottom)
       5.0E-2_field_r,    &   !< parameter 12  - [m] z0 roughness length for momentum
       0.17_field_r,      &   !< parameter 13  - [-] albedo
       0.93_field_r       &   !< parameter 14  - [-] emissivity
    /)

 END SUBROUTINE slurb_default_pars
 !--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Process parameters dependent on the building/pavement type and properties.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE process_surface_parameters
    use modglobal, only: i1, j1, i2, j2
    implicit none
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  type_tmp  !< array to contain building type temporarily
    integer i,j,k


    ALLOCATE( type_tmp(i2, j2) )
    type_tmp(:,:) = 2 ! TODOSELF buiding type staat voor nu gewoon even vast op 2

    CALL slurb_default_pars

    do j=2,j1
      do i=2,i1
        slurb_tile%f_win(i,j) = building_pars_slurb(0,type_tmp(i,j))

        IF ( n_layers_roofs == 4 )  THEN
            slurb_tile%dz_roof(1,i,j) = building_pars_slurb(1,type_tmp(i,j))
            slurb_tile%dz_roof(2,i,j) = building_pars_slurb(2,type_tmp(i,j))
            slurb_tile%dz_roof(3,i,j) = building_pars_slurb(3,type_tmp(i,j))
            slurb_tile%dz_roof(4,i,j) = building_pars_slurb(4,type_tmp(i,j))

            slurb_tile%c_roof(1,i,j) = building_pars_slurb(5,type_tmp(i,j))
            slurb_tile%c_roof(2,i,j) = building_pars_slurb(6,type_tmp(i,j))
            slurb_tile%c_roof(3,i,j) = building_pars_slurb(7,type_tmp(i,j))
            slurb_tile%c_roof(4,i,j) = building_pars_slurb(8,type_tmp(i,j))

            slurb_tile%lambda_roof(1,i,j) = building_pars_slurb(9,type_tmp(i,j))
            slurb_tile%lambda_roof(2,i,j) = building_pars_slurb(10,type_tmp(i,j))
            slurb_tile%lambda_roof(3,i,j) = building_pars_slurb(11,type_tmp(i,j))
            slurb_tile%lambda_roof(4,i,j) = building_pars_slurb(12,type_tmp(i,j))
        ENDIF

        slurb_tile%z0_roof(i,j)     = building_pars_slurb(13,type_tmp(i,j))
        slurb_tile%z0h_roof(i,j)    = building_pars_slurb(13,type_tmp(i,j)) * 1.0E-2
        slurb_tile%albedo_roof(i,j) = building_pars_slurb(14,type_tmp(i,j))
        slurb_tile%emiss_roof(i,j)  = building_pars_slurb(15,type_tmp(i,j))

        IF ( n_layers_walls == 4 )  THEN
            slurb_tile%dz_wall(1,i,j) = building_pars_slurb(16,type_tmp(i,j))
            slurb_tile%dz_wall(2,i,j) = building_pars_slurb(17,type_tmp(i,j))
            slurb_tile%dz_wall(3,i,j) = building_pars_slurb(18,type_tmp(i,j))
            slurb_tile%dz_wall(4,i,j) = building_pars_slurb(19,type_tmp(i,j))

            slurb_tile%c_wall(1,i,j) = building_pars_slurb(20,type_tmp(i,j))
            slurb_tile%c_wall(2,i,j) = building_pars_slurb(21,type_tmp(i,j))
            slurb_tile%c_wall(3,i,j) = building_pars_slurb(22,type_tmp(i,j))
            slurb_tile%c_wall(4,i,j) = building_pars_slurb(23,type_tmp(i,j))

            slurb_tile%lambda_wall(1,i,j) = building_pars_slurb(24,type_tmp(i,j))
            slurb_tile%lambda_wall(2,i,j) = building_pars_slurb(25,type_tmp(i,j))
            slurb_tile%lambda_wall(3,i,j) = building_pars_slurb(26,type_tmp(i,j))
            slurb_tile%lambda_wall(4,i,j) = building_pars_slurb(27,type_tmp(i,j))
        ENDIF

        slurb_tile%z0_wall(i,j)     = building_pars_slurb(28,type_tmp(i,j))
        slurb_tile%albedo_wall(i,j) = building_pars_slurb(29,type_tmp(i,j))
        slurb_tile%emiss_wall(i,j)  = building_pars_slurb(30,type_tmp(i,j))

        IF ( n_layers_windows == 4 )  THEN
            slurb_tile%dz_win(1,i,j) = building_pars_slurb(31,type_tmp(i,j))
            slurb_tile%dz_win(2,i,j) = building_pars_slurb(32,type_tmp(i,j))
            slurb_tile%dz_win(3,i,j) = building_pars_slurb(33,type_tmp(i,j))
            slurb_tile%dz_win(4,i,j) = building_pars_slurb(34,type_tmp(i,j))

            slurb_tile%c_win(1,i,j) = building_pars_slurb(35,type_tmp(i,j))
            slurb_tile%c_win(2,i,j) = building_pars_slurb(36,type_tmp(i,j))
            slurb_tile%c_win(3,i,j) = building_pars_slurb(37,type_tmp(i,j))
            slurb_tile%c_win(4,i,j) = building_pars_slurb(38,type_tmp(i,j))

            slurb_tile%lambda_win(1,i,j) = building_pars_slurb(39,type_tmp(i,j))
            slurb_tile%lambda_win(2,i,j) = building_pars_slurb(40,type_tmp(i,j))
            slurb_tile%lambda_win(3,i,j) = building_pars_slurb(41,type_tmp(i,j))
            slurb_tile%lambda_win(4,i,j) = building_pars_slurb(42,type_tmp(i,j))
        ENDIF

        slurb_tile%transmissivity_win(i,j) = building_pars_slurb(43,type_tmp(i,j))

        slurb_tile%albedo_win(i,j) = building_pars_slurb(44,type_tmp(i,j))
        slurb_tile%emiss_win(i,j)  = building_pars_slurb(45,type_tmp(i,j))

        enddo
    enddo

!
!-- Process pavement type.
    type_tmp(:,:) = 2 !TODOSELF pavement type staat nu ook gewoon constant op 2

    do j=2,j1
      do i=2,i1
        IF ( n_layers_roads == 4 )  THEN
            slurb_tile%dz_road(1,i,j) = pavement_pars_slurb(0,type_tmp(i,j))
            slurb_tile%dz_road(2,i,j) = pavement_pars_slurb(1,type_tmp(i,j))
            slurb_tile%dz_road(3,i,j) = pavement_pars_slurb(2,type_tmp(i,j))
            slurb_tile%dz_road(4,i,j) = pavement_pars_slurb(3,type_tmp(i,j))

            slurb_tile%c_road(1,i,j) = pavement_pars_slurb(4,type_tmp(i,j))
            slurb_tile%c_road(2,i,j) = pavement_pars_slurb(5,type_tmp(i,j))
            slurb_tile%c_road(3,i,j) = pavement_pars_slurb(6,type_tmp(i,j))
            slurb_tile%c_road(4,i,j) = pavement_pars_slurb(7,type_tmp(i,j))

            slurb_tile%lambda_road(1,i,j) = pavement_pars_slurb(8,type_tmp(i,j))
            slurb_tile%lambda_road(2,i,j) = pavement_pars_slurb(9,type_tmp(i,j))
            slurb_tile%lambda_road(3,i,j) = pavement_pars_slurb(10,type_tmp(i,j))
            slurb_tile%lambda_road(4,i,j) = pavement_pars_slurb(11,type_tmp(i,j))

            slurb_tile%z0_road(i,j)     = pavement_pars_slurb(12,type_tmp(i,j))
            slurb_tile%z0h_road(i,j)    = pavement_pars_slurb(12,type_tmp(i,j)) * 1.0E-2
            slurb_tile%albedo_road(i,j) = pavement_pars_slurb(13,type_tmp(i,j))
            slurb_tile%emiss_road(i,j)  = pavement_pars_slurb(14,type_tmp(i,j))
        ENDIF
        enddo
    enddo
    DEALLOCATE( type_tmp )

!
!-- Process material layer information such as thickness, heat capacities, if given.
!-- By default, use information provided on building type.
    ! CALL get_grid_variable_1d_real( 'albedo_roof', slurb_tile%albedo_roof )
    ! CALL check_grid_variable_1d_real( 'albedo_roof', slurb_tile%albedo_roof, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'dz_roof', slurb_tile%dz_roof )
    ! CALL check_grid_variable_2d_real( 'dz_roof', slurb_tile%dz_roof, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
    ! CALL get_grid_variable_1d_real( 'emiss_roof', slurb_tile%emiss_roof )
    ! CALL check_grid_variable_1d_real( 'emiss_roof', slurb_tile%emiss_roof, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'c_roof', slurb_tile%c_roof )
    ! CALL check_grid_variable_2d_real( 'c_roof', slurb_tile%c_roof, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_roof by dz_roof.

    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_roof(k,i,j) = slurb_tile%c_roof(k,i,j) * slurb_tile%dz_roof(k,i,j)
        enddo
      enddo
    enddo
    ! CALL get_grid_variable_1d_real( 'z0_roof', slurb_tile%z0_roof )
    ! CALL check_grid_variable_1d_real( 'z0_roof', slurb_tile%z0_roof,                                     &
    !                                   TINY( 1.0_field_r ), 0.5_field_r * MINVAL( slurb_tile%z_mo ) )
    ! CALL get_grid_variable_1d_real( 'z0h_roof', slurb_tile%z0h_roof )
    ! CALL check_grid_variable_1d_real( 'z0h_roof', slurb_tile%z0h_roof,                                   &
    !                                   TINY( 1.0_field_r ), 0.5_field_r * MINVAL( slurb_tile%z_mo ) )
    ! CALL get_grid_variable_2d_real( 'lambda_roof', slurb_tile%lambda_roof )
    ! CALL check_grid_variable_2d_real( 'lambda_roof', slurb_tile%lambda_roof,                             &
    !                                   TINY( 1.0_field_r ), HUGE( 1.0_field_r )  )
    ! CALL get_grid_variable_1d_real( 'albedo_wall', slurb_tile%albedo_wall )
    ! CALL check_grid_variable_1d_real( 'albedo_wall', slurb_tile%albedo_wall, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'dz_wall', slurb_tile%dz_wall )
    ! CALL check_grid_variable_2d_real( 'dz_wall', slurb_tile%dz_wall, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
    ! CALL get_grid_variable_1d_real( 'emiss_wall', slurb_tile%emiss_wall )
    ! CALL check_grid_variable_1d_real( 'emiss_wall', slurb_tile%emiss_wall, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'c_wall', slurb_tile%c_wall )
    ! CALL check_grid_variable_2d_real( 'c_wall', slurb_tile%c_wall, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_wall by dz_wall.
    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_wall(k,i,j) = slurb_tile%c_wall(k,i,j) * slurb_tile%dz_wall(k,i,j)
        enddo
      enddo
    enddo
    ! CALL get_grid_variable_1d_real( 'z0_wall', slurb_tile%z0_wall )
    ! CALL check_grid_variable_1d_real( 'z0_wall', slurb_tile%z0_wall, TINY( 1.0_field_r ), 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'lambda_wall', slurb_tile%lambda_wall )
    ! CALL check_grid_variable_2d_real( 'lambda_wall', slurb_tile%lambda_wall,                             &
    !                                   TINY( 1.0_field_r ), HUGE( 1.0_field_r )  )
    ! CALL get_grid_variable_1d_real( 'albedo_window', slurb_tile%albedo_win )
    ! CALL check_grid_variable_1d_real( 'albedo_window', slurb_tile%albedo_win, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'dz_window', slurb_tile%dz_win )
    ! CALL check_grid_variable_2d_real( 'dz_window', slurb_tile%dz_win, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
    ! CALL get_grid_variable_1d_real( 'emiss_window', slurb_tile%emiss_win )
    ! CALL check_grid_variable_1d_real( 'emiss_window', slurb_tile%emiss_win, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'c_window', slurb_tile%c_win )
    ! CALL check_grid_variable_2d_real( 'c_window', slurb_tile%c_win, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity, so multiply c_win
!-- by dz_win.
    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_win(k,i,j) = slurb_tile%c_win(k,i,j) * slurb_tile%dz_win(k,i,j)
        enddo
      enddo
    enddo
    ! slurb_tile%c_win = slurb_tile%c_win * slurb_tile%dz_win
    ! CALL get_grid_variable_2d_real( 'lambda_window', slurb_tile%lambda_win )
    ! CALL check_grid_variable_2d_real( 'lambda_window', slurb_tile%lambda_win,                            &
    !                                   TINY( 1.0_field_r ), HUGE( 1.0_field_r )  )
    ! CALL get_grid_variable_1d_real( 'transmissivity_window', slurb_tile%transmissivity_win )
    ! CALL check_grid_variable_1d_real( 'transmissivity_window', slurb_tile%transmissivity_win,            &
    !                                   0.0_field_r, 1.0_field_r )

    ! CALL get_grid_variable_1d_real( 'window_fraction', slurb_tile%f_win, window_fraction )
    ! CALL check_grid_variable_1d_real( 'window_fraction', slurb_tile%f_win, 0.0_field_r, 1.0_field_r )


    ! CALL get_grid_variable_1d_real( 'albedo_road', slurb_tile%albedo_road )
    ! CALL check_grid_variable_1d_real( 'albedo_road', slurb_tile%albedo_road, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'dz_road', slurb_tile%dz_road )
    ! CALL check_grid_variable_2d_real( 'dz_road', slurb_tile%dz_road, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
    ! CALL get_grid_variable_1d_real( 'emiss_road', slurb_tile%emiss_road )
    ! CALL check_grid_variable_1d_real( 'emiss_road', slurb_tile%emiss_road, 0.0_field_r, 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'c_road', slurb_tile%c_road )
    ! CALL check_grid_variable_2d_real( 'c_road', slurb_tile%c_road, TINY( 1.0_field_r ), HUGE( 1.0_field_r ) )
!
!-- SLUrb uses the total layer heat capacity instead of specific heat capacity,
!-- so multiply c_road by dz_road.

    do j=2,j1
      do i=2,i1
        do k=1,4
            slurb_tile%c_road(k,i,j) = slurb_tile%c_road(k,i,j) * slurb_tile%dz_road(k,i,j)
        enddo
      enddo
    enddo
    ! CALL get_grid_variable_1d_real( 'z0_road', slurb_tile%z0_road )
    ! CALL check_grid_variable_1d_real( 'z0_road', slurb_tile%z0_road, TINY( 1.0_field_r ), 1.0_field_r )
    ! CALL check_grid_variable_1d_real( 'z0h_road', slurb_tile%z0h_road, TINY( 1.0_field_r ), 1.0_field_r )
    ! CALL get_grid_variable_2d_real( 'lambda_road', slurb_tile%lambda_road )
    ! CALL check_grid_variable_2d_real( 'lambda_road', slurb_tile%lambda_road,                             &
    !                                   TINY( 1.0_field_r ), HUGE( 1.0_field_r )  )
!
!-- Compute weighted wall-window albedo.


    call slurb_read_RADnamelist
    if ( namemiss_road /= -9999 ) THEN
        slurb_tile%emiss_road(:,:) = namemiss_road
    ENDIF 
    if ( namemiss_wall /= -9999 ) THEN
        slurb_tile%emiss_wall(:,:) = namemiss_wall
    ENDIF
    if ( namemiss_win /= -9999 ) THEN
        slurb_tile%emiss_win(:,:) = namemiss_win
    ENDIF
    if ( namemiss_roof /= -9999 ) THEN
        slurb_tile%emiss_roof(:,:) = namemiss_roof
    ENDIF
    if ( namalbedo_road /= -9999 ) THEN
        slurb_tile%albedo_road(:,:) = namalbedo_road
    ENDIF
    if ( namalbedo_wall /= -9999 ) THEN
        slurb_tile%albedo_wall(:,:) = namalbedo_wall
    ENDIF
    if ( namalbedo_win /= -9999 ) THEN
        slurb_tile%albedo_win(:,:) = namalbedo_win
    ENDIF
    if ( namalbedo_roof /= -9999 ) THEN
        slurb_tile%albedo_roof(:,:) = namalbedo_roof
    ENDIF
    if ( namlambda_road /= -9999 ) THEN
        slurb_tile%lambda_road(:,:,:) = namlambda_road
    ENDIF
    if ( namlambda_wall /= -9999 ) THEN
        slurb_tile%lambda_wall(:,:,:) = namlambda_wall
    ENDIF
    if ( namlambda_win /= -9999 ) THEN
        slurb_tile%lambda_win(:,:,:) = namlambda_win
    ENDIF
    if ( namlambda_roof /= -9999 ) THEN
        slurb_tile%lambda_roof(:,:,:) = namlambda_roof
    ENDIF


    do j=2,j1
      do i=2,i1
       slurb_tile%albedo_wall_win(i,j) = ( 1.0_field_r - slurb_tile%f_win(i,j) ) * slurb_tile%albedo_wall(i,j) +                &
                                 slurb_tile%f_win(i,j) * slurb_tile%albedo_win(i,j)
      enddo
    enddo
!
!-- Compute the cumulative layer thickness zw for windows.
    do j=2,j1
        do i=2,i1
            slurb_tile%zw_win(nzt_win,i,j) = slurb_tile%dz_win(nzt_win,i,j)
            DO  k = nzt_win+1, nzb_win
                slurb_tile%zw_win(k,i,j) = slurb_tile%zw_win(k-1,i,j) + slurb_tile%dz_win(k,i,j)
            enddo
        enddo
    enddo

 END SUBROUTINE process_surface_parameters



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes SLUrb model variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_slurb_variables
    use modfields, only : thl0, ql0, qt0, u0, v0, exnf
    use modglobal, only : cp, rlv, cu, cv, i1, j1, ep
    use modsurface, only : ps
    REAL(field_r) ::  bc_atm  !< initial atmospheric boundary condition for temperature
    REAL(field_r) ::  e_s     !< initial water vapor saturation pressure
    real du,dv
    integer i,j,k_atm,k_topo

    do j=2,j1
        do i=2,i1

        k_atm = 1 !TODOSELF
        k_topo = 1

        ! TODOSELF is ql0 the correct liquid water?
        ! in PALM pt=liquid water potential temperature
        ! this implies slurb_tile%pt1 = pt+L/cpexn ql0, slurb_tile%q1 = q - ql, vpt1 = pt1 * (1+0.61q1)
        ! write (6,*) 'thl0: ', thl0(i, j, k_atm) 
        ! write (6,*) 'ql0: ', ql0(i,j,k_atm)
        slurb_tile%pt1(i,j)  = thl0(i, j, k_atm) + (rlv/(cp * exnf(k_atm)))  * ql0(i,j,k_atm)
        slurb_tile%q1(i,j)   = qt0(i, j, k_atm) - ql0(i, j, k_atm)
        slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )


        !
        !--    Calculate the pt, vpt and q for atmosphere depending on what modules are enabled.
        ! IF ( bulk_cloud_model )  THEN
        !     slurb_tile%pt1(i,j)  = pt(k_atm,j,i) + lv_d_cp * (1 / exnf(k_atm)) * ql(k_atm,j,i)
        !     slurb_tile%q1(i,j)   = q(k_atm,j,i) - ql(k_atm,j,i)
        !     slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )
        ! ELSEIF ( cloud_droplets )  THEN
        !     slurb_tile%pt1(i,j)  = pt(k_atm,j,i) + lv_d_cp * (1 / exnf(k_atm)) * ql(k_atm,j,i)
        !     slurb_tile%q1(i,j)   = q(k_atm,j,i)
        !     slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )
        ! ELSE
        !     slurb_tile%pt1(i,j) = pt(k_atm,j,i)
        !     IF ( moist_physics )  THEN
        !         slurb_tile%q1(i,j)   = q(k_atm,j,i)
        !         slurb_tile%vpt1(i,j) = slurb_tile%pt1(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q1(i,j) )
        !     ENDIF
        ! ENDIF



        

        ! slurb_tile%uv_abs1(i,j) = SQRT( ( 0.5 * ( u(k_atm,j,i) + u(k_atm,j,i+1) ) )**2 +                    &
        !                         ( 0.5 * ( v(k_atm,j,i) + v(k_atm,j+1,i) ) )**2 )
        du = 0.5*(u0(i,j,1) + u0(i+1,j,1)) + cu
        dv = 0.5*(v0(i,j,1) + v0(i,j+1,1)) + cv
        slurb_tile%uv_abs1(i,j) = sqrt(du**2 + dv**2)
        ! slurb_tile%uv_abs1(i,j) = max(0.1, sqrt(du**2 + dv**2)) DALES VERSION
        slurb_tile%uv_eff1(i,j) = slurb_tile%uv_abs1(i,j)

    !
    !--    Initialize the prognostic model variables and all the other variables where the value from
    !--    the previous time step is used. For restart runs or if spinup data is available, these are
    !--    read from the restart files, and should not be overwritten here.
        ! IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.  .NOT. read_spinup_data )   &
        ! THEN
    !
    !--       If spinup is enabled for current run, use diurnal mean spinup pt as the initial
    !--       atmospheric boundary condition. Otherwise, use the first atmospheric grid level.
        ! IF ( spinup )  THEN
        !     bc_atm = spinup_pt_mean * exnf(k_topo)
        ! ELSE
        bc_atm = slurb_tile%pt1(i,j) * exnf(k_topo)
        ! ENDIF

        slurb_tile%us_urb(i,j)  = 1.0_field_r
        slurb_tile%us_can(i,j)  = 1.0_field_r
        slurb_tile%us_roof(i,j) = 1.0_field_r
        slurb_tile%us_road(i,j) = 1.0_field_r

        slurb_tile%uv_abs_can(i,j) = slurb_tile%uv_abs_can_coef(i,j) * slurb_tile%uv_abs1(i,j)
        slurb_tile%uv_eff_can(i,j) = slurb_tile%uv_abs_can(i,j)

        slurb_tile%shf_urb(i,j)  = 0.0_field_r
        slurb_tile%qsws_urb(i,j) = 0.0_field_r

        slurb_tile%t_can(i,j)   = bc_atm
        slurb_tile%t_can_p(i,j) = slurb_tile%t_can(i,j)

        slurb_tile%t_indoor(i,j) = building_indoor_temperature
        slurb_tile%t_soil(i,j) = deep_soil_temperature
        slurb_tile%shf_external(i,j) = shf_external
        slurb_tile%qsws_external(i,j) = qsws_external

!
!--       For subsurface temps, a steady-state 1D heat equation solution will be used as the
!--       initial temperature profile. This might or might not speed up the spinup process.
!--       In case of windowless facade, set window temps to fill value to prevent meaningless
!--       output values. Vice versa for the opposite case.
        IF ( slurb_tile%f_win(i,j) < 1.0_field_r )  THEN
            slurb_tile%t_wall_a(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_wall_a, 1 ), bc_atm,         &
                                                        slurb_tile%t_indoor(i,j),                         &
                                                        slurb_tile%conductivity_wall(:,i,j) )
            slurb_tile%t_wall_a_p(:,i,j) = slurb_tile%t_wall_a(:,i,j)
            slurb_tile%t_wall_b(:,i,j)   = slurb_tile%t_wall_a(:,i,j)
            slurb_tile%t_wall_b_p(:,i,j) = slurb_tile%t_wall_a(:,i,j)
        ELSE
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%t_wall_a(:,i,j)   = output_fill_value
                slurb_tile%t_wall_a_p(:,i,j) = output_fill_value
            ENDIF
        ENDIF

        IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
            slurb_tile%t_win_a(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_win_a, 1 ), bc_atm,           &
                                                        slurb_tile%t_indoor(i,j),                          &
                                                        slurb_tile%conductivity_win(:,i,j) )
            slurb_tile%t_win_a_p(:,i,j) = slurb_tile%t_win_a(:,i,j)
            slurb_tile%t_win_b(:,i,j)   = slurb_tile%t_win_a(:,i,j)
            slurb_tile%t_win_b_p(:,i,j) = slurb_tile%t_win_a(:,i,j)
        ELSE
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%t_win_a(:,i,j)   = output_fill_value
                slurb_tile%t_win_a_p(:,i,j) = output_fill_value
            ENDIF
        ENDIF

        slurb_tile%t_roof(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_roof, 1 ), bc_atm,                &
                                                    slurb_tile%t_indoor(i,j), slurb_tile%conductivity_roof(:,i,j) )
        slurb_tile%t_roof_p(:,i,j) = slurb_tile%t_roof(:,i,j)

        slurb_tile%t_road(:,i,j) = calc_1d_heat_equation( SIZE( slurb_tile%t_road, 1 ), bc_atm,                &
                                                    slurb_tile%t_soil(i,j), slurb_tile%conductivity_road(:,i,j) )
        slurb_tile%t_road_p(:,i,j) = slurb_tile%t_road(:,i,j)

        IF ( moist_physics )  THEN
            slurb_tile%vpt_can(i,j) = 0.0_field_r

            slurb_tile%q_can(i,j)        = slurb_tile%q1(i,j)
            slurb_tile%q_can_p(i,j)      = slurb_tile%q_can(i,j)
            slurb_tile%m_liq_roof(i,j)   = 0.0_field_r
            slurb_tile%m_liq_roof_p(i,j) = slurb_tile%m_liq_roof(i,j)
            slurb_tile%m_liq_road(i,j)   = 0.0_field_r
            slurb_tile%m_liq_road_p(i,j) = slurb_tile%m_liq_road(i,j)

            slurb_tile%q_roof(i,j) = slurb_tile%q1(i,j)
            slurb_tile%q_road(i,j) = slurb_tile%q_can(i,j)
        ENDIF

        slurb_tile%ol_roof(i,j) = slurb_tile%z_mo(i,j)     / zeta_min
        slurb_tile%ol_road(i,j) = slurb_tile%z_mo_can(i,j) / zeta_min
        slurb_tile%ol_can(i,j)  = slurb_tile%z_mo(i,j)     / zeta_min
        slurb_tile%ol_urb(i,j)  = slurb_tile%z_mo(i,j)     / zeta_min
        ! ENDIF

    !
    !--    Init potential temperatures and virtual potential temperatures. These need to be computed
    !--    also for the restart case, as d_exner is not yet available when rrd routines are called.
        slurb_tile%pt_can(i,j) = slurb_tile%t_can(i,j) / exnf(k_topo)

        IF ( slurb_tile%f_win(i,j) < 1.0_field_r )  THEN
            slurb_tile%pt_wall_a(i,j) = slurb_tile%t_wall_a(nzt_wall,i,j) / exnf(k_topo)
            slurb_tile%pt_wall_b(i,j) = slurb_tile%t_wall_b(nzt_wall,i,j) / exnf(k_topo)
        ELSE
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%pt_wall_a(i,j) = output_fill_value
                slurb_tile%pt_wall_b(i,j) = output_fill_value
            ENDIF
        ENDIF

        IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
            slurb_tile%pt_win_a(i,j) = slurb_tile%t_win_a(nzt_win,i,j) / exnf(k_topo)
            slurb_tile%pt_win_b(i,j) = slurb_tile%t_win_b(nzt_win,i,j) / exnf(k_topo)
        ELSE
            IF ( .NOT. data_output_raw )  THEN
                slurb_tile%pt_win_a(i,j) = output_fill_value
                slurb_tile%pt_win_b(i,j) = output_fill_value
            ENDIF
        ENDIF

        slurb_tile%pt_roof(i,j) = slurb_tile%t_roof(nzt_roof,i,j) / exnf(k_topo)
        slurb_tile%pt_road(i,j) = slurb_tile%t_road(nzt_road,i,j) / exnf(k_topo)

        IF ( moist_physics )  THEN
            slurb_tile%vpt_can(i,j)  = slurb_tile%pt_can(i,j)  * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_can(i,j)  )
            slurb_tile%vpt_roof(i,j) = slurb_tile%pt_roof(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_roof(i,j) )
            slurb_tile%vpt_road(i,j) = slurb_tile%pt_road(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_road(i,j) )
        ENDIF

    !
    !--    Initialize tendencies to zero.
        slurb_tile%tt_can(i,j)      = 0.0_field_r
        slurb_tile%tt_wall_a(:,i,j) = 0.0_field_r
        slurb_tile%tt_wall_b(:,i,j) = 0.0_field_r
        slurb_tile%tt_win_a(:,i,j)  = 0.0_field_r
        slurb_tile%tt_win_b(:,i,j)  = 0.0_field_r
        slurb_tile%tt_roof(:,i,j)   = 0.0_field_r
        slurb_tile%tt_road(:,i,j)   = 0.0_field_r
        IF ( moist_physics )  THEN
            slurb_tile%tq_can(i,j)      = 0.0_field_r
            slurb_tile%tm_liq_roof(i,j) = 0.0_field_r
            slurb_tile%tm_liq_road(i,j) = 0.0_field_r
        ENDIF

    !
    !--    Initialize model variables which are not used prior to an assignment in the model itself.
    !--    Thus these initializations should not end up being used in code, but as this is not
    !--    guaranteed with e.g. future changes, initialize them nevertheless. For the same reason,
    !--    these are not included in the restart data. But if in the future there is an usage prior to
    !--    proper assignment by the model, the respective variable should be added to restart routines,
    !--    and given a proper intialization.
        slurb_tile%albedo_urb(i,j)     = 0.0_field_r
        slurb_tile%emiss_urb(i,j)      = 1.0_field_r
        slurb_tile%rad_lw_in_urb(i,j)  = 0.0_field_r
        slurb_tile%rad_lw_out_urb(i,j) = 0.0_field_r
        slurb_tile%rad_sw_in_urb(i,j)  = 0.0_field_r
        slurb_tile%rad_sw_out_urb(i,j) = 0.0_field_r
        slurb_tile%ram_urb(i,j)        = 1E3_field_r
        slurb_tile%rib_urb(i,j)        = 0.0_field_r
        slurb_tile%t_2m_urb(i,j)       = 0.0_field_r
        slurb_tile%t_c_urb(i,j)        = 0.0_field_r
        slurb_tile%t_h_urb(i,j)        = 0.0_field_r
        slurb_tile%t_rad_urb(i,j)      = 0.0_field_r
        slurb_tile%usws_urb(i,j)       = 0.0_field_r
        slurb_tile%vsws_urb(i,j)       = 0.0_field_r

        slurb_tile%shf_can(i,j)    = 0.0_field_r
        slurb_tile%shf_road(i,j)   = 0.0_field_r
        slurb_tile%shf_roof(i,j)   = 0.0_field_r
        slurb_tile%shf_wall_a(i,j) = 0.0_field_r
        slurb_tile%shf_wall_b(i,j) = 0.0_field_r
        slurb_tile%shf_win_a(i,j)  = 0.0_field_r
        slurb_tile%shf_win_b(i,j)  = 0.0_field_r

        slurb_tile%ghf_road(i,j)   = 0.0_field_r
        slurb_tile%ghf_roof(i,j)   = 0.0_field_r
        slurb_tile%ghf_wall_a(i,j) = 0.0_field_r
        slurb_tile%ghf_wall_b(i,j) = 0.0_field_r
        slurb_tile%ghf_win_a(i,j)  = 0.0_field_r
        slurb_tile%ghf_win_b(i,j)  = 0.0_field_r

        slurb_tile%rad_lw_net_can(i,j)    = 0.0_field_r
        slurb_tile%rad_lw_net_road(i,j)   = 0.0_field_r
        slurb_tile%rad_lw_net_roof(i,j)   = 0.0_field_r
        slurb_tile%rad_lw_net_urb(i,j)    = 0.0_field_r
        slurb_tile%rad_lw_net_wall_a(i,j) = 0.0_field_r
        slurb_tile%rad_lw_net_wall_b(i,j) = 0.0_field_r
        slurb_tile%rad_lw_net_win_a(i,j)  = 0.0_field_r
        slurb_tile%rad_lw_net_win_b(i,j)  = 0.0_field_r
        slurb_tile%rad_sw_in_road(i,j)    = 0.0_field_r
        slurb_tile%rad_sw_in_win_a(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_in_win_b(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_road(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_roof(i,j)   = 0.0_field_r
        slurb_tile%rad_sw_net_urb(i,j)    = 0.0_field_r
        slurb_tile%rad_sw_net_wall_a(i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_wall_b(i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_win_a (i,j) = 0.0_field_r
        slurb_tile%rad_sw_net_win_b (i,j) = 0.0_field_r

        slurb_tile%rib_can(i,j) = 0.0_field_r
        slurb_tile%rib_road(i,j) = 0.0_field_r
        slurb_tile%rib_roof(i,j) = 0.0_field_r

        slurb_tile%rah_can(i,j)  = 1E3_field_r
        slurb_tile%rah_road(i,j) = 1E3_field_r
        slurb_tile%rah_roof    = 1E3_field_r

        IF ( facade_rah_doe )  THEN
            slurb_tile%rah_wall_a(i,j) = 1E3_field_r
            slurb_tile%rah_wall_b(i,j) = 1E3_field_r
            slurb_tile%rah_win_a(i,j)  = 1E3_field_r
            slurb_tile%rah_win_b(i,j)  = 1E3_field_r
        ELSE
            slurb_tile%rah_facade(i,j) = 1E3_field_r
        ENDIF

        IF ( moist_physics )  THEN
            slurb_tile%qsws_can(i,j)      = 0.0_field_r
            slurb_tile%qsws_liq_road(i,j) = 0.0_field_r
            slurb_tile%qsws_liq_roof(i,j) = 0.0_field_r
            slurb_tile%qsws_road(i,j)     = 0.0_field_r
            slurb_tile%qsws_roof(i,j)     = 0.0_field_r

            e_s = 0.01_field_r * magnus( MIN( slurb_tile%t_road(nzt_road,i,j), 333.15_field_r ) )
            slurb_tile%qs_road(i,j) = ep * e_s / ( ps - e_s )
            e_s = 0.01_field_r * magnus( MIN( slurb_tile%t_roof(nzt_roof,i,j), 333.15_field_r ) )
            slurb_tile%qs_roof(i,j) = ep * e_s / ( ps - e_s )

            slurb_tile%c_liq_road(i,j)    = MIN( 1.0_field_r, ( slurb_tile%m_liq_road(i,j) / m_liq_max_road )**0.67 )
            slurb_tile%c_liq_roof(i,j)    = MIN( 1.0_field_r, ( slurb_tile%m_liq_roof(i,j) / m_liq_max_roof )**0.67 )
        ENDIF


        enddo
    enddo

        !
    !-- Calculate logarithms of ratio z/z0.
    !>  TODO: Since the ratios do not change during the simulation, they can be stored once at the
    !>        and stored in surf_slurb, like it is done for the other surface types, too.
            

    do j=2,j1
      do i=2,i1
       ln_z_z0_roof(i,j)  = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_roof(i,j)  )
       ln_z_z0h_roof(i,j) = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0h_roof(i,j) )
       ln_z_z0_urb(i,j)   = LOG( slurb_tile%z_mo(i,j) / slurb_tile%z0_urb(i,j)   )
      enddo
    enddo

 END SUBROUTINE init_slurb_variables

 !--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the magnus formula (Press et al., 1992).
!> The magnus formula is needed to calculate the saturation vapor pressure.
!--------------------------------------------------------------------------------------------------!
 FUNCTION magnus( t )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    REAL(field_r), INTENT(IN) ::  t  !< temperature (K)

    REAL(field_r) ::  magnus

!
!-- Saturation vapor pressure for a specific temperature:
    magnus =  611.2_field_r * EXP( 17.62_field_r * ( t - 273.15_field_r ) / ( t - 29.65_field_r  ) )

 END FUNCTION magnus

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes a steady-state solution for 1D heat equation using Gauss-Seidel iteration.
!> This is used to initialize the material temperatures for roofs, walls, windows and roads,
!> shortening the time required for the spinup. For windows, SW absorption is not considered.
!--------------------------------------------------------------------------------------------------!
 !TODOSELF PURE
 FUNCTION calc_1d_heat_equation( result_size , t_bc_1, t_bc_2, lambda ) RESULT( t_result )

    INTEGER, INTENT(IN) ::  result_size  !< output target size

    REAL(field_r), INTENT(IN) ::  t_bc_1  !< outer t boundary condition
    REAL(field_r), INTENT(IN) ::  t_bc_2  !< inner t boundary condition

    REAL(field_r), DIMENSION(:), INTENT(IN) ::  lambda  !< total layer heat conductivity

    INTEGER ::  ix  !< iteration counter
    INTEGER ::  kx  !< layer running index

    REAL(field_r), PARAMETER ::  omega = 1.0_field_r   !< relaxation to control convergence
    REAL(field_r), PARAMETER ::  tol = 1.0E-6_field_r  !< maximum residual for convergence

    REAL(field_r) ::  res    !< iteration residual for convergence check
    REAL(field_r) ::  t_old  !< previous t of layer for convergence check

    REAL(field_r), DIMENSION(result_size) ::  t_result  !< result t profile
    REAL(field_r), DIMENSION(1:result_size+1) ::  t     !< intermediate t array containing also the BCs


!
!-- Set boundary conditions for the iteration array. The layer against the atmosphere will have a
!-- constant boundary condition and the other boundary is treated similarly to the inner layer
!-- boundary condition in prognostic equations. Thus, an extra layer is neeeded for the inner
!-- temperature array for iteration.
    t(LBOUND( t, 1 )) = t_bc_1
    t(UBOUND( t, 1 )) = t_bc_2

!
!-- Set initial guess for temperature for subsurface layers.
    t(LBOUND( t, 1 )+1:UBOUND( t, 1 )-1) = (t_bc_1 + t_bc_2) / 2.0_field_r
 
!
!-- Gauss-Seidel iteration.
    DO  ix = 1, 1000
       DO  kx = LBOUND( t, 1 )+1, UBOUND( t, 1 )-1
          t_old = t(kx)
          res = 0.0_field_r
          t(kx) = ( lambda(kx) * ( t(kx+1) - t(kx) ) + lambda(kx-1) * ( t(kx-1) - t(kx) ) ) /      &
                  ( lambda(kx) + lambda(kx-1) ) * omega + t(kx)
          res = MAX( res, ABS( t(kx) - t_old ) )
       ENDDO
!
!--    Check for convergence using the stored maximum residual.
       IF ( res < tol )  EXIT
    ENDDO
!
!-- Return the solution.
    t_result(:) = t(LBOUND( t, 1 ):UBOUND( t, 1 )-1)

 END FUNCTION calc_1d_heat_equation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Swap timelevel of the SLUrb model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_swap_timelevel ( mod_count )

    INTEGER, INTENT(IN) ::  mod_count

    SELECT CASE ( mod_count )

       CASE ( 0 )
          slurb_tile%q_can => q_can_1; slurb_tile%q_can_p => q_can_2
          slurb_tile%t_can => t_can_1; slurb_tile%t_can_p => t_can_2
          slurb_tile%m_liq_road => m_liq_road_1; slurb_tile%m_liq_road_p => m_liq_road_2
          slurb_tile%m_liq_roof => m_liq_roof_1; slurb_tile%m_liq_roof_p => m_liq_roof_2
          slurb_tile%t_wall_a => t_wall_a_1; slurb_tile%t_wall_a_p => t_wall_a_2
          slurb_tile%t_wall_b => t_wall_b_1; slurb_tile%t_wall_b_p => t_wall_b_2
          slurb_tile%t_win_a => t_win_a_1; slurb_tile%t_win_a_p => t_win_a_2
          slurb_tile%t_win_b => t_win_b_1; slurb_tile%t_win_b_p => t_win_b_2
          slurb_tile%t_road => t_road_1; slurb_tile%t_road_p => t_road_2
          slurb_tile%t_roof => t_roof_1; slurb_tile%t_roof_p => t_roof_2

       CASE ( 1 )
          slurb_tile%q_can => q_can_2; slurb_tile%q_can_p => q_can_1
          slurb_tile%t_can => t_can_2; slurb_tile%t_can_p => t_can_1
          slurb_tile%m_liq_road => m_liq_road_2; slurb_tile%m_liq_road_p => m_liq_road_1
          slurb_tile%m_liq_roof => m_liq_roof_2; slurb_tile%m_liq_roof_p => m_liq_roof_1
          slurb_tile%t_wall_a => t_wall_a_2; slurb_tile%t_wall_a_p => t_wall_a_1
          slurb_tile%t_wall_b => t_wall_b_2; slurb_tile%t_wall_b_p => t_wall_b_1
          slurb_tile%t_win_a => t_win_a_2; slurb_tile%t_win_a_p => t_win_a_1
          slurb_tile%t_win_b => t_win_b_2; slurb_tile%t_win_b_p => t_win_b_1
          slurb_tile%t_roof => t_roof_2; slurb_tile%t_roof_p => t_roof_1
          slurb_tile%t_road => t_road_2; slurb_tile%t_road_p => t_road_1

    END SELECT

 END SUBROUTINE slurb_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes latent variables which can be inferred from the inputs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE precompute_latent_variables

    use modfields, only: rhof
    use modglobal, only: rhow, rlv, i1, j1, boltz
    use modmpi, only: comm3d

    use modmpi, only: comm3d, mpierr,mpi_min, D_MPI_ALLREDUCE
    
    implicit none

    REAL(field_r) ::  emiss_facade       !< aggregated facade emissivity
    REAL(field_r) ::  f_wall             !< wall fraction
    REAL(field_r) ::  f_win              !< window fraction
    REAL(field_r) ::  wake               !< wake parameter for U_can parametrization in SURFEX
    REAL(field_r) ::  win_nonrefl_1side  !< 1-side nonreflected radiation (for windows)
    REAL(field_r) ::  win_absorp         !< window absorption coefficient
    real :: dt_slurb_individual

    integer i,j,k
!
!-- Precompute model constants.
    rho_lv = rlv * rhof(1)
    drho_l_lv = 1.0_field_r / (rhow * rlv)


    !TODOSELF
    do j=2,j1
        do i=2,i1
            slurb_tile%hw_can(i,j) = 0.3
        enddo
    ENDDO
    
!
!-- Check street canyon height-to-width ratio.
    
    ! do j=2,j1
    !     do i=2,i1
    !         IF ( slurb_tile%h_bld(i,j) == 0.0_field_r  .OR.  slurb_tile%hw_can(i,j) == 0.0_field_r )  THEN
    !             WRITE( message_string, * ) 'Building height or street canyon aspect ratio ' //           &
    !                                         'set to zero at (j,i)=', j, i, 'with non-zero ' //            &
    !                                         'urban_fraction = ', fr_urb(j,i), '.'
    !             CALL message( 'slurb_init', 'SLU0024', 2, 2, myid, 6, 0 )
    !         ENDIF
    !     enddo
    ! ENDDO

!
!-- Precompute layer total conductivities from layer thicknesses and thermal conductivities.
    do j=2,j1
        do i=2,i1
            DO  k = nzt_roof, nzb_roof-1
                slurb_tile%conductivity_roof(k,i,j) = 2.0_field_r / ( slurb_tile%dz_roof(k,i,j)   / slurb_tile%lambda_roof(k,i,j) +   &
                                                        slurb_tile%dz_roof(k+1,i,j) / slurb_tile%lambda_roof(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_roof(nzb_roof,i,j) = 2.0_field_r * slurb_tile%lambda_roof(nzb_roof,i,j) /                &
                                                    slurb_tile%dz_roof(nzb_roof,i,j)

            DO  k = nzt_wall, nzb_wall-1
                slurb_tile%conductivity_wall(k,i,j) = 2.0_field_r / ( slurb_tile%dz_wall(k,i,j)   / slurb_tile%lambda_wall(k,i,j) +   &
                                                        slurb_tile%dz_wall(k+1,i,j) / slurb_tile%lambda_wall(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_wall(nzb_wall,i,j) = 2.0_field_r * slurb_tile%lambda_wall(nzb_wall,i,j) /                &
                                                    slurb_tile%dz_wall(nzb_wall,i,j)

            DO  k = nzt_win, nzb_win-1
                slurb_tile%conductivity_win(k,i,j) = 2.0_field_r / ( slurb_tile%dz_win(k,i,j)   / slurb_tile%lambda_win(k,i,j) +      &
                                                        slurb_tile%dz_win(k+1,i,j) / slurb_tile%lambda_win(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_win(nzb_wall,i,j) = 2.0_field_r * slurb_tile%lambda_wall(nzb_win,i,j) /                  &
                                                slurb_tile%dz_wall(nzb_win,i,j)

        !
        !--    For the road, the last conductance depends on the soil conductance, so we need to
        !--    compute it during time-stepping (as it depends on soil moisture).
            DO  k = nzt_road, nzb_road-1
                slurb_tile%conductivity_road(k,i,j) = 2.0_field_r / ( slurb_tile%dz_road(k,i,j)   / slurb_tile%lambda_road(k,i,j) +   &
                                                        slurb_tile%dz_road(k+1,i,j) / slurb_tile%lambda_road(k+1,i,j) )
            ENDDO
            slurb_tile%conductivity_road(nzb_road,i,j) = 2.0_field_r * slurb_tile%lambda_road(nzb_road,i,j) /                &
                                                    slurb_tile%dz_road(nzb_road,i,j)
        enddo
    ENDDO

!
!-- Precompute sky-view factors.
    slurb_tile%svf_road(:,:) = 0.0_field_r
    slurb_tile%svf_wall(:,:) = 0.0_field_r

    do j=2,j1
        do i=2,i1
            slurb_tile%svf_road(i,j) = SQRT( (slurb_tile%hw_can(i,j))**2 + 1.0_field_r ) - slurb_tile%hw_can(i,j)

            slurb_tile%svf_wall(i,j) = ( slurb_tile%hw_can(i,j) + 1.0_field_r - SQRT( slurb_tile%hw_can(i,j)**2 + 1.0_field_r ) ) /       &
                                ( 2.0_field_r * slurb_tile%hw_can(i,j) )
        enddo
    ENDDO

!
!-- Precompute urban emissivity based on SVFs.
    do j=2,j1
        do i=2,i1
            slurb_tile%emiss_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%emiss_roof(i,j)                                      &
                                + ( 1.0_field_r - slurb_tile%f_bld(i,j) )                                            &
                                    * ( slurb_tile%svf_road(i,j) * slurb_tile%emiss_road(i,j)                            &
                                        + slurb_tile%svf_wall(i,j) * slurb_tile%hw_can(i,j)                                &
                                        * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) * 2.0_field_r * slurb_tile%emiss_wall(i,j)     &
                                            + slurb_tile%f_win(i,j) * 2.0_field_r * slurb_tile%emiss_win(i,j) ) )
        enddo
    ENDDO

!
!-- Preompute the longwave interaction coefficients for surface elements as these are
!-- static in time. Based on Johnson et al. (1991) general formula. Absorption from reflected
!-- radiation is taken into account only after first reflection. The first reflections contribute
!-- around 5% of the total LW budget, while higher order reflections would contribute only <0.5%.
!-- Coefficients are grouped per variable, so they can be effectively used in time-stepping
!-- without wasting too much computational time or memory.
    do j=2,j1
        do i=2,i1
        !
        !--    Compute aggregated facade emissivity for simplification of reflections.
            f_win  = slurb_tile%f_win(i,j)
            f_wall = ( 1.0_field_r - f_win )
            emiss_facade = f_wall * slurb_tile%emiss_wall(i,j) + f_win * slurb_tile%emiss_win(i,j)
        !
        !--    Roof.
        !--    To be multiplied by t_roof**4 in the LW budget:
            slurb_tile%lw_roof_coef(1,i,j) = -slurb_tile%emiss_roof(i,j) * boltz
        !
        !--    To be multiplied by lw_rad_in_urb in the LW budget:
            slurb_tile%lw_roof_coef(2,i,j) = slurb_tile%emiss_roof(i,j)
        !
        !--    Roads.
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_road_coef(1,i,j) = ( - slurb_tile%emiss_road(i,j)                                             &
                                        + slurb_tile%emiss_road(i,j)**2 * ( 1 - emiss_facade )                   &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * slurb_tile%svf_wall(i,j)             &
                                        ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_road_coef(2,i,j) = slurb_tile%emiss_road(i,j) * slurb_tile%svf_road(i,j)                              &
                                    - slurb_tile%emiss_road(i,j) * ( 1.0_field_r - emiss_facade ) * slurb_tile%svf_wall(i,j) &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )
        !
        !--    To be multiplied by (t_wall_a**4 + t_wall_b**4) in the LW budget:
            slurb_tile%lw_road_coef(3,i,j) = ( slurb_tile%emiss_road(i,j) * slurb_tile%emiss_wall(i,j)                          &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + slurb_tile%emiss_road(i,j) * slurb_tile%emiss_wall(i,j)                        &
                                            * ( 1.0_field_r - emiss_facade )                                    &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                       &
                                        ) * 0.5_field_r * f_wall * boltz
            slurb_tile%lw_road_coef(4,i,j) = ( slurb_tile%emiss_road(i,j) * slurb_tile%emiss_win(i,j)                           &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + slurb_tile%emiss_road(i,j) * slurb_tile%emiss_win(i,j)                         &
                                            * ( 1.0_field_r - emiss_facade )                                    &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                       &
                                        ) * 0.5_field_r * f_win * boltz
        !
        !--    Walls.
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_wall_coef(1,i,j) = ( - slurb_tile%emiss_wall(i,j)                                             &
                                        + 0.5_field_r * f_wall * slurb_tile%emiss_wall(i,j)**2                         &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) ) * slurb_tile%svf_wall(i,j)            &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_wall(i,j)**2                                  &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                        ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_wall_coef(2,i,j) = slurb_tile%emiss_wall(i,j) * slurb_tile%svf_wall(i,j)                              &
                                        + slurb_tile%emiss_wall(i,j) * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )             &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                            &
                                        + slurb_tile%emiss_wall(i,j)                                               &
                                        * ( 1.0_field_r - emiss_facade )                                      &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                            &
                                        + slurb_tile%emiss_wall(i,j)                                               &
                                        * ( 1.0_field_r - emiss_facade )                                      &
                                        * slurb_tile%svf_wall(i,j)                                               &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        !
        !--    To be multiplied by t_wall_b**4 in the LW budget:
            slurb_tile%lw_wall_coef(3,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_wall(i,j)**2                          &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                                &
                                        * slurb_tile%svf_wall(i,j)                                               &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                  &
                                            + f_wall * slurb_tile%emiss_wall(i,j)**2                               &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                     &
                                        ) * boltz
        !
        !--    To be multiplied by t_win_a**4 in the LW budget:
            slurb_tile%lw_wall_coef(4,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_wall(i,j)                              &
                                            * slurb_tile%emiss_win(i,j)                                            &
                                            * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                              &
                                            * slurb_tile%svf_wall(i,j)                                             &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                &
                                        + f_win * slurb_tile%emiss_wall(i,j)                                     &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - emiss_facade )                                   &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                   &
                                        ) * boltz
        !
        !--    To be multiplied by t_win_b**4 in the LW budget:
            slurb_tile%lw_wall_coef(5,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_wall(i,j)                              &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                             &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                            * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                               &
                                        + f_win * slurb_tile%emiss_wall(i,j)                                     &
                                            * slurb_tile%emiss_win(i,j)                                           &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                      &
                                        ) * boltz
        !
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_wall_coef(6,i,j) = ( slurb_tile%emiss_wall(i,j) * slurb_tile%emiss_road(i,j)                          &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                        + slurb_tile%emiss_wall(i,j) * slurb_tile%emiss_road(i,j)                        &
                                            * ( 1.0_field_r - emiss_facade )                                   &
                                            * slurb_tile%svf_wall(i,j)                                            &
                                            * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                      &
                                        ) * boltz
        !
        !--    Windows.
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_win_coef(1,i,j) = ( - slurb_tile%emiss_win(i,j)                                               &
                                        + 0.5_field_r * f_win * slurb_tile%emiss_win(i,j)**2                           &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_win * slurb_tile%emiss_win(i,j)**2                                    &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                    ) * boltz
        !
        !--    To be multiplied by lw_rad_in in the LW budget:
            slurb_tile%lw_win_coef(2,i,j) = slurb_tile%emiss_win(i,j) * slurb_tile%svf_wall(i,j)                                &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                                 &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                             &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - emiss_facade )                                       &
                                        * slurb_tile%svf_wall(i,j) * slurb_tile%svf_road(i,j)                             &
                                    + slurb_tile%emiss_win(i,j)                                                 &
                                        * ( 1.0_field_r - emiss_facade )                                       &
                                        * slurb_tile%svf_wall(i,j)                                                &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        !
        !--    To be multiplied by t_win_b**4 in the LW budget:
            slurb_tile%lw_win_coef(3,i,j) = ( 0.5_field_r * f_win * slurb_tile%emiss_win(i,j)**2                             &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_win * slurb_tile%emiss_win(i,j)**2                                    &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        !
        !--    To be multiplied by t_wall_a**4 in the LW budget:
            slurb_tile%lw_win_coef(4,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_win(i,j)                               &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_win(i,j)                                      &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )**2                     &
                                    ) * boltz
        !
        !--    To be multiplied by t_wall_b**4 in the LW budget:
            slurb_tile%lw_win_coef(5,i,j) = ( 0.5_field_r * f_wall * slurb_tile%emiss_win(i,j)                               &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - slurb_tile%emiss_road(i,j) )                               &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                                 &
                                        + f_wall * slurb_tile%emiss_win(i,j)                                      &
                                        * slurb_tile%emiss_wall(i,j)                                            &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        !
        !--    To be multiplied by t_road**4 in the LW budget:
            slurb_tile%lw_win_coef(6,i,j) = ( slurb_tile%emiss_win(i,j) * slurb_tile%emiss_road(i,j)                            &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        + slurb_tile%emiss_win(i,j) * slurb_tile%emiss_road(i,j)                          &
                                        * ( 1.0_field_r - emiss_facade )                                     &
                                        * slurb_tile%svf_wall(i,j)                                              &
                                        * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )                        &
                                    ) * boltz
        enddo
    ENDDO

!
!-- Precompute shortwave radiation reflection denominator.
    do j=2,j1
        do i=2,i1
            slurb_tile%sw_ref_denom(i,j) = 1.0_field_r - slurb_tile%albedo_road(i,j)                                         &
                                    * slurb_tile%albedo_wall_win(i,j)                                             &
                                    * slurb_tile%svf_wall(i,j) * ( 1.0_field_r - slurb_tile%svf_road(i,j) )                    &
                                    - slurb_tile%albedo_wall_win(i,j) * ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) )
        enddo
    ENDDO

!
!-- Compute window layer shortwave absorption based on USM documentation.
!-- @todo This computation needs checking. Now sw_transmitted is not simply equal to
!-- sw_net_win*transmissivity, as 1.0_field_r - SUM(absorption(:,i,j)) != transmissivity(i,j). This is
!-- mitigated for now at the output side. No side effects for the model, as the transmitted
!-- radiation is purely an output.
    do j=2,j1
        do i=2,i1
            win_nonrefl_1side = 1.0 - (slurb_tile%albedo_win(i,j) + slurb_tile%transmissivity_win(i,j)                  &
                                + 1.0_field_r  - SQRT( ( slurb_tile%albedo_win(i,j)                                  &
                                + slurb_tile%transmissivity_win(i,j) + 1.0_field_r )**2                              &
                                - 4.0_field_r * slurb_tile%albedo_win(i,j) ) ) / 2.0_field_r

            win_absorp = -LOG( ( slurb_tile%transmissivity_win(i,j) + slurb_tile%albedo_win(i,j)                        &
                                    - 1.0_field_r + win_nonrefl_1side ) / win_nonrefl_1side                     &
                                ) / slurb_tile%zw_win(nzb_win,i,j)

            DO  k = nzt_win, nzb_win
                IF ( k /= nzt_win)  THEN
        !
        !--          The absorbed fraction is difference between cumulative absorption over the layer.
                    slurb_tile%absorption_win(k,i,j) = win_nonrefl_1side                                          &
                                                * ( EXP( -win_absorp * slurb_tile%zw_win(k-1,i,j) ) -              &
                                                    EXP( -win_absorp * slurb_tile%zw_win(k,i,j)   ) )
                ELSE
        !
        !--          For the first layer, it is the cumulative absorption so far.
                    slurb_tile%absorption_win(k,i,j) = win_nonrefl_1side *                                        &
                                                ( 1.0_field_r - EXP( -win_absorp * slurb_tile%zw_win(k,i,j) ) )
                ENDIF
            ENDDO
        enddo
    ENDDO

!
!-- Coefficient for the canyon wind speed Krayenhoff & Voogt (2007) Eq. (9).
    IF ( uv_can_factor_kray )  THEN
        do j=2,j1
            do i=2,i1
                slurb_tile%uv_abs_can_coef(i,j) = LOG( slurb_tile%h_bld(i,j)  / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /          &
                                    LOG( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) / 3.0_field_r ) / slurb_tile%z0_urb(i,j) ) * &
                                    EXP( -slurb_tile%f_bld_frn(i,j) / ( 2.0_field_r * ( 1.0_field_r - slurb_tile%f_bld(i,j) ) ) )
            enddo
       ENDDO
!
!-- Coefficient for the canyon windspeed as derived in Masson (2000) (original TEB)
    ELSEIF ( uv_can_factor_masson )  THEN
        do j=2,j1
            do i=2,i1
                slurb_tile%uv_abs_can_coef(i,j) = LOG( slurb_tile%h_bld(i,j)  / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /          &
                                    LOG( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) / 3.0_field_r ) / slurb_tile%z0_urb(i,j) ) * &
                                    EXP( -slurb_tile%hw_can(i,j) / 4.0_field_r )
            enddo
       ENDDO
!
!-- Coefficient for the canyon windspeed as implemented in SURFEX v8.1
    ELSEIF ( uv_can_factor_surfex )  THEN
        do j=2,j1
            do i=2,i1
                wake = 1.0_field_r + ( 2.0_field_r / pi - 1.0_field_r ) * 2.0 * ( slurb_tile%hw_can(i,j) - 0.5_field_r )
                wake = MAX( MIN( wake, 1.0_field_r ), 2.0_field_r / pi )
                slurb_tile%uv_abs_can_coef(i,j) = wake * EXP( - slurb_tile%hw_can(i,j) / 4.0_field_r ) *                      &
                                            LOG( 2.0_field_r * slurb_tile%h_bld(i,j) / ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) ) /  &
                                            LOG( ( slurb_tile%z_mo(i,j) + 2.0_field_r * slurb_tile%h_bld(i,j) ) /               &
                                                ( 3.0_field_r * slurb_tile%z0_urb(i,j) ) )
        enddo
       ENDDO
    ELSE
        write(*,*) 'WARNING: no canyon coefficient set for calculating absolute canyon velocity, this might have unintended consequences.'
        call abort
    endif

!
!-- Compute minimum timestep based on SLUrb internal diffusivities.
    do j=2,j1
        do i=2,i1

        !
        !--    Criterion based on subsurface heat diffusion. Heat capacities are already multiplied with
        !--    layer thickness and conductivity is already divided with it. Thus, no need to multiply with
        !--    dz**2 like in usm and lsm. Dimension analysis:
        !--    c [J m-2 K-1] and layer conductivity [W m-2 K-1] -> [J/W] -> [s]
            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_roof(:,i,j) / slurb_tile%conductivity_roof(:,i,j) ) )

            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_road(:,i,j) / slurb_tile%conductivity_road(:,i,j) ) )

            slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                       &
                                    MINVAL( slurb_tile%c_wall(:,i,j) / slurb_tile%conductivity_wall(:,i,j) ) )

            IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
                slurb_tile%dt_max(i,j) = MIN( slurb_tile%dt_max(i,j),                                                    &
                                        MINVAL( slurb_tile%c_win(:,i,j) / slurb_tile%conductivity_win(:,i,j) ) )
            ENDIF
        enddo

    ENDDO
!
!-- Consider a pre-factor (1/8) for the diffusion criterion.
    dt_slurb_individual = MINVAL( slurb_tile%dt_max(2:i1,2:j1) ) * 0.125_field_r
    ! IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL D_MPI_ALLREDUCE( dt_slurb, dt_slurb_individual, 1, mpi_min, comm3d, mpierr )
    write(*,*) 'max_dt for slurb'
    write(*,*) dt_slurb
! #endif

 END SUBROUTINE precompute_latent_variables



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Surface and subsurface energy balance computations of roofs, walls, windows and roads.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_energy_balance_model
   use modglobal, only : i1, j1, cp, rlv, rhow, rk3step, rdt, ep
   use modfields, only : exnf, rhobf, ql0
   use modmicrodata, only : precep, imicro !TODOSELF TEST MICRO
   use modsurface,  only : ps
   implicit none
    INTEGER ::  i       !< loop index (x-direction)
    INTEGER ::  j       !< loop index (y-direction)
    INTEGER ::  k_topo  !< k index of topography
    INTEGER ::  k_atm   !< k index of the first atmospheric level
    INTEGER ::  m       !< running SLUrb tile index

    LOGICAL ::  runge_l  !< flag to vectorize timestep scheme switch

    real :: rho_cp !< cp * rho (J m^-3 K^-1)


    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_energy_balance_model'
    !    CALL debug_message( debug_string, 'start' )
    ! ENDIF

    ! runge_l = ( timestep_scheme(1:5) == 'runge' )
    runge_l = .true.
    k_topo = 1
    k_atm = 1
    rho_cp = cp * rho_air_zw(k_topo) !TODOSELF
   do j=2,j1
      do i=2,i1
      !  k_topo = topo_top_ind(j,i,0)
      !  k_atm = topo_top_ind(j,i,0) + 1


!
!--    Call specific models for all the facets.
       CALL roof_model
       CALL wall_model
       IF ( slurb_tile%f_win(i,j) /= 0.0_field_r )  CALL window_model
       CALL road_model
      enddo
   enddo

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_energy_balance_model'
    !    CALL debug_message( debug_string, 'end' )
    ! ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the new surface prognostic temperature for current time step using RK3.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_surf_t_p ( t, t_p, tt_current, coef_1, coef_2, c )

    REAL(field_r), INTENT(IN) ::  c       !< total layer heat capacity (J m^-2 K^-1)
    REAL(field_r), INTENT(IN) ::  coef_1  !< coefficient A in the prognostic equation (W m^-2)
    REAL(field_r), INTENT(IN) ::  coef_2  !< coefficient B in the prognostic equation (W m^-2 K^-1)
    REAL(field_r), INTENT(IN) ::  t       !< current layer temperature (K)

    REAL(field_r), INTENT(OUT) ::  t_p  !< new layer temperature (K)

    REAL(field_r), INTENT(INOUT) ::  tt_current  !< current temperature tendency (K s^-1)

    REAL(field_r) ::  tt_new  !< new temperature tendency (K s^-1)

!
!-- Compute the prognostic temperature without RK weighting.
    !K = (W m^-2 s + J m^-2 K^-1 K) / (J m^-2 K^-1 + W m^-2 K^-1 s)
    !K = (J m^-2 + J m^-2) / (J m^-2 K^-1 + J m^-2 K^-1)
    t_p = ( coef_1 * rdt * tsc(2) + c * t )  / ( c + coef_2 * rdt * tsc(2) )

!
!-- Compute the RK3 tendency for next time step.
    IF ( c /= 0.0_field_r )  THEN

       t_p    = t_p + rdt * tsc(3) * tt_current
       tt_new = ( t_p - t - rdt * tsc(3) * tt_current ) / ( rdt  * tsc(2) )

       CALL calc_rk3_tend( tt_current, tt_new )

    ENDIF

 END SUBROUTINE calc_surf_t_p


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the new layer prognostic temperature by solving the Fourier diffusion equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_heat_diffusion ( t, t_p, tt_current, c, lambda, t_bc, sw_in, phi )

    REAL(field_r), INTENT(IN) ::  t_bc  !< temperature boundary condition (K)

    REAL(field_r), INTENT(IN), OPTIONAL ::  sw_in  !< incoming shortwave radiation for windows

    REAL(field_r), DIMENSION(:), INTENT(IN) ::  c       !< total heat capacity of the layer (J m^-2 K^-1)
    REAL(field_r), DIMENSION(:), INTENT(IN) ::  lambda  !< total heat conductivity between layers (W m^-2 K^-1)
    REAL(field_r), DIMENSION(:), INTENT(IN) ::  t       !< current time level temperature (K)

    REAL(field_r), DIMENSION(:), INTENT(IN), OPTIONAL ::  phi  !< fraction of incoming shortwave radiation absorbed at window layer

    REAL(field_r), DIMENSION(:), INTENT(OUT) ::  t_p  !< new layer temperature (K)

    REAL(field_r), DIMENSION(:), INTENT(INOUT) ::  tt_current  !< current temperature tendency (K s^-1)

    INTEGER ::  k  !< material layer loop index

    REAL(field_r) ::  tt_new  !<  new temperature tendency (K s^-1)


!
!-- Loop through non-boundary layers of the material.
!-- @todo Split loop into three to move IFs out for better vecotrization.
    DO  k = LBOUND( t, 1 ) + 1, UBOUND( t, 1 )
!
!--    New prognostic layer temperature.
!--    Compute the diffusion between neighbouring layers.
       IF ( k /= UBOUND( t , 1 ) )  THEN
          tt_new = ( 1.0_field_r / c(k) ) * ( lambda(k) * ( t(k+1) - t(k) ) +                           &
                   lambda(k-1) * ( t(k-1) - t(k) ) )
       ELSE
!
!--    Use a constant value boundary condition (skin temperature) for the innermost layer.
            ! (J^-1 m^2 K) (J s^-1 m^-2 K^-1) * K
            ! K s^-1
          tt_new = ( 1.0_field_r / c(k) ) * ( lambda(k) * ( t_bc - t(k) ) +                             &
                   lambda(k-1) * ( t(k-1) - t(k) ) )
       ENDIF
!
!--    Add tendency from absorbed shortwave radiation.
       IF ( PRESENT( sw_in ) )  THEN
          tt_new = tt_new + ( 1.0_field_r / c(k) ) * sw_in * phi(k)
       ENDIF
!
!--    Compute the prognostic temperature and RK3 tendency for next time step.
       t_p(k) = t(k) + rdt * ( tsc(2) * tt_new + tsc(3) * tt_current(k) )

       CALL calc_rk3_tend( tt_current(k), tt_new )
    ENDDO

 END SUBROUTINE calc_heat_diffusion


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the weighted RK3 tendency for the next timestep
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rk3_tend ( tend_current, tend_new )

    REAL(field_r), INTENT(IN) ::  tend_new  !< new RK3 tendency

    REAL(field_r), INTENT(INOUT) ::  tend_current  !< current time level tendency


    IF ( runge_l )  THEN
       IF ( rk3step == 1 )  THEN
          tend_current = tend_new
       ELSEIF ( rk3step < 3 )  THEN !TODOSELF CHECK OF HET WEL KLOPT
          tend_current = -9.5625_field_r * tend_new + 5.3125_field_r * tend_current
       ENDIF
    ENDIF

 END SUBROUTINE calc_rk3_tend


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for roofs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE roof_model

    REAL(field_r) ::  coef_1      !< coefficient A of the prognostic equation
    REAL(field_r) ::  coef_2      !< coefficient B of the prognostic equation
    REAL(field_r) ::  dq_s_dt     !< water vapour mixing ratio tendency
    REAL(field_r) ::  e_s         !< saturation water vapour pressure
    REAL(field_r) ::  e_s_dt      !< saturation water vapour pressure tendency
    REAL(field_r) ::  f_shf       !< factor for the roof sensible heat flux (W m^-2 K^-1)
    REAL(field_r) ::  f_qsws_liq  !< factor for the latent heat flux from/to liquid water reservoir
    REAL(field_r) ::  tm_new      !< new liquid water reservoir tendency


!
!-- Surface sensible heat flux factor.
    f_shf = rho_cp / slurb_tile%rah_roof(i,j)

!
!-- Compute the nominator and denominator coefficients in
!-- the prognostic equation for the moist case.
    IF ( moist_physics )  THEN
!
!--    Computation of factor for the latent heat flux due to
!--    liquid water reservoir evaporation/condensation.
       e_s = ps * slurb_tile%qs_roof(i,j) / ( slurb_tile%qs_roof(i,j) + ep )

!
!--    In case of evaporation, evaporate only for the liquid water coverage area,
!--    in case of condensation, use the total surface.
       IF ( slurb_tile%qs_roof(i,j) > slurb_tile%q1(i,j) )  THEN
          f_qsws_liq = rho_lv * slurb_tile%c_liq_roof(i,j) / slurb_tile%rah_roof(i,j)
       ELSE
          f_qsws_liq = rho_lv / slurb_tile%rah_roof(i,j)
       ENDIF

       e_s_dt = e_s * ( 17.62_field_r / ( slurb_tile%t_roof(nzt_roof,i,j) -  29.65_field_r ) -                       &
                        17.62_field_r * ( slurb_tile%t_roof(nzt_roof,i,j) - 273.15_field_r ) /                       &
                        ( slurb_tile%t_roof(nzt_roof,i,j) - 29.65_field_r )**2                                  &
                      )

       dq_s_dt = ep * e_s_dt / ( ps - e_s_dt )

!
!--    The coefficients for the moist prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_roof(i,j) + slurb_tile%rad_lw_net_roof(i,j)                                  &
                - 3.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**4                     &
                + f_shf * slurb_tile%pt1(i,j)                                                              &
                + f_qsws_liq * ( slurb_tile%q1(i,j) - slurb_tile%qs_roof(i,j)                                      &
                                 + dq_s_dt * slurb_tile%t_roof(nzt_roof,i,j) )                             &
                + slurb_tile%conductivity_roof(nzt_roof,i,j) * slurb_tile%t_roof(nzt_roof+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**3                      &
                + f_shf * (1 / exnf(k_topo))                                                          &
                + f_qsws_liq * dq_s_dt                                                             &
                + slurb_tile%conductivity_roof(nzt_roof,i,j)

    ELSE
!
!-- The coefficients for the dry prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_roof(i,j) + slurb_tile%rad_lw_net_roof(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**4                      &
                + f_shf * slurb_tile%pt1(i,j)                                                              &
                + slurb_tile%conductivity_roof(nzt_roof,i,j) * slurb_tile%t_roof(nzt_roof+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**3                      &
                + f_shf * (1 / exnf(k_topo))                                                          &
                + slurb_tile%conductivity_roof(nzt_roof,i,j)
    ENDIF

    CALL calc_surf_t_p( slurb_tile%t_roof(nzt_roof,i,j), slurb_tile%t_roof_p(nzt_roof,i,j),                        &
                        slurb_tile%tt_roof(nzt_roof,i,j), coef_1, coef_2, slurb_tile%c_roof(nzt_roof,i,j) )

!
!-- Explicit solution of the Fourier heat equation for the subsurface layers.
    CALL calc_heat_diffusion( slurb_tile%t_roof(:,i,j), slurb_tile%t_roof_p(:,i,j), slurb_tile%tt_roof(:,i,j),             &
                              slurb_tile%c_roof(:,i,j), slurb_tile%conductivity_roof(:,i,j), slurb_tile%t_indoor(i,j) )

!
!-- Compute the diagnostic fluxes for the roof surface.
    slurb_tile%ghf_roof(i,j) = slurb_tile%conductivity_roof(nzb_roof,i,j) *                                        &
                       ( slurb_tile%t_roof_p(nzb_roof,i,j) - slurb_tile%t_indoor(i,j) )

    slurb_tile%pt_roof(i,j) = slurb_tile%t_roof_p(nzt_roof,i,j) * (1 / exnf(k_topo))

    slurb_tile%shf_roof(i,j) = -f_shf * ( slurb_tile%pt1(i,j) - slurb_tile%pt_roof(i,j) )

!
!-- Update longwave radiative flux following linearization.
    slurb_tile%rad_lw_net_roof(i,j) = slurb_tile%rad_lw_net_roof(i,j)                                              &
                              + slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**4                &
                              - 4.0_field_r * slurb_tile%lw_roof_coef(1,i,j) * slurb_tile%t_roof(nzt_roof,i,j)**3       &
                              * ( slurb_tile%t_roof(nzt_roof,i,j) - slurb_tile%t_roof_p(nzt_roof,i,j) )

!
!-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    IF ( moist_physics )  THEN
       slurb_tile%qsws_liq_roof(i,j) = -f_qsws_liq * ( slurb_tile%q1(i,j) - slurb_tile%qs_roof(i,j) +                      &
                                               dq_s_dt * slurb_tile%t_roof(nzt_roof,i,j) -                 &
                                               dq_s_dt * slurb_tile%t_roof_p(nzt_roof,i,j)                 &
                                             )

       slurb_tile%qsws_roof(i,j) = slurb_tile%qsws_liq_roof(i,j)
!
!
!--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
!--    is assumed to be drained into the drainage system (liquid water is not conserved).
!--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
       if (imicro == 0 .or. imicro == 1) then
            slurb_tile%qsws_liq_roof(i,j) = slurb_tile%qsws_roof(i,j)
            else
          IF ( slurb_tile%m_liq_roof(i,j) < m_liq_max_roof )  THEN
             slurb_tile%qsws_liq_roof(i,j) = slurb_tile%qsws_roof(i,j) -                                           &
                                     precep(i,j,k_atm) * rhobf(k_atm) * rhow * rlv
          ENDIF

          !todoself even morme assume precipitation
       ENDIF
          !todoself assume precipitation
!
!--    Compute the total latent heat flux.
       slurb_tile%qsws_roof(i,j) = slurb_tile%qsws_roof(i,j) / rlv
!
!--    Compute the prognostic liquid water reservoir.
       tm_new = - slurb_tile%qsws_liq_roof(i,j) * drho_l_lv
       slurb_tile%m_liq_roof_p(i,j) = slurb_tile%m_liq_roof(i,j) +                                                 &
                              rdt * ( tsc(2) * tm_new + tsc(3) * slurb_tile%tm_liq_roof(i,j) )
!
!--    Check if the liquid water reservoir is overfull. If so, drain excess to the
!--    assumed drainage system (water is not conserved here).
       slurb_tile%m_liq_roof_p(i,j) = MIN( slurb_tile%m_liq_roof_p(i,j), m_liq_max_roof )
!
!--    Check for negative water reservoir. @todo store the removed water as runoff for output.
       slurb_tile%m_liq_roof_p(i,j) = MAX( slurb_tile%m_liq_roof_p(i,j), 0.0_field_r )
!
!--    Compute RK3 tendency.
       CALL calc_rk3_tend( slurb_tile%tm_liq_roof(i,j), tm_new )
!
!--    Compute the new liquid water coverage.
       slurb_tile%c_liq_roof(i,j) = MIN( 1.0_field_r, ( slurb_tile%m_liq_roof_p(i,j) / m_liq_max_roof )**0.67 )
!
!--    Compute new saturation mixing ratio.
       e_s = 0.01_field_r * magnus( MIN( slurb_tile%t_roof_p(nzt_roof,i,j), 333.15_field_r ) )
       slurb_tile%qs_roof(i,j) = ep * e_s / ( ps - e_s )
!
!--    Calculate new mixing ratio and vpt at roof surface.
       slurb_tile%q_roof(i,j) = q_surf( slurb_tile%qs_roof(i,j), slurb_tile%rah_roof(i,j), slurb_tile%q1(i,j), f_qsws_liq )
       slurb_tile%vpt_roof(i,j) = slurb_tile%pt_roof(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_roof(i,j) )

    ENDIF

 END SUBROUTINE roof_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for roads.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE road_model

    REAL(field_r) ::  coef_1      !< coefficient A of the prognostic equation
    REAL(field_r) ::  coef_2      !< coefficient B of the prognostic equation
    REAL(field_r) ::  dq_s_dt     !< water vapour mixing ratio tendency
    REAL(field_r) ::  e_s         !< saturation water vapour pressure
    REAL(field_r) ::  e_s_dt      !< saturation water vapour pressure tendency
    REAL(field_r) ::  f_shf       !< factor for the sensible heat flux from roads (W m^-2 K^-1)
    REAL(field_r) ::  f_qsws_liq  !< factor for the latent heat flux from/to liquid water reservoir
    REAL(field_r) ::  tm_new      !< new liquid water reservoir tendency


!
!-- Surface sensible heat flux factor.
    f_shf = rho_cp / slurb_tile%rah_road(i,j)
!
!-- Compute the nominator and denominator coefficients in
!-- the prognostic equation for the moist case.
    IF ( moist_physics )  THEN
!
!--    Computation of factor for the latent heat flux due to
!--    liquid water reservoir evaporation/condensation.
       e_s = ps * slurb_tile%qs_road(i,j) / ( slurb_tile%qs_road(i,j) + ep )

!
!--    In case of evaporation, evaporate only for the liquid water coverage area,
!--    in case of condensation, use the total surface.
       IF ( slurb_tile%qs_road(i,j) > slurb_tile%q_can(i,j) )  THEN
          f_qsws_liq = rho_lv * slurb_tile%c_liq_road(i,j) / slurb_tile%rah_road(i,j)
       ELSE
          f_qsws_liq = rho_lv / slurb_tile%rah_road(i,j)
       ENDIF

       e_s_dt = e_s * ( 17.62_field_r / ( slurb_tile%t_road(nzt_road,i,j) - 29.65_field_r ) -                        &
                        17.62_field_r * ( slurb_tile%t_road(nzt_road,i,j) - 273.15_field_r ) /                       &
                        ( slurb_tile%t_road(nzt_road,i,j) - 29.65_field_r)**2                                   &
                      )

       dq_s_dt = ep * e_s_dt / ( ps - e_s_dt )

!
!--    The coefficients for the moist prognostic equation for temperature. For the longwave balance,
!--    both direct emission and the effect of backreflection are linearized.
       coef_1 = slurb_tile%rad_sw_net_road(i,j) + slurb_tile%rad_lw_net_road(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**4                      &
                + f_shf * slurb_tile%t_can(i,j)                                                            &
                + f_qsws_liq * ( slurb_tile%q_can(i,j) - slurb_tile%qs_road(i,j)                                   &
                                 + dq_s_dt * slurb_tile%t_road(nzt_road,i,j) )                             &
                + slurb_tile%conductivity_road(nzt_road,i,j) * slurb_tile%t_road(nzt_road+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**3                      &
                + f_shf                                                                            &
                + f_qsws_liq * dq_s_dt                                                             &
                + slurb_tile%conductivity_road(nzt_road,i,j)

    ELSE
!
!--    The coefficients for the dry prognostic equation for temperature.
       coef_1 = slurb_tile%rad_sw_net_road(i,j) + slurb_tile%rad_lw_net_road(i,j)                                  &
                -3.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**4                      &
                + f_shf * slurb_tile%t_can(i,j)                                                            &
                + slurb_tile%conductivity_road(nzt_road,i,j) * slurb_tile%t_road(nzt_road+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**3                      &
                + f_shf                                                                            &
                + slurb_tile%conductivity_road(nzt_road,i,j)
    ENDIF

    CALL calc_surf_t_p( slurb_tile%t_road(nzt_road,i,j), slurb_tile%t_road_p(nzt_road,i,j),                        &
                        slurb_tile%tt_road(nzt_road,i,j), coef_1, coef_2, slurb_tile%c_road(nzt_road,i,j) )

!
!-- Heat diffusion through subsurface layers.
    CALL calc_heat_diffusion( slurb_tile%t_road(:,i,j), slurb_tile%t_road_p(:,i,j), slurb_tile%tt_road(:,i,j),             &
                              slurb_tile%c_road(:,i,j), slurb_tile%conductivity_road(:,i,j), slurb_tile%t_soil(i,j) )


    ! write(*,*) "f_shf"
    ! write (*,*) i,j,f_shf
    ! write(*,*) "slurb_tile%t_road_p(nzt_road,i,j)"
    ! write (*,*) i,j,slurb_tile%t_road_p(nzt_road,i,j)
    ! write(*,*) "slurb_tile%t_can(i,j)"
    ! write (*,*) i,j,slurb_tile%t_can(i,j)
    slurb_tile%shf_road(i,j) = -f_shf * ( slurb_tile%t_can(i,j) - slurb_tile%t_road_p(nzt_road,i,j) )

    slurb_tile%pt_road(i,j)  = slurb_tile%t_road_p(nzt_road,i,j) * (1 / exnf(k_topo))

    slurb_tile%ghf_road(i,j) = slurb_tile%conductivity_road(nzb_road,i,j) *                                        &
                       ( slurb_tile%t_road_p(nzb_road,i,j) - slurb_tile%t_soil(i,j) )

!
!-- Update longwave radiative flux following linearization.
    slurb_tile%rad_lw_net_road(i,j) = slurb_tile%rad_lw_net_road(i,j)                                              &
                              + slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**4                &
                              - 4.0_field_r * slurb_tile%lw_road_coef(1,i,j) * slurb_tile%t_road(nzt_road,i,j)**3       &
                              * ( slurb_tile%t_road(nzt_road,i,j) - slurb_tile%t_road_p(nzt_road,i,j) )

!
!-- Compute the water vapor flux from/to liquid water reservoir and the prognostic reservoir level.
    IF ( moist_physics )  THEN
       slurb_tile%qsws_liq_road(i,j) = -f_qsws_liq * ( slurb_tile%q_can(i,j) - slurb_tile%qs_road(i,j) +                   &
                                               dq_s_dt * slurb_tile%t_road(nzt_road,i,j) -                 &
                                               dq_s_dt * slurb_tile%t_road_p(nzt_road,i,j)                 &
                                             )

       slurb_tile%qsws_road(i,j) = slurb_tile%qsws_liq_road(i,j)
!
!--    Modification due to precipitiation. If the liquid reservoir is full, the liquid water
!--    is assumed to be drained into the drainage system (liquid water is not conserved).
!--    The precipitation flux is not included in the surface-atmosphere latent heat flux (qsws).
    !    IF ( precipitation )  THEN

    !    ENDIF
       if (imicro == 0 .or. imicro == 1) then
            slurb_tile%qsws_liq_roof(i,j) = slurb_tile%qsws_roof(i,j)
            else
          IF ( slurb_tile%m_liq_road(i,j) < m_liq_max_road )  THEN
             slurb_tile%qsws_road(i,j) = slurb_tile%qsws_road(i,j) -                                               &
                                 precep(i,j,k_atm) * rhobf(k_atm) * rhow * rlv
          ENDIF

          !todoself even morme assume precipitation
       ENDIF
       ! liquid water reservoir is in m^3/m^2, rain rate in m/s (m^3/m^2 /s)
!
!--    Compute the total latent heat flux.
       slurb_tile%qsws_road(i,j) = slurb_tile%qsws_road(i,j) / rlv
!
!--    Compute the prognostic liquid water reservoir.
       tm_new = - slurb_tile%qsws_liq_road(i,j) * drho_l_lv
       slurb_tile%m_liq_road_p(i,j) = slurb_tile%m_liq_road(i,j) +                                                 &
                              rdt * ( tsc(2) * tm_new + tsc(3) * slurb_tile%tm_liq_road(i,j) )
!
!--    Check if the liquid water reservoir is overfull. If so, drain excess to the
!--    assumed drainage system (water is not conserved here).
       slurb_tile%m_liq_road_p(i,j) = MIN( slurb_tile%m_liq_road_p(i,j), m_liq_max_road )
!
!--    Check for negative water reservoir. Should we adjust qsws_road accordingly?
       slurb_tile%m_liq_road_p(i,j) = MAX( slurb_tile%m_liq_road_p(i,j), 0.0_field_r )
!
!--    Compute RK3 tendency
       CALL calc_rk3_tend( slurb_tile%tm_liq_road(i,j), tm_new )
!
!--    Compute the new liquid water coverage.
       slurb_tile%c_liq_road(i,j) = MIN( 1.0_field_r, ( slurb_tile%m_liq_road_p(i,j) / m_liq_max_road )**0.67 )
!
!--    Compute new saturation mixing ratio.
       e_s = 0.01_field_r * magnus( MIN( slurb_tile%t_road_p(nzt_road,i,j), 333.15_field_r ) )
       slurb_tile%qs_road(i,j) = ep * e_s / ( ps - e_s )
!
!--    Calculate new mixing ratio and vpt at road surface.
       slurb_tile%q_road(i,j) = q_surf( slurb_tile%qs_road(i,j), slurb_tile%rah_road(i,j), slurb_tile%q1(i,j), f_qsws_liq )
       slurb_tile%vpt_road(i,j) = slurb_tile%pt_road(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_road(i,j) )

    ENDIF

 END SUBROUTINE road_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance and subsurface heat diffusion for both walls.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wall_model

   REAL(field_r) ::  coef_1   !< coefficient A of the prognostic equation
   REAL(field_r) ::  coef_2   !< coefficient B of the prognostic equation
   REAL(field_r) ::  f_shf_a  !< factor for the wall surface heat flux (W m^-2 K^-1)
   REAL(field_r) ::  f_shf_b  !< factor for the wall surface heat flux (W m^-2 K^-1)


    IF ( facade_rah_doe )  THEN
       f_shf_a = rho_cp / slurb_tile%rah_wall_a(i,j)
       IF ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = rho_cp / slurb_tile%rah_wall_b(i,j)
    ELSE
       f_shf_a = rho_cp / slurb_tile%rah_facade(i,j)
       IF ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = f_shf_a
    ENDIF

!
!-- The coefficients for the moist prognostic equation for temperature. For the longwave balance,
!-- both direct emission and the effect of backreflection are linearized. The linearization depends
!-- if the canyon is isotropic or not as an average backreflection is used for isotropic canyons.
!-- We consider the walls are dry in all cases, so moist physical processes are not considered.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
       coef_1 = slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_a(i,j)                              &
                - 3.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4                   &
                + f_shf_a * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_a(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**3                    &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       CALL calc_surf_t_p(slurb_tile%t_wall_a(nzt_wall,i,j), slurb_tile%t_wall_a_p(nzt_wall,i,j),                  &
                          slurb_tile%tt_wall_a(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )

       coef_1 = slurb_tile%rad_sw_net_wall_b(i,j) + slurb_tile%rad_lw_net_wall_b(i,j)                              &
                - 3.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4                   &
                + f_shf_b * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_b(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**3                    &
                + f_shf_b                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       CALL calc_surf_t_p( slurb_tile%t_wall_b(nzt_wall,i,j), slurb_tile%t_wall_b_p(nzt_wall,i,j),                 &
                           slurb_tile%tt_wall_b(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )
    ELSE
!
!--    In case of isotropic canyon, wall A and B temperatures are averaged, and thus the prognostic
!--    equation for t_wall_a is representative of both of the walls. Thus both the terms for
!--    t_wall_a as well as for t_wall_b in the longwave radiation balance has a dependency on
!--    the surface temperature.
       coef_1 = slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_a(i,j)                              &
                - 3.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )                     &
                   * slurb_tile%t_wall_a(nzt_wall,i,j)**4                                                  &
                + f_shf_a * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j) * slurb_tile%t_wall_a(nzt_wall+1,i,j)

       coef_2 = -4.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )                      &
                   * slurb_tile%t_wall_a(nzt_wall,i,j)**3                                                  &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_wall(nzt_wall,i,j)

       CALL calc_surf_t_p( slurb_tile%t_wall_a(nzt_wall,i,j), slurb_tile%t_wall_a_p(nzt_wall,i,j),                 &
                           slurb_tile%tt_wall_a(nzt_wall,i,j), coef_1, coef_2, slurb_tile%c_wall(nzt_wall,i,j) )
    ENDIF

    slurb_tile%pt_wall_a(i,j)  = slurb_tile%t_wall_a_p(nzt_wall,i,j) * (1 / exnf(k_topo))
    slurb_tile%shf_wall_a(i,j) = -f_shf_a * ( slurb_tile%t_can(i,j) - slurb_tile%t_wall_a_p(nzt_wall,i,j) )

!
!-- Heat diffusion through subsurface layers.
    CALL calc_heat_diffusion( slurb_tile%t_wall_a(:,i,j), slurb_tile%t_wall_a_p(:,i,j), slurb_tile%tt_wall_a(:,i,j),       &
                              slurb_tile%c_wall(:,i,j), slurb_tile%conductivity_wall(:,i,j), slurb_tile%t_indoor(i,j) )

    slurb_tile%ghf_wall_a(i,j) = slurb_tile%conductivity_wall(nzb_wall,i,j) *                                      &
                         ( slurb_tile%t_wall_a_p(nzb_wall,i,j) - slurb_tile%t_indoor(i,j) )

!
!-- Same treatment for wall B if this is an anisotropic canyon, otherwise copy.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
       slurb_tile%pt_wall_b(i,j)  = slurb_tile%t_wall_b_p(nzt_wall,i,j) * (1 / exnf(k_topo))
       slurb_tile%shf_wall_b(i,j) = -f_shf_b * ( slurb_tile%t_can(i,j) - slurb_tile%t_wall_b_p(nzt_wall,i,j) )

       CALL calc_heat_diffusion( slurb_tile%t_wall_b(:,i,j), slurb_tile%t_wall_b_p(:,i,j), slurb_tile%tt_wall_b(:,i,j),    &
                                 slurb_tile%c_wall(:,i,j), slurb_tile%conductivity_wall(:,i,j), slurb_tile%t_indoor(i,j) )

       slurb_tile%ghf_wall_b(i,j) = slurb_tile%conductivity_wall(nzb_wall,i,j) *                                   &
                            ( slurb_tile%t_wall_b_p(nzb_wall,i,j) - slurb_tile%t_indoor(i,j) )

!
!--    Update longwave radiative fluxes following linearization.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)                                       &
                                   + slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4         &
                                   - 4.0_field_r * slurb_tile%lw_wall_coef(1,i,j)                               &
                                      * slurb_tile%t_wall_a(nzt_wall,i,j)**3                               &
                                   * ( slurb_tile%t_wall_a(nzt_wall,i,j) - slurb_tile%t_wall_a_p(nzt_wall,i,j) )

       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_b(i,j)                                       &
                                   + slurb_tile%lw_wall_coef(1,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4         &
                                   - 4.0_field_r * slurb_tile%lw_wall_coef(1,i,j)                               &
                                      * slurb_tile%t_wall_b(nzt_wall,i,j)**3                               &
                                   * ( slurb_tile%t_wall_b(nzt_wall,i,j) - slurb_tile%t_wall_b_p(nzt_wall,i,j) )
    ELSE
!
!--    Copy all layers including the surface for wall B.
       slurb_tile%t_wall_b_p(:,i,j) = slurb_tile%t_wall_a_p(:,i,j)
       slurb_tile%tt_wall_b(:,i,j)  = slurb_tile%tt_wall_a(:,i,j)
       slurb_tile%pt_wall_b(i,j)    = slurb_tile%pt_wall_a(i,j)
       slurb_tile%shf_wall_b(i,j)   = slurb_tile%shf_wall_a(i,j)
       slurb_tile%ghf_wall_b(i,j)   = slurb_tile%ghf_wall_a(i,j)

!
!--    For longwave radiative flux, we need to add terms for both the t_wall_a and the t_wall_b
!--    in the longwave balance, thus coefficients 1 and 3 are summed here.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)                                       &
                                   + ( slurb_tile%lw_wall_coef(1,i,j) + slurb_tile%lw_wall_coef(3,i,j) )           &
                                      * slurb_tile%t_wall_a(nzt_wall,i,j)**4                               &
                                   - 4.0_field_r * ( slurb_tile%lw_wall_coef(1,i,j)                             &
                                                + slurb_tile%lw_wall_coef(3,i,j) )                         &
                                   * slurb_tile%t_wall_a(nzt_wall,i,j)**3                                  &
                                   * ( slurb_tile%t_wall_a(nzt_wall,i,j) - slurb_tile%t_wall_a_p(nzt_wall,i,j) )
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)
    ENDIF

 END SUBROUTINE wall_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Models the surface energy balance, SW transmission and subsurface heat diffusion for windows.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE window_model

    REAL(field_r) ::  coef_1   !< coefficient A of the prognostic equation
    REAL(field_r) ::  coef_2   !< coefficient B of the prognostic equation
    REAL(field_r) ::  f_shf_a  !< factor for the window surface heat flux
    REAL(field_r) ::  f_shf_b  !< factor for the window surface heat flux


    IF ( facade_rah_doe )  THEN
       f_shf_a = rho_cp / slurb_tile%rah_win_a(i,j)
       IF ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = rho_cp / slurb_tile%rah_win_b(i,j)
    ELSE
       f_shf_a = rho_cp / slurb_tile%rah_facade(i,j)
       IF ( slurb_tile%anisotropic_canyon(i,j) )  f_shf_b = f_shf_a
    ENDIF

!
!-- Computation of the prognostic equation similarly to the walls, with exception of added
!-- shortwave transmission component for surface and subsurface layers. Explanatory comments
!-- are not repeated from the wall model, comments reflect differences specific to windows.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
!
!--    For windows, some of the incoming shortwave radiation is transmitted through the material.
       coef_1 = slurb_tile%rad_sw_net_win_a(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_a(i,j)                                                         &
                - 3.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4                      &
                + f_shf_a * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_a(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**3                       &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       CALL calc_surf_t_p( slurb_tile%t_win_a(nzt_win,i,j), slurb_tile%t_win_a_p(nzt_win,i,j),                     &
                           slurb_tile%tt_win_a(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )

       coef_1 = slurb_tile%rad_sw_net_win_b(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_b(i,j)                                                         &
                - 3.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4                      &
                + f_shf_b * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_b(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**3                       &
                + f_shf_b                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       CALL calc_surf_t_p( slurb_tile%t_win_b(nzt_win,i,j), slurb_tile%t_win_b_p(nzt_win,i,j),                     &
                           slurb_tile%tt_win_b(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )

    ELSE
       coef_1 = slurb_tile%rad_sw_net_win_a(i,j) * slurb_tile%absorption_win(nzt_win,i,j)                          &
                + slurb_tile%rad_lw_net_win_a(i,j)                                                         &
                - 3.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )                       &
                   * slurb_tile%t_win_a(nzt_win,i,j)**4                                                    &
                + f_shf_a * slurb_tile%t_can(i,j)                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j) * slurb_tile%t_win_a(nzt_win+1,i,j)

       coef_2 = -4.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )                        &
                   * slurb_tile%t_win_a(nzt_win,i,j)**3                                                    &
                + f_shf_a                                                                          &
                + slurb_tile%conductivity_win(nzt_win,i,j)

       CALL calc_surf_t_p( slurb_tile%t_win_a(nzt_win,i,j), slurb_tile%t_win_a_p(nzt_win,i,j),                     &
                           slurb_tile%tt_win_a(nzt_win,i,j), coef_1, coef_2, slurb_tile%c_win(nzt_win,i,j) )
    ENDIF

    slurb_tile%pt_win_a(i,j)  = slurb_tile%t_win_a_p(nzt_win,i,j) * (1 / exnf(k_topo))
    slurb_tile%shf_win_a(i,j) = -f_shf_a * ( slurb_tile%t_can(i,j) - slurb_tile%t_win_a_p(nzt_win,i,j) )

!
!-- The transmitted shortwave radiation is included also in the prognostic equations for material
!-- subsurface temperatures.
    CALL calc_heat_diffusion( slurb_tile%t_win_a(:,i,j), slurb_tile%t_win_a_p(:,i,j),                              &
                              slurb_tile%tt_win_a(:,i,j), slurb_tile%c_win(:,i,j),                                 &
                              slurb_tile%conductivity_win(:,i,j), slurb_tile%t_indoor(i,j),                        &
                              slurb_tile%rad_sw_net_win_a(i,j), slurb_tile%absorption_win(:,i,j) )

    slurb_tile%ghf_win_a(i,j) = slurb_tile%conductivity_win(nzb_win,i,j) *                                         &
                        ( slurb_tile%t_win_a_p(nzb_win,i,j) - slurb_tile%t_indoor(i,j) )

    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
       slurb_tile%pt_win_b(i,j)  = slurb_tile%t_win_b_p(nzt_win,i,j) * (1 / exnf(k_topo))
       slurb_tile%shf_win_b(i,j) = -f_shf_b * ( slurb_tile%t_can(i,j) - slurb_tile%t_win_b_p(nzt_win,i,j) )

       CALL calc_heat_diffusion( slurb_tile%t_win_b(:,i,j), slurb_tile%t_win_b_p(:,i,j),                           &
                                 slurb_tile%tt_win_b(:,i,j), slurb_tile%c_win(:,i,j),                              &
                                 slurb_tile%conductivity_win(:,i,j), slurb_tile%t_indoor(i,j),                     &
                                 slurb_tile%rad_sw_net_win_b(i,j), slurb_tile%absorption_win(:,i,j) )

       slurb_tile%ghf_win_b(i,j) = slurb_tile%conductivity_win(nzb_win,i,j) *                                      &
                           ( slurb_tile%t_win_b_p(nzb_win,i,j) - slurb_tile%t_indoor(i,j) )

       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j)                                         &
                                  + slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4             &
                                  - 4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**3    &
                                     * ( slurb_tile%t_win_a(nzt_win,i,j) - slurb_tile%t_win_a_p(nzt_win,i,j) )

       slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%rad_lw_net_win_b(i,j)                                         &
                                  + slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4             &
                                  - 4.0_field_r * slurb_tile%lw_win_coef(1,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**3    &
                                     * ( slurb_tile%t_win_b(nzt_win,i,j) - slurb_tile%t_win_b_p(nzt_win,i,j) )
    ELSE
       slurb_tile%t_win_b_p(:,i,j) = slurb_tile%t_win_a_p(:,i,j)
       slurb_tile%tt_win_b(:,i,j)  = slurb_tile%tt_win_a(:,i,j)
       slurb_tile%pt_win_b(i,j)    = slurb_tile%pt_win_a(i,j)
       slurb_tile%shf_win_b(i,j)   = slurb_tile%shf_win_a(i,j)
       slurb_tile%ghf_win_b(i,j)   = slurb_tile%ghf_win_a(i,j)

       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j)                                         &
                                  + ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )              &
                                     * slurb_tile%t_win_a(nzt_win,i,j)**4                                  &
                                  - 4.0_field_r * ( slurb_tile%lw_win_coef(1,i,j) + slurb_tile%lw_win_coef(3,i,j) )     &
                                  * slurb_tile%t_win_a(nzt_win,i,j)**3                                     &
                                  * ( slurb_tile%t_win_a(nzt_win,i,j) - slurb_tile%t_win_a_p(nzt_win,i,j) )
       slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%rad_lw_net_win_a(i,j)
    ENDIF

 END SUBROUTINE window_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface mixing ratio using resistance weighting.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION q_surf( q_s, rah, q_a, f_qsws )

    REAL(field_r), INTENT(IN) ::  f_qsws  !< factor for the latent heat flux
    REAL(field_r), INTENT(IN) ::  q_a     !< mixing ratio of adjacent air
    REAL(field_r), INTENT(IN) ::  q_s     !< saturation mixing ratio at the surface
    REAL(field_r), INTENT(IN) ::  rah     !< aerodynamic resistance for heat (and for water vapor)

    REAL(field_r) ::  q_surf  !< mixing ratio for the surface
    REAL(field_r) ::  res     !< total surface resistance


!
!-- Total surface resistance.
    res = rah / ( rah + ABS( rho_lv / ( f_qsws + 1.0E-20_field_r ) - rah ) )

!
!-- Use the newly calculated total surface resistance to compute weighted surface mixing ratio.
    ! IF ( bulk_cloud_model )  THEN
!
!--    Assume equal liquid water content in canyon as in air above.
       q_surf = res * q_s + ( 1.0_field_r - res ) * ( q_a - ql0(i,j,k_atm) )
    ! ELSE
    !    q_surf = res * q_s + ( 1.0_field_r - res ) * q_a
    ! ENDIF

 END FUNCTION q_surf

 END SUBROUTINE slurb_energy_balance_model

 !--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the dynamical conditions (wind speed, pt, q, vpt) in the street canyon.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE slurb_canyon_model
    use modglobal, only : i1, j1, cp, rlv, rhow, rk3step, rdt, ep
    INTEGER ::  i       !< loop index (x-direction)
    INTEGER ::  j       !< loop index (y-direction)
    INTEGER ::  k_topo  !< k index of topography
    INTEGER ::  k_atm   !< k index of the first atmospheric level
    INTEGER ::  m       !< running SLUrb tile index

    LOGICAL  ::  runge_l  !< timestep scheme switch for vectorization

    REAL(field_r) ::  c          !< total heat capacity of canyon air column per square metre (J K^-1 m^-2)
    REAL(field_r) ::  coef_1     !< coefficient A for the prognostic equation (W m^-2)
    REAL(field_r) ::  coef_2     !< coefficient B for the prognostic equation (W m^-2 K^-1)
    REAL(field_r) ::  f_shf      !< factor for the sensible heat flux  (W m^-2 K^-1)
    REAL(field_r) ::  f_qsws     !< factor for the latent heat flux
    REAL(field_r) ::  qsws_surf  !< aggregated latent heat flux from canyon surfaces per unit area (W m^-2)
    REAL(field_r) ::  shf_surf   !< aggregated sensible heat flux from canyon surfaces per unit area (W m^-2)
    REAL(field_r) ::  tq_new     !< mixing ratio tendency for the new RK3 time step
    REAL(field_r) ::  tt_new     !< temperature tendency for the new RK3 time step
    REAL(field_r) ::  vtws       !< buoyancy flux (m K s^-1)
    REAL(field_r) ::  ws         !< free-convection scale (m/s)
    real :: rho_cp !< cp * rho (J m^-3 K^-1)
    k_topo = 1
    k_atm = 1
    rho_cp = cp * rho_air_zw(k_topo)

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_canyon_model'
    !    CALL debug_message( debug_string, 'start' )
    ! ENDIF

    ! runge_l = ( timestep_scheme(1:5) == 'runge' )
    runge_l = .true.

        do j=2,j1
      do i=2,i1

   !
   !--    Index offset of surface element point with respect to adjoining atmospheric grid point.
         !  k_topo = topo_top_ind(j,i,0)
         !  k_atm  = topo_top_ind(j,i,0) + 1


         f_shf  = rho_cp / slurb_tile%rah_can(i,j)
   !
   !--    Consider total air mass column within the street canyon.
         c = rho_cp * slurb_tile%h_bld(i,j)
        !  write(*,*) "slurb_tile%hw_can(i,j)"
        !  write (*,*) i,j,slurb_tile%hw_can(i,j)
        !  write(*,*) "slurb_tile%f_win(i,j)"
        !  write (*,*) i,j,slurb_tile%f_win(i,j)
        !  write(*,*) "slurb_tile%shf_road(i,j)"
        !  write (*,*) i,j,slurb_tile%shf_road(i,j)
        !  write(*,*) "slurb_tile%shf_wall_a(i,j)"
        !  write (*,*) i,j,slurb_tile%shf_wall_a(i,j)
        !  write(*,*) "slurb_tile%shf_wall_b(i,j)"
        !  write (*,*) i,j,slurb_tile%shf_wall_b(i,j)
        !  write(*,*) "slurb_tile%shf_win_a(i,j)"
        !  write (*,*) i,j,slurb_tile%shf_win_a(i,j)
        !  write(*,*) "slurb_tile%shf_win_b(i,j)"
        !  write (*,*) i,j,slurb_tile%shf_win_b(i,j)
   !
   !--    In canyon temperature prognostic equation, we use already computed fluxes from surfaces
   !--    in order to ensure consistency and conservation of energy. Thus, only the fluxes between
   !--    canyon air and the atmosphere are linearized.
   !
   !--    Aggregated sensible heat flux from canyon surfaces (per unit area).
         shf_surf = slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                                  &
                                       ( slurb_tile%shf_wall_a(i,j) + slurb_tile%shf_wall_b(i,j) ) +                 &
                                       slurb_tile%f_win(i,j) * ( slurb_tile%shf_win_a(i,j) + slurb_tile%shf_win_b(i,j) )     &
                                    ) + slurb_tile%shf_road(i,j)
   !
   !--    Aggregated flux doesn't contain c_p yet.
         shf_surf = shf_surf

        !  write(*,*) "f_shf"
        !  write (*,*) i,j,f_shf
        !  write(*,*) "slurb_tile%pt1(i,j)"
        !  write (*,*) i,j,slurb_tile%pt1(i,j)
        !  write(*,*) "shf_surf"
        !  write (*,*) i,j,shf_surf
   !
   !--    Coefficients for the prognostic equation of street canyon temperature.
         coef_1 = f_shf * slurb_tile%pt1(i,j) + (shf_surf)
         coef_2 = f_shf * (1 / exnf(k_topo))

        !  write(*,*) "coef_1"
        !  write (*,*) i,j,coef_1
        !  write(*,*) "coef_2"
        !  write (*,*) i,j,coef_2
        ! K = (w m^-2 s + J k^-1 m^-2 K) / (J K^-1 m^-2 + W m^-2 K^-1 s)
        ! K = (J m^-2 + J m^-2) / (J K^-1 m^-2 + J m^-2 K^-1)
         slurb_tile%t_can_p(i,j) = ( coef_1 * rdt * tsc(2) + c * slurb_tile%t_can(i,j) ) /                         &
                           ( c + coef_2 * rdt * tsc(2) )

         slurb_tile%t_can_p(i,j) = slurb_tile%t_can_p(i,j) + rdt * tsc(3) * slurb_tile%tt_can(i,j)

         tt_new = ( slurb_tile%t_can_p(i,j) - slurb_tile%t_can(i,j) - rdt * tsc(3) * slurb_tile%tt_can(i,j) ) /            &
                  ( rdt  * tsc(2) )
   !
   !--    Compute the weighted RK3 tendency to be used in next time step.
         IF ( runge_l )  THEN
            IF ( rk3step == 1 )  THEN
               slurb_tile%tt_can(i,j) = tt_new
            ELSEIF ( rk3step < 3 )  THEN
               slurb_tile%tt_can(i,j) = -9.5625_field_r * tt_new + 5.3125_field_r * slurb_tile%tt_can(i,j)
            ENDIF
         ENDIF

        ! write(*,*) "t_can_p"
        !  write (*,*) i,j,slurb_tile%t_can_p(i,j)
        ! write(*,*) "exnf(k_topo)"
        !  write(*,*) exnf(k_topo)
        !  write(*,*) "pt1"
        !  write (*,*) i,j,slurb_tile%pt1(i,j)
        !  write(*,*) "f_shf"
        !  write (*,*) i,j,f_shf
         !
         !--    Calculate new pt and shf from canyon to atmosphere.
         slurb_tile%pt_can(i,j) = slurb_tile%t_can_p(i,j) * (1 / exnf(k_topo))
         slurb_tile%shf_can(i,j) = -f_shf * ( slurb_tile%pt1(i,j) - slurb_tile%pt_can(i,j) )
         
        !  write(*,*) "pt_can"
        !  write (*,*) i,j,slurb_tile%pt_can(i,j)
   !
   !--    Compute prognostic street canyon mixing ratio.
         IF ( moist_physics )  THEN

            f_qsws = rho_lv / slurb_tile%rah_can(i,j)
   !
   !--       Here our "latent heat capacity" is the canyon air column total mass.
            c = rho_lv * slurb_tile%h_bld(i,j)
   !
   !--       Same for the latent heat flux. Currently only the roads, walls are always dry.
   !--       This is a placeholder aggregation for street canyon vegetation,
   !--       e.g. green walls, low vegetation etc.
            qsws_surf = slurb_tile%qsws_road(i,j)
   !
   !--       Aggregated flux doesn't contain l_v yet.
            qsws_surf = slurb_tile%qsws_road(i,j) * rlv
   !
   !--       Compute new prognostic canyon mixing ratio.
            coef_1 = f_qsws * slurb_tile%q1(i,j) + qsws_surf
            coef_2 = f_qsws

            slurb_tile%q_can_p(i,j) = ( coef_1 * rdt * tsc(2) + c * slurb_tile%q_can(i,j) ) /                      &
                              ( c + coef_2 * rdt * tsc(2) )

            slurb_tile%q_can_p(i,j) = slurb_tile%q_can_p(i,j) + rdt * tsc(3) * slurb_tile%tq_can(i,j)
   !
   !--       Prevent negative mixing ratios due to temporal discretization. This is done before
   !--       the computation of tq_new in order to conserve energy.
            ! write(*,*) "q_can_p"
            ! write(*,*) i,j,slurb_tile%q_can_p(i,j)
            IF ( slurb_tile%q_can_p(i,j) < 0.0_field_r )  slurb_tile%q_can_p(i,j) = 0.0_field_r

            tq_new = ( slurb_tile%q_can_p(i,j) - slurb_tile%q_can(i,j) - rdt * tsc(3) * slurb_tile%tq_can(i,j) ) /         &
                     ( rdt  * tsc(2) )
   !
   !--       Compute the weighted RK3 tendency to be used in next time step.
            IF ( runge_l )  THEN
               IF ( rk3step == 1 )  THEN
                  slurb_tile%tq_can(i,j) = tq_new
               ELSEIF ( rk3step < 3 )  THEN
                  slurb_tile%tq_can(i,j) = -9.5625_field_r * tq_new + 5.3125_field_r * slurb_tile%tq_can(i,j)
               ENDIF
            ENDIF

            slurb_tile%q_can(i,j) = slurb_tile%q_can_p(i,j)
            slurb_tile%vpt_can(i,j) = slurb_tile%pt_can(i,j) * ( 1.0_field_r + 0.61_field_r * slurb_tile%q_can_p(i,j) )
            slurb_tile%qsws_can(i,j) = - f_qsws * ( slurb_tile%q1(i,j) - slurb_tile%q_can_p(i,j) ) / rlv

         ENDIF
   !
   !--    Compute the canyon horizontal wind speed, Eq. (9), Krayenhoff & Voogt (2007).
         slurb_tile%uv_abs_can(i,j) = slurb_tile%uv_abs_can_coef(i,j) * slurb_tile%uv_abs1(i,j)
   !
   !--    Calculate the canyon effective wind speed taking into account turbulent processes
   !--    Lemonsu et al. (2004) Eqs. (2-3).
   !
   !--    Free convection scale (wstar) at building roof height (for unstable cases)
   !--    In case of moist physics, use virtual temperature (buoyancy) flux in free-convection scale.
        !  write(*,*) "shf_can"
        !  write (*,*) i,j,slurb_tile%shf_can(i,j)
        !  write(*,*) "qsws_can"
        !  write (*,*) i,j,slurb_tile%qsws_can(i,j)
         IF ( moist_physics )  THEN
            vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j) + (rlv/cp) * slurb_tile%qsws_can(i,j)
         ELSE
            vtws =  (1/(rho_cp)) * slurb_tile%shf_can(i,j)
         ENDIF
   !
   !--    No scaling for stable cases:
        !  write(*,*) "VTWS"
        !  write (*,*) i,j,vtws
         vtws = MERGE( vtws, 0.0_field_r, vtws > 0.0_field_r )
         ws = ( g / slurb_tile%pt_can(i,j) * slurb_tile%z_mo_can(i,j) * vtws )**( 1.0_field_r / 3.0_field_r )
   !
   !--    Canyon effective wind speed taking account both mean and turbulent wind.
         slurb_tile%uv_eff_can(i,j) = SQRT( slurb_tile%uv_abs_can(i,j)**2 + ( slurb_tile%us_can(i,j) + ws )**2 )

      enddo
    enddo

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_canyon_model'
    !    CALL debug_message( debug_string, 'end' )
    ! ENDIF

 END SUBROUTINE slurb_canyon_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SLUrb's internal model to model urban surface - atmosphere coupling.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_urban_aggregation_model
    use modglobal, only : i1, j1, rlv, cp
    use modfields, only : u0, v0, rhof, qt0, thl0
    INTEGER ::  i       !< running index
    INTEGER ::  j       !< running index
    INTEGER ::  k_topo  !< k index of topography
    INTEGER ::  k_atm   !< k index of the first atmospheric level
    INTEGER ::  m       !< running index of surface tiles

    LOGICAL ::  runge_l  !< flag for timestep scheme to allow vectorization
    real :: rhocp_i, rholv_i

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_urban_aggregation_model'
    !    CALL debug_message( debug_string, 'start' )
    ! ENDIF

    ! runge_l = ( timestep_scheme(1:5) == 'runge' )
    runge_l = .true.

        do j=2,j1
      do i=2,i1
      !  k_topo = topo_top_ind(j,i,0)
      !  k_atm = topo_top_ind(j,i,0) + 1
       k_topo = 1
       k_atm = 1

        !
        !--    For shf and qsws, use direct aggregation.
       ! must be W m^-2
       slurb_tile%shf_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%shf_roof(i,j) +                                        &
                         ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%shf_can(i,j) + slurb_tile%shf_external(i,j)

       IF ( moist_physics )  THEN
          slurb_tile%qsws_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%qsws_roof(i,j) +                                   &
                             ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%qsws_can(i,j) + slurb_tile%qsws_external(i,j)
       ENDIF

!
!--    Calculate momentum flux for horizontal wind components.
       slurb_tile%usws_urb(i,j) = -u0(i,j,k_atm) / slurb_tile%ram_urb(i,j) * rho_air_zw(k_topo)
       slurb_tile%vsws_urb(i,j) = -v0(i,j,k_atm) / slurb_tile%ram_urb(i,j) * rho_air_zw(k_topo)

!
!--    Aggregate radiative fluxes. Note that this aggregation is done here rather than in the
!--    slurb_radiation_model on purpose to include the longwave term dependent on
!--    the surface temperature of given surface.
!
!--    Compute the net LW radiation flux at canyon top and for urban surface.
       slurb_tile%rad_lw_net_can(i,j) = slurb_tile%rad_lw_net_road(i,j) + slurb_tile%hw_can(i,j) *                         &
                                (   ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                                   &
                                    ( slurb_tile%rad_lw_net_wall_a(i,j) + slurb_tile%rad_lw_net_wall_b(i,j) )      &
                                  + slurb_tile%f_win(i,j) *                                                &
                                    ( slurb_tile%rad_lw_net_win_a(i,j)  + slurb_tile%rad_lw_net_win_b(i,j) )       &
                                )

       slurb_tile%rad_lw_net_urb(i,j) = slurb_tile%f_bld(i,j) * slurb_tile%rad_lw_net_roof(i,j) +                          &
                                ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * slurb_tile%rad_lw_net_can(i,j)
!
!--    Outgoing LW flux.
       slurb_tile%rad_lw_out_urb(i,j) = slurb_tile%rad_lw_in_urb(i,j) - slurb_tile%rad_lw_net_urb(i,j)
!
!--    Calculate urban aggregated surface temperatures.
       CALL calc_urban_aggregated_temperatures

        rhocp_i = 1. / (rhof(1) * cp)
        rholv_i = 1. / (rhof(1) * rlv)
        ! Calculate surface values
        slurb_tile%thlskin(i,j) = thl0(i,j,1) + (slurb_tile%shf_urb(i,j) * rhocp_i) * slurb_tile%ram_urb(i,j)
        slurb_tile%qtskin (i,j) = qt0(i,j,1) + (slurb_tile%qsws_urb(i,j) * rholv_i) * slurb_tile%ram_urb(i,j)

      enddo
    enddo

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_urban_aggregation_model'
    !    CALL debug_message( debug_string, 'end' )
    ! ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Update the aggregated urban surface temperatures. These are diagnostic outputs and are not
!> prognostic model variables. Four aggregated urban surface temperatures are computed:
!> 1) Effective surface temperature T_H derived from conservation of heat flux contributions
!> 2) Radiative surface temperature T_rad derived from the outgoing LW radiation
!> 3) Complete surface temperature T_C which is an area-weighted temperature of all facets
!> 4) Theoretical temperature at 2 m height extrapolated using stability-corrected log profile
!> For 1-3 formulations of Kanda et al. 2005, adapted for SLUrb configuration, are used.
!< Note that prognostic temperatures (suffix _p) are used. These are the ones that are output
!> for current time step, as the timelevel is swapped right after the prognostic equation calls.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_urban_aggregated_temperatures
    use modglobal, only : boltz, rlv, cp
    REAL(field_r) ::  c_h_roof    !< bulk heat transfer coefficient for roof (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_wall_a  !< bulk heat transfer coefficient for wall a (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_wall_b  !< bulk heat transfer coefficient for wall b (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_win_a   !< bulk heat transfer coefficient for window a (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_win_b   !< bulk heat transfer coefficient for window b (J kg^-1 K^-1)
    REAL(field_r) ::  c_h_road    !< bulk heat transfer coefficient for road (J kg^-1 K^-1)
    REAL(field_r) ::  ts          !< scaling temperature (K)
    REAL(field_r) ::  vtws        !< virtual potential temperature flux (buoyancy flux) (m K s^-1)
    real :: rho_cp !< (J m^-3 K^-1)
    real :: rho !< (kg m^-3)
    rho_cp = cp * rho_air_zw(k_topo)
    rho = rho_air_zw(k_topo)
!
!-- 1) Effective surface temperature T_H.
!-- First, compute the bulk heat transfer coefficients.
    IF ( calc_t_h )  THEN
       c_h_roof = ABS( slurb_tile%shf_roof(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                             &
                       ( slurb_tile%t_roof_p(nzt_roof,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_wall_a = ABS( slurb_tile%shf_wall_a(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                         &
                         ( slurb_tile%t_wall_a_p(nzt_wall,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_wall_b = ABS( slurb_tile%shf_wall_b(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                         &
                         ( slurb_tile%t_wall_b_p(nzt_wall,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_win_a = ABS( slurb_tile%shf_win_a(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                           &
                        ( slurb_tile%t_win_a_p(nzt_win,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_win_b = ABS( slurb_tile%shf_win_b(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                           &
                        ( slurb_tile%t_win_b_p(nzt_win,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       c_h_road = ABS( slurb_tile%shf_road(i,j) / ( rho * slurb_tile%uv_eff1(i,j) *                             &
                       ( slurb_tile%t_road_p(nzt_road,i,j) - slurb_tile%pt1(i,j) * exnf(k_atm) ) ) )

       slurb_tile%t_h_urb(i,j) = ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * (                                                    &
                                                ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                                ( c_h_wall_a * slurb_tile%t_wall_a_p(nzt_wall,i,j)         &
                                                + c_h_wall_b * slurb_tile%t_wall_b_p(nzt_wall,i,j) )       &
                                              + slurb_tile%f_win(i,j) *                                    &
                                                ( c_h_win_a * slurb_tile%t_win_a_p(nzt_win,i,j)            &
                                                + c_h_win_b * slurb_tile%t_win_b_p(nzt_win,i,j) )          &
                                              )                                                    &
                           + c_h_road * slurb_tile%t_road_p(nzt_road,i,j)                                  &
                           )                                                                       &
                         + slurb_tile%f_bld(i,j) * c_h_roof * slurb_tile%t_roof_p(nzt_roof,i,j)                    &
                         ) /                                                                       &
                         ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * (                                                    &
                                                ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                                ( c_h_wall_a + c_h_wall_b )                        &
                                              + slurb_tile%f_win(i,j) *                                    &
                                                ( c_h_win_a + c_h_win_b )                          &
                                              )                                                    &
                           + c_h_road                                                              &
                           )                                                                       &
                         + slurb_tile%f_bld(i,j) * c_h_roof + 1E-10_field_r                                     &
                         )
    ENDIF

!
!-- 2) Radiative surface temperature T_rad.
    slurb_tile%t_rad_urb(i,j) = SQRT( SQRT( slurb_tile%rad_lw_out_urb(i,j) / ( slurb_tile%emiss_urb(i,j) * boltz ) ) )

!
!-- 3) Complete surface temperature T_C, similarly to T_H but without the C_h weighting.
    IF ( calc_t_c )  THEN
       slurb_tile%t_c_urb(i,j) = ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                            &
                           ( slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                       &
                                     ( slurb_tile%t_wall_a_p(nzt_wall,i,j) + slurb_tile%t_wall_b_p(nzt_wall,i,j) ) &
                                     + slurb_tile%f_win(i,j) *                                             &
                                     ( slurb_tile%t_win_a_p(nzt_win,i,j)   + slurb_tile%t_win_b_p(nzt_win,i,j)   ) &
                                              )                                                    &
                           + slurb_tile%t_road_p(nzt_road,i,j)                                             &
                           )                                                                       &
                           + slurb_tile%f_bld(i,j) * slurb_tile%t_roof_p(nzt_roof,i,j)                             &
                         ) /                                                                       &
                         ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) * ( 2.0_field_r * slurb_tile%hw_can(i,j) + 1.0_field_r )       &
                           + slurb_tile%f_bld(i,j)                                                         &
                         )
    ENDIF

!
!-- 4) Theoretical 2 m temperature extrapolated using MOST.
    IF ( calc_t_2m )  THEN
       IF ( moist_physics )  THEN
          vtws =  (1/(rho*cp)) * slurb_tile%shf_can(i,j) + (rlv/cp) * slurb_tile%qsws_can(i,j)
       ELSE
          vtws =  (1/(rho*cp)) * slurb_tile%shf_can(i,j)
       ENDIF
       ts = -vtws * (1 / rhof(k_atm)) / slurb_tile%us_urb(i,j)

       slurb_tile%t_2m_urb(i,j) = ts / kappa *                                                             &
                          ( LOG( 2.0_field_r / ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) ) ) -                     &
                            psi_h( 2.0_field_r / slurb_tile%ol_urb(i,j) ) +                                     &
                            psi_h( ( slurb_tile%z_mo(i,j) + slurb_tile%h_bld(i,j) ) / slurb_tile%ol_urb(i,j) )             &
                          ) + slurb_tile%pt1(i,j) * exnf(k_atm)
    ENDIF





 END SUBROUTINE calc_urban_aggregated_temperatures

 END SUBROUTINE slurb_urban_aggregation_model

 !--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Shortwave and longwave radiation parametrisations of the model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE slurb_radiation_model
    use modglobal, only : i1,j1,xtime,rtimee,xday,xlat,xlon
    use modraddata, only : zenith_lon_lat
    INTEGER ::  day_of_year  !< day of year for the current day
    INTEGER ::  i            !< loop index
    INTEGER ::  j            !< loop index
    INTEGER ::  k_topo       !< k index of topography top
    INTEGER ::  k_atm        !< k index of the first atmospheric level
    INTEGER ::  m            !< running index of surface tiles

    REAL(field_r) ::  azimuth        !< solar azimuth angle
    REAL(field_r) ::  second_of_day  !< second of the current day
    REAL(field_r) ::  tan_zenith     !< tangent of the solar zenith angle
    REAL(field_r) ::  zenith         !< solar zenith angle
    real :: sun_dir_lon, sun_dir_lat, cos_zenith


    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_radiation_model'
    !    CALL debug_message( debug_string, 'start' )
    ! ENDIF

!
!-- Calculate solar angles if not already done by RTM.
    ! IF ( .NOT. radiation_interactions )  THEN
    !    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                  &
    !                        second_of_day = second_of_day )
        ! assume for now zenith is magically already calculated.... TODOSELF
    !    CALL calc_zenith( day_of_year, second_of_day )
    ! ENDIF

    call zenith_lon_lat(xtime*3600 + rtimee,xday,xlat,xlon,zenith, sun_dir_lon, sun_dir_lat)
    azimuth = ATAN2( sun_dir_lon, sun_dir_lat )
    cos_zenith = COS(zenith)


!
!-- Split the incoming SW radiation into direct and diffuse parts.
!-- Direct-diffuse SW split is quite weirdly done in the radiation mod if radiation
!-- interactions are enabled. However, we do need it here even without interactions.
   !  IF ( cos_zenith > 0.0_field_r )  CALL radiation_calc_diffusion_radiation

   do j=2,j1
      do i=2,i1

       k_topo = 1
       k_atm = 1

!
!--    Update SLUrb internal radiative fluxes based on the new surface temperatures
!--    Compute the internal longwave radiation interactions at every timestep.
       CALL calc_rad_lw

!
!--    Compute the SW radiation fluxesd.
!--    Do this only if the radiation model has updated SW fluxes at previous timestep,
!--    as otherwise the computation would just yield the same fluxes.
    !    IF ( radiation_called  .OR.  first_call )  CALL calc_rad_sw TODOSELF
       call calc_rad_sw
    ENDDO
   enddo

    ! IF ( debug_output_timestep )  THEN
    !    WRITE( debug_string, * ) 'slurb_radiation_model'
    !    CALL debug_message( debug_string, 'end' )
    ! ENDIF

!
!-- Private functions and subroutines of slurb_radiation_model.
    CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the LW radiative fluxes and their differentials for the time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rad_lw
   use modglobal, only : boltz
   use modraddata, only : lwd
    REAL(field_r) ::  t_rad_sky  !< Radiative temperature of the sky


!
!-- Compute the effective radiative temperature of the incoming LW radiation.
    slurb_tile%rad_lw_in_urb(i,j) = abs(lwd(i,j,1)) !TODOSELF
    ! lwd is positive in DALES
    t_rad_sky = SQRT( SQRT( slurb_tile%rad_lw_in_urb(i,j) / boltz ) )

!
!-- Computation of net LW fluxes based on Lemonsu et al. 2012 Eqs. (1-3) (+ windows).
!-- Note that these are NOT YET the final net longwave fluxes for the surfaces, as the term
!-- dependent on the surface's own surface temperature (coef=1) is omitted at this stage.
!-- This term is added after computing the prognostic equation for the surface temperature,
!-- as it is included in the prognostic equations in an linearized form.
    slurb_tile%rad_lw_net_roof(i,j) = slurb_tile%lw_roof_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j)


    slurb_tile%rad_lw_net_road(i,j) = slurb_tile%lw_road_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                     &
                              slurb_tile%lw_road_coef(3,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4 +              &
                              slurb_tile%lw_road_coef(3,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4 +              &
                              slurb_tile%lw_road_coef(4,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4 +                &
                              slurb_tile%lw_road_coef(4,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4

!
!-- The term dependent on t_wall_b is omitted at this stage, as for isotropic canyons the mean wall
!-- temperature is used, including both wall A and B interactions. Thus, the terms for both
!-- t_wall_a and t_wall_b have to be included in linearization. For anisotropic canyons there is
!-- no direct dependence, so it can be directly added (see below).
    slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%lw_wall_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                   &
                                slurb_tile%lw_wall_coef(4,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4 +              &
                                slurb_tile%lw_wall_coef(5,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4 +              &
                                slurb_tile%lw_wall_coef(6,i,j) * slurb_tile%t_road(nzt_road,i,j)**4

    IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
       slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%lw_win_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                  &
                                  slurb_tile%lw_win_coef(4,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4 +           &
                                  slurb_tile%lw_win_coef(5,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4 +           &
                                  slurb_tile%lw_win_coef(6,i,j) * slurb_tile%t_road(nzt_road,i,j)**4
    ENDIF

!
!-- Inverse for facade B, if anisotropic canyons are used. If not, copy.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
!
!--    In case of anisotropic canyons, t_wall_b doesn't have dependency on t_wall_a in the
!--    prognostic equation, and thus it's contribution to longwave balance can be directly added
!--    to the net longwave radiation before prognostic equations. Vice versa for t_wall_b.
       slurb_tile%rad_lw_net_wall_a(i,j) = slurb_tile%rad_lw_net_wall_a(i,j) +                                     &
                                   slurb_tile%lw_wall_coef(3,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4

!
!--    Note that for wall (and window) B the coefficients 4 and 5 are also swapped.
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%lw_wall_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +                &
                                   slurb_tile%lw_wall_coef(3,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4 +         &
                                   slurb_tile%lw_wall_coef(4,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4 +           &
                                   slurb_tile%lw_wall_coef(5,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4 +           &
                                   slurb_tile%lw_wall_coef(6,i,j) * slurb_tile%t_road(nzt_road,i,j)**4

       IF ( slurb_tile%f_win(i,j) > 0.0_field_r )  THEN
          slurb_tile%rad_lw_net_win_a(i,j) = slurb_tile%rad_lw_net_win_a(i,j) +                                    &
                                     slurb_tile%lw_win_coef(3,i,j) * slurb_tile%t_win_b(nzt_win,i,j)**4

          slurb_tile%rad_lw_net_win_b(i,j) = slurb_tile%lw_win_coef(2,i,j) * slurb_tile%rad_lw_in_urb(i,j) +               &
                                     slurb_tile%lw_win_coef(3,i,j) * slurb_tile%t_win_a(nzt_win,i,j)**4 +          &
                                     slurb_tile%lw_win_coef(4,i,j) * slurb_tile%t_wall_b(nzt_wall,i,j)**4 +        &
                                     slurb_tile%lw_win_coef(5,i,j) * slurb_tile%t_wall_a(nzt_wall,i,j)**4 +        &
                                     slurb_tile%lw_win_coef(6,i,j) * slurb_tile%t_road(nzt_road,i,j)**4
       ENDIF
    ELSE
       slurb_tile%rad_lw_net_wall_b(i,j) = slurb_tile%rad_lw_net_wall_a(i,j)
       slurb_tile%rad_lw_net_win_b(i,j)  = slurb_tile%rad_lw_net_win_a(i,j)
    ENDIF

 END SUBROUTINE calc_rad_lw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the SW radiative fluxes for the time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_rad_sw
   use modraddata, only : swdir, swdif

    REAL(field_r) ::  rad_sw_diff_road      !< incoming diffuse shortwave radiation on road
    REAL(field_r) ::  rad_sw_diff_wall_a    !< incoming diffuse shortwave radiation on wall A
    REAL(field_r) ::  rad_sw_diff_wall_b    !< incoming diffuse shortwave radiation on wall B
    REAL(field_r) ::  rad_sw_dir_road       !< incoming direct shortwave radiation on road
    REAL(field_r) ::  rad_sw_dir_wall_a     !< incoming direct shortwave radiation on wall A
    REAL(field_r) ::  rad_sw_dir_wall_b     !< incoming direct shortwave radiation on wall B
    REAL(field_r) ::  rad_sw_ref_nomin      !< nominator of the sum of reflections at infinity.
    REAL(field_r) ::  rad_sw_wall_modifier  !< modifier term for anisotropic walls
    REAL(field_r) ::  theta0                !< critical canyon orientation for road illumination
    REAL(field_r) ::  w_inf                 !< mean wall reflection at infinity



!
!-- Check if there is any shortwave radiation to take care of in the first place.
    IF ( .NOT. ( cos_zenith > 0.0_field_r ) )  THEN
       slurb_tile%rad_sw_in_urb(i,j)     = 0.0_field_r
       slurb_tile%rad_sw_net_urb(i,j)    = 0.0_field_r
       slurb_tile%rad_sw_net_roof(i,j)   = 0.0_field_r
       slurb_tile%rad_sw_net_road(i,j)   = 0.0_field_r
       slurb_tile%rad_sw_net_wall_a(i,j) = 0.0_field_r
       slurb_tile%rad_sw_net_wall_b(i,j) = 0.0_field_r
       slurb_tile%albedo_urb(i,j)        = 0.1_field_r
       RETURN
    ENDIF

    ! whatever radiation model we use, shortwave DOWN will always be positive, so we ensure that by taking absolute value.
    slurb_tile%rad_sw_in_urb(i,j) = abs(swdir(i,j,1)) + abs(swdif(i,j,1))

!
!-- Compute the net shortwave radiation for roofs, which is the simplest case.
    slurb_tile%rad_sw_net_roof(i,j) = ( 1.0_field_r - slurb_tile%albedo_roof(i,j) ) * slurb_tile%rad_sw_in_urb(i,j)

!
!-- Next, compute then et shortwave radiation within the street canyon. This is quite complex,
!-- including the effect of shading and within-canyon reflections. See Lemonsu et al. (2012)
!-- for reference.

!
!-- Calculate tangent of the zenith angle, with limiters and safety margins applied to prevent
!-- floating point overflows and division by zero. Shouldn't affect the physics too much.
    IF ( ABS( 0.5_field_r * pi - zenith ) < 1.0E-6_field_r )  THEN
       IF ( 0.5_field_r * pi - zenith >  0.0_field_r )  tan_zenith = TAN( 0.5_field_r * pi - 1.0E-6_field_r )
       IF ( 0.5_field_r * pi - zenith <= 0.0_field_r )  tan_zenith = TAN( 0.5_field_r * pi + 1.0E-6_field_r )
    ELSEIF ( ABS( zenith ) < 1.0E-6_field_r )  THEN
       tan_zenith = SIGN(1.0_field_r, zenith) * TAN( 1.0E-6_field_r )
    ELSE
       tan_zenith = TAN( zenith )
    ENDIF

!
!-- Direct SW radiation received by the walls (and windows), the road and vegetation.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
!
!--    Lemonsu et al. (2012) Eq. (A1)
!--    @note There is an error in this equation in the article. It should be that
!--    the direct radiation on road should decrease when difference between the sun azimuth
!--    angles increase, not vice versa.
       rad_sw_dir_road = abs(swdir(i,j,1)) * MAX( 0.0_field_r, 1.0_field_r - slurb_tile%hw_can(i,j) *               &
                         tan_zenith *  SIN( ABS( azimuth - slurb_tile%theta_can(i,j) ) ) )

!
!--    Lemonsu et al. (2012) Eqs. (A2-A4)
       rad_sw_dir_wall_a = ( abs(swdir(i,j,1)) - rad_sw_dir_road ) * 0.5_field_r / slurb_tile%hw_can(i,j)

       IF ( SIN( azimuth - slurb_tile%theta_can(i,j) ) > 0.0_field_r )  THEN
          rad_sw_dir_wall_a = 2.0_field_r * rad_sw_dir_wall_a
          rad_sw_dir_wall_b = 0.0_field_r
       ELSE
          rad_sw_dir_wall_b = 2.0_field_r * rad_sw_dir_wall_a
          rad_sw_dir_wall_a = 0.0_field_r
       ENDIF

    ELSE
!
!--    Revert to the anisotropic integrated solution by Masson (2000).
!
!--    Calculate the critical canyon orientation theta0 for anisotropic street canyons.
       theta0 = ASIN( MIN( 1.0_field_r / ( tan_zenith * slurb_tile%hw_can(i,j) ), 1.0_field_r ) )

!
!--    Masson (2000) Eqs. (13-15)
       rad_sw_dir_road = abs(swdir(i,j,1)) * ( 2.0_field_r * theta0 / pi -                             &
                         2.0_field_r * tan_zenith / pi * slurb_tile%hw_can(i,j) * ( 1.0_field_r - COS( theta0 ) ) )

       rad_sw_dir_wall_a = ( abs(swdir(i,j,1)) - rad_sw_dir_road ) * 0.5_field_r / slurb_tile%hw_can(i,j)

       rad_sw_dir_wall_b = rad_sw_dir_wall_a

   ENDIF

!
!-- Diffuse (from sky) solar radiation received by the surfaces.
    rad_sw_diff_road   = abs(swdif(i,j,1)) * slurb_tile%svf_road(i,j)
    rad_sw_diff_wall_a = abs(swdif(i,j,1)) * slurb_tile%svf_wall(i,j)
    rad_sw_diff_wall_b = rad_sw_diff_wall_a

!
!-- Canyon internal scattering based on both Masson (2000) Eqs. (16-20) and
!-- Lemonsu et al. (2012) Appendix A2. This has been modified to include windows: the weighted
!-- average reflection from walls and windows is taken into account by using weighted average
!-- albedo. The wall and window surfaces are assumed to be uniformly distributed.

!
!-- Nominator of the sum of reflections at infinity.
    rad_sw_ref_nomin = slurb_tile%albedo_wall_win(i,j) * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +        &
                       rad_sw_dir_wall_b + rad_sw_diff_wall_b ) / 2.0_field_r +                         &
                       slurb_tile%albedo_wall_win(i,j) * slurb_tile%svf_wall(i,j) * slurb_tile%albedo_road(i,j) *          &
                       rad_sw_dir_road

!
!-- Sum of refelctions at infinity.
    w_inf = rad_sw_ref_nomin / slurb_tile%sw_ref_denom(i,j)

!
!-- Total solar radiation absorbed after infinite reflections.
    slurb_tile%rad_sw_in_road(i,j) = rad_sw_dir_road + rad_sw_diff_road +                                  &
                             ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf
    slurb_tile%rad_sw_net_road(i,j) = ( 1.0_field_r - slurb_tile%albedo_road(i,j) ) * slurb_tile%rad_sw_in_road(i,j)

    slurb_tile%rad_sw_net_wall_a(i,j) = ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) *                                 &
                                ( 0.5_field_r * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a +              &
                                             rad_sw_dir_wall_b + rad_sw_diff_wall_b )              &
                                + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                         &
                                  ( rad_sw_dir_road + rad_sw_diff_road )                           &
                                + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                         &
                                  ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf                            &
                                + ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) * w_inf                   &
                                )

    slurb_tile%rad_sw_net_wall_b(i,j) = slurb_tile%rad_sw_net_wall_a(i,j)

    IF ( slurb_tile%f_win(i,j) /= 0.0_field_r  )  THEN
       slurb_tile%rad_sw_in_win_a(i,j) =   0.5_field_r * ( rad_sw_dir_wall_a + rad_sw_diff_wall_a               &
                                            + rad_sw_dir_wall_b + rad_sw_diff_wall_b )             &
                                 + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                        &
                                   ( rad_sw_dir_road + rad_sw_diff_road )                          &
                                 + slurb_tile%albedo_road(i,j) * slurb_tile%svf_wall(i,j) *                        &
                                   ( 1.0_field_r - slurb_tile%svf_road(i,j) ) * w_inf                           &
                                 + ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) * w_inf

       slurb_tile%rad_sw_net_win_a(i,j) = ( 1.0_field_r - slurb_tile%albedo_win(i,j) ) * slurb_tile%rad_sw_in_win_a(i,j)

       slurb_tile%rad_sw_in_win_b(i,j)  = slurb_tile%rad_sw_in_win_a(i,j)
       slurb_tile%rad_sw_net_win_b(i,j) = slurb_tile%rad_sw_net_win_a(i,j)
    ENDIF

!
!-- Modification of reflected solar radiation for anisotropic street canyons.
    IF ( slurb_tile%anisotropic_canyon(i,j) )  THEN
       rad_sw_wall_modifier = ( 1.0_field_r + slurb_tile%albedo_wall_win(i,j) *                                 &
                                ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) /                           &
                                ( 1.0_field_r + slurb_tile%albedo_wall_win(i,j) *                               &
                                  ( 1.0_field_r - 2.0_field_r * slurb_tile%svf_wall(i,j) ) )                         &
                              ) *                                                                  &
                              0.5_field_r * ( ( rad_sw_dir_wall_a + rad_sw_diff_wall_a )                &
                                       - ( rad_sw_dir_wall_b + rad_sw_diff_wall_b ) )

       slurb_tile%rad_sw_net_wall_a(i,j) = slurb_tile%rad_sw_net_wall_a(i,j) +                                     &
                                   ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) * rad_sw_wall_modifier

       slurb_tile%rad_sw_net_wall_b(i,j) = slurb_tile%rad_sw_net_wall_b(i,j) -                                     &
                                   ( 1.0_field_r - slurb_tile%albedo_wall(i,j) ) * rad_sw_wall_modifier

       IF ( slurb_tile%f_win(i,j) /= 0.0_field_r )  THEN
          slurb_tile%rad_sw_in_win_a(i,j)  = slurb_tile%rad_sw_in_win_a(i,j) + rad_sw_wall_modifier
          slurb_tile%rad_sw_net_win_a(i,j) = slurb_tile%rad_sw_in_win_a(i,j) * ( 1.0_field_r - slurb_tile%albedo_win(i,j) )
          slurb_tile%rad_sw_in_win_b(i,j)  = slurb_tile%rad_sw_in_win_b(i,j) - rad_sw_wall_modifier
          slurb_tile%rad_sw_net_win_b(i,j) = slurb_tile%rad_sw_in_win_b(i,j) * ( 1.0_field_r - slurb_tile%albedo_win(i,j) )
       ENDIF
    ENDIF

!
!-- The upward shortwave radiation is computed as residual of absorbed radiation per uniturban
!-- area. Aggregated effective albedo of urban surface is computed so that the raditaiton models end
!-- up with the same figure for outgoing shortwave radiation.
    slurb_tile%rad_sw_out_urb(i,j) = slurb_tile%rad_sw_in_urb(i,j) -                                               &
                             ( ( 1.0_field_r - slurb_tile%f_bld(i,j) ) *                                        &
                               ( slurb_tile%hw_can(i,j) * ( ( 1.0_field_r - slurb_tile%f_win(i,j) ) *                   &
                                       ( slurb_tile%rad_sw_net_wall_a(i,j) + slurb_tile%rad_sw_net_wall_b(i,j) )   &
                                       + slurb_tile%f_win(i,j) *                                           &
                                       ( slurb_tile%rad_sw_net_win_a(i,j)  + slurb_tile%rad_sw_net_win_b(i,j)  )   &
                                                  )                                                &
                               + slurb_tile%rad_sw_net_road(i,j)                                           &
                               )                                                                   &
                             + slurb_tile%f_bld(i,j) * slurb_tile%rad_sw_net_roof(i,j)                             &
                             )

!
!-- Compute the net SW flux for diagnostics and output.
    slurb_tile%rad_sw_net_urb(i,j) = slurb_tile%rad_sw_in_urb(i,j) - slurb_tile%rad_sw_out_urb(i,j)

!
!-- Save effective albedo for the radiation model.
    if (slurb_tile%rad_sw_in_urb(i,j) /= 0.0) then
     slurb_tile%albedo_urb(i,j) = slurb_tile%rad_sw_out_urb(i,j) / slurb_tile%rad_sw_in_urb(i,j)
    endif ! TODOSELF WHY IS SWD 0?

 END SUBROUTINE calc_rad_sw

 END SUBROUTINE slurb_radiation_model

end module modslurb