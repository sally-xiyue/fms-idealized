module mixed_layer_mod

!
! Implementation of mixed layer boundary condition
!

use                  fms_mod, only: set_domain, write_version_number, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING

use                  fms_mod, only: stdlog, check_nml_error, close_file,&
                                    open_namelist_file, stdout, file_exist, &
                                    read_data, write_data, open_file, &
                                    nullify_domain

use            constants_mod, only: HLV, PI, RHO_CP, CP_AIR, CP_OCEAN, OMEGA, RADIUS, TFREEZE ! TFREEZE is FW freezing pt in Kelvin, XZ 02/2018

use         diag_manager_mod, only: register_diag_field, send_data

use       time_manager_mod,   only: time_type

use           transforms_mod, only: get_deg_lat, grid_domain, get_lat_max, get_deg_lon ! get_deg_lon added by ZTAN 06/07/2011

use            vert_diff_mod, only: surf_diff_type

use          mpp_domains_mod, only: mpp_global_field

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.f90 $'

character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, replicate, mrdnl_gradient, gauss_smooth_field, mixed_layer, mixed_layer_end

!=================================================================================================================================

logical :: evaporation = .true.
logical :: ekman_layer = .false.
logical :: gauss_smooth = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
real    :: depth = 40.0
real    :: depth_land = 1.0
logical :: load_qflux = .false.
! Zonal asymmetric and Walker-circ ocean heating by ZTAN
real    :: zonasym_width = 16.0
real    :: zonasym_amp   = 0.0
real    :: zonasym_peaklat = 16.0
real    :: walker_e_lon    = 90.0
real    :: walker_w_lon    = 270.0
real    :: walker_n_lat    = 0.0
real    :: walker_s_lat    = 0.0
real    :: walker_n_amp    = 0.0
real    :: walker_s_amp    = 0.0
real    :: walker_width_ew = 30.0
real    :: walker_width_ns = 7.0
! Sea ice parameters by Ian Eisenman, XZ 02/2018
logical :: sea_ice_in_mixed_layer = .true.  ! include sea ice model
logical :: ice_as_albedo_only = .false.     ! no sea ice thermodynamics, but change albedo when ML temp = TFREEZE
logical :: sfc_melt_from_file = .false.     ! whether to compute ice surface melt based on input surface temperature file
character(len=64) :: sfc_melt_file = 'INPUT/t_sfc.nc' ! (optional) specified t_surf for ice surface melt
integer  :: num_input_times = 72   ! number of times during year in input t_surf file
! parameters for sea ice model, XZ 02/2018
real :: L_ice = 3e8 ! latent heat of fusion (J/m^3)
real, parameter  :: k_ice = 2 ! conductivity of ice (W/m/K)
real, parameter  :: ice_basal_flux_const = 120 ! linear coefficient for heat flux from ocean to ice base (W/m^2/K)
real, parameter  :: t_ice_base = tfreeze ! temperature at base of ice, taken to be freshwater freezing point

namelist/mixed_layer_nml/ evaporation, qflux_amp, depth, depth_land, qflux_width, load_qflux, ekman_layer, &
						zonasym_amp, zonasym_width, zonasym_peaklat, &
           				walker_e_lon, walker_w_lon, walker_n_lat, walker_s_lat, &
           				walker_n_amp, walker_s_amp, walker_width_ew, walker_width_ns, &
                  sea_ice_in_mixed_layer, ice_as_albedo_only, sfc_melt_from_file ! Namelist option added XZ 02/2018

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux
     id_flux_oceanq_sym,   &   ! oceanic Q flux  (zonal_symmetric part)
     id_flux_zonasym,      &   ! zonal asymmetric oceanic flux, ZTAN 06/2011
     id_flux_walker,       &   ! Walker-like oceanic flux by T.Merlis, ZTAN 05/03/2013
     id_flux_u,            &   ! surface stress
     id_flux_t                 ! sensible heat flux at surface
     ! Sea ice by Ian Eisenman, XZ 02/2018
     id_h_ice,             &   ! sea ice thickness
     id_a_ice,             &   ! sea ice fractional area
     id_t_ml,              &   ! mixed layer temperature

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux
     ocean_qflux_sym,       &   ! Q-flux (zonal_symmetric part)
     zon_asym_flux,         &   ! The extra flux due to zonal asymmetry, ZTAN 06/2011
     walker_flux,           &   ! The Walker flux by T.Merlis, ZTAN 05/03/2013
     rad_lat_2d,            &   ! latitude in radians
     rad_long_2d                ! longitude in radians
     ! Sea ice by Ian Eisenman, XZ 02/2018
     t_surf_for_melt            ! if sfc_melt_from_file, current time t_surf from input file (otherwise, =t_surf)

real, allocatable, dimension(:)   ::     &
     deg_lat, deg_long  ! long added 06/07/2011

real, allocatable, dimension(:,:)   ::                                        &
     gamma_t,               &   ! Used to calculate the implicit
     gamma_q,               &   ! correction to the diffusion in
     fn_t,                  &   ! the lowest layer
     fn_q,                  &   !
     en_t,                  &   !
     en_q,                  &   !
     alpha_t,               &   !
     alpha_q,               &   !
     alpha_lw,              &   !
     beta_t,                &   !
     beta_q,                &   !
     beta_lw,               &   !
     t_surf_dependence,     &   !
     corrected_flux,        &   !
     eff_heat_capacity,     &   ! Effective heat capacity
     delta_t_surf,          &   ! Increment in surface temperature
     depth_map                  ! 2d depth for vaiable heat capacity
     ! Sea ice by Ian Eisenman, XZ 02/2018
     dFdt_surf,             &   ! d(corrected_flux)/d(t_surf), for calculation of ice sfc temperature
     delta_t_ml,            &   ! Increment in mixed layer temperature
     delta_h_ice,           &   ! Increment in sea ice thickness
     delta_t_ice                ! Increment in ice (steady-state) surface temperature

real, allocatable, dimension(:,:)   ::                                        &
     rad_lat_replicate,         &
     rad_lat_diff_replicate,    &
     flux_u_replicate,          &
     flux_m_replicate,          &
     t_surf_replicate,          &
     rad_lat_global,            &
     rad_lat_diff_global,       &
     coriolis_global,           &
     flux_u_global,             &
     flux_m_global,             &
     t_surf_global,             &
     d_t_surf_d_lat_global,     &
     flux_h_integrand_global,   &
     flux_h_global,             &
     d_flux_h_cos_d_lat_global, &
     div_flux_h_global,         &
     smooth_div_flux_h_global

real, allocatable, dimension(:)   ::                                          &
     flux_h_1d,        &
     flux_m_1d,        &
     flux_u_1d,        &
     t_surf_1d

real, allocatable, dimension(:,:)  ::                                         &
     rad_lat_del,      &
     rad_lat_diff

real, allocatable, dimension(:,:,:) ::                                        &
     smooth_factor,    &
     smooth_field_X,   &
     ! Sea ice by Ian Eisenman, XZ 02/2018
     input_t_sfc            ! Sfc temp read from file for ice sfc melt


real inv_cp_air


!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, h_ice, a_ice, t_ml, q_surf, u_surf, v_surf, gust, &
           bucket_depth, ocean_mask, axes, Time)

type(time_type), intent(in)       :: Time
real, intent(inout), dimension(:,:) :: t_surf, q_surf, u_surf, v_surf, gust, h_ice, a_ice, t_ml
real, intent(inout), dimension(:,:,:) :: bucket_depth
real, intent(inout), dimension(:,:) :: ocean_mask
integer, intent(in), dimension(4) :: axes
integer, intent(in) :: is, ie, js, je, num_levels

integer :: i,j
real    :: rad_qwidth
integer:: ierr, io, unit
integer :: lat_max

if(module_is_initialized) return

call write_version_number(version, tagname)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=mixed_layer_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'mixed_layer_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=mixed_layer_nml)

call get_lat_max(lat_max)

allocate(rad_lat_2d              (is:ie, js:je))
allocate(ocean_qflux             (is:ie, js:je))
allocate(ocean_qflux_sym         (is:ie, js:je))
allocate(rad_long_2d             (is:ie, js:je))  ! Added by ZTAN 06/08/2011
allocate(zon_asym_flux           (is:ie, js:je))  ! Added by ZTAN 06/08/2011
allocate(walker_flux             (is:ie, js:je))  ! Added by ZTAN 05/06/2013

allocate(deg_lat                 (js:je))
allocate(deg_long                 (is:ie))! Added by ZTAN 06/08/2011

allocate(gamma_t                 (is:ie, js:je))
allocate(gamma_q                 (is:ie, js:je))
allocate(en_t                    (is:ie, js:je))
allocate(en_q                    (is:ie, js:je))
allocate(fn_t                    (is:ie, js:je))
allocate(fn_q                    (is:ie, js:je))
allocate(alpha_t                 (is:ie, js:je))
allocate(alpha_q                 (is:ie, js:je))
allocate(alpha_lw                (is:ie, js:je))
allocate(beta_t                  (is:ie, js:je))
allocate(beta_q                  (is:ie, js:je))
allocate(beta_lw                 (is:ie, js:je))
allocate(delta_t_surf            (is:ie, js:je))
allocate(eff_heat_capacity       (is:ie, js:je))
allocate(corrected_flux          (is:ie, js:je))
allocate(t_surf_dependence       (is:ie, js:je))
allocate(depth_map               (is:ie, js:je))
! Sea ice by Ian Eisenman, XZ 02/2018
allocate(t_surf_for_melt         (is:ie, js:je))
allocate(dFdt_surf               (is:ie, js:je))
allocate(delta_t_ml              (is:ie, js:je))
allocate(delta_h_ice             (is:ie, js:je))
allocate(delta_t_ice             (is:ie, js:je))
allocate(input_t_sfc             (ie-is+1, je-js+1, num_input_times+2))

if (ekman_layer) then
allocate(flux_u_1d                  (js:je))
allocate(flux_m_1d                  (js:je))
allocate(t_surf_1d                  (js:je))
allocate(rad_lat_replicate          (is:ie, js:je))
allocate(rad_lat_diff_replicate     (is:ie, js:je))
allocate(flux_u_replicate           (is:ie, js:je))
allocate(flux_m_replicate           (is:ie, js:je))
allocate(t_surf_replicate           (is:ie, js:je))

allocate(rad_lat_global             (is:ie, 1:lat_max))
allocate(rad_lat_diff_global        (is:ie, 1:lat_max))
allocate(coriolis_global            (is:ie, 1:lat_max))
allocate(flux_u_global              (is:ie, 1:lat_max))
allocate(flux_m_global              (is:ie, 1:lat_max))
allocate(t_surf_global              (is:ie, 1:lat_max))
allocate(d_t_surf_d_lat_global      (is:ie, 1:lat_max))
allocate(flux_h_integrand_global    (is:ie, 1:lat_max))
allocate(flux_h_global              (is:ie, 1:lat_max))
allocate(d_flux_h_cos_d_lat_global  (is:ie, 1:lat_max))
allocate(div_flux_h_global          (is:ie, 1:lat_max))
allocate(smooth_div_flux_h_global   (is:ie, 1:lat_max))

allocate(rad_lat_del                (is:ie, 1:lat_max))
allocate(rad_lat_diff               (is:ie, 1:lat_max))
allocate(smooth_factor              (is:ie, 1:lat_max, 1:lat_max))
allocate(smooth_field_X             (is:ie, 1:lat_max, 1:lat_max))
endif

! Sea ice by Ian Eisenman, XZ 02/2018
! read seasonally-varying input_sfc_melt(lat,lon,time)
if ( sfc_melt_from_file .and. file_exist(trim(sfc_melt_file)) ) then
  do j = 1, num_input_times
    call read_data( trim(sfc_melt_file), 't_surf', input_t_sfc(:,:,j+1), &
       & grid_domain, timelevel=j )
  end do
  ! wrap edges of seasonal cycle for interpolation
  input_t_sfc(:,:,1)=input_t_sfc(:,:,num_input_times+1)
  input_t_sfc(:,:,num_input_times+2)=input_t_sfc(:,:,2)
endif


!
!see if restart file exists for the surface temperature
!
if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
! added
   call read_data(trim('INPUT/mixed_layer.res'), 'q_surf',   q_surf, grid_domain)
   call read_data(trim('INPUT/mixed_layer.res'), 'u_surf',   u_surf, grid_domain)
   call read_data(trim('INPUT/mixed_layer.res'), 'v_surf',   v_surf, grid_domain)
   call read_data(trim('INPUT/mixed_layer.res'), 'gust',     gust,   grid_domain)
! end of addition
   call read_data(trim('INPUT/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)
else
  call error_mesg('mixed_layer','mixed_layer restart file not found', WARNING)
endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')
id_flux_oceanq_sym = register_diag_field(mod_name, 'flux_oceanq_sym',        &
                                 axes(1:2), Time, 'oceanic Q-flux, symmetric part','watts/m2')
id_flux_zonasym = register_diag_field(mod_name, 'flux_zonasym',        &
                                 axes(1:2), Time, 'Equivalent flux of zonal asymmetry','watts/m2')
id_flux_walker  = register_diag_field(mod_name, 'flux_walker',        &
                                 axes(1:2), Time, 'Equivalent flux of Walker ocean flux','watts/m2')
! sea ice by Ian Eisenman, XZ 02/2018
id_h_ice = register_diag_field(mod_name, 'h_ice',        &
                                axes(1:2), Time, 'sea ice thickness','m')
id_a_ice = register_diag_field(mod_name, 'a_ice',        &
                                axes(1:2), Time, 'sea ice area','fraction of grid box')
id_t_ml = register_diag_field(mod_name, 't_ml',        &
                                axes(1:2), Time, 'mixed layer tempeature','K')

! set up depth_map for spatially varying heat capacity  ! Added by LJJ
where (ocean_mask > 0.0)
  depth_map(:,:) = depth
elsewhere
  depth_map(:,:) = depth_land
endwhere

! End LJJ addition

! latitude will be needed for oceanic q flux
call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo

call get_deg_lon(deg_long)
do i=is,ie
  rad_long_2d(i,:) = deg_long(i)*PI/180.
enddo


! calculate ocean Q flux
rad_qwidth = qflux_width*PI/180.
ocean_qflux_sym = qflux_amp*(1-2.*rad_lat_2d**2/rad_qwidth**2) * &
        exp(- ((rad_lat_2d)**2/(rad_qwidth)**2))
zon_asym_flux =  zonasym_amp* sin(rad_long_2d)*&
        exp(- ((abs(rad_lat_2d) - zonasym_peaklat*PI/180. )**2/(zonasym_width*PI/180.)**2))
walker_flux = walker_n_amp*exp( -( mod(rad_long_2d+PI - walker_e_lon*PI/180, 2*PI)- PI)**2/(walker_width_ew*PI/180.)**2  - &
                (rad_lat_2d-walker_n_lat*PI/180.)**2/(walker_width_ns*PI/180.)**2) - &
              walker_n_amp*exp( -( mod(rad_long_2d+PI - walker_w_lon*PI/180, 2*PI)- PI)**2/(walker_width_ew*PI/180.)**2  - &
                (rad_lat_2d-walker_n_lat*PI/180.)**2/(walker_width_ns*PI/180.)**2) + &
              walker_s_amp*exp( -( mod(rad_long_2d+PI - walker_e_lon*PI/180, 2*PI)- PI)**2/(walker_width_ew*PI/180.)**2  - &
                (rad_lat_2d-walker_s_lat*PI/180.)**2/(walker_width_ns*PI/180.)**2) - &
              walker_s_amp*exp( -( mod(rad_long_2d+PI - walker_w_lon*PI/180, 2*PI)- PI)**2/(walker_width_ew*PI/180.)**2  - &
                (rad_lat_2d-walker_s_lat*PI/180.)**2/(walker_width_ns*PI/180.)**2)
ocean_qflux = ocean_qflux_sym + walker_flux !+ zon_asym_flux

! zero Q flux over land         ! Added by LJJ
where (ocean_mask > 0.0)
     ocean_qflux(:,:) = ocean_qflux(:,:)
elsewhere
     ocean_qflux(:,:) = 0.0
endwhere

! End LJJ addition

! load Q flux
if (load_qflux) then
  call read_data('INPUT/ocean_qflux.nc', 'ocean_qflux',  ocean_qflux)
endif

inv_cp_air = 1.0 / CP_AIR

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine replicate(ni, nf, field_2D, field_1D)

   integer :: i

   integer, intent(in) ::                                               &
       ni,                                                              &
       nf

    real, dimension(:), intent(in) ::                                   &
       field_1D             ! 1d input field

    real, dimension(:,:), intent(inout) ::                              &
       field_2D             ! 2d output field

    field_2D = 0
    do i = 1, (nf-ni+1)
       field_2D(i,:) = field_1D(:)
    enddo

end subroutine replicate

subroutine mrdnl_gradient(field_X, rad_lat, d_field_X_d_lat)

    integer :: l, k

    real, dimension(:,:), intent(in) ::                                 &
         field_X             ! 2d input field

    real, dimension(:,:), intent(in) ::                                 &
         rad_lat             ! 2d latitude in radian

    real, dimension(:,:), intent(out) ::                                &
         d_field_X_d_lat     ! meridional gradient of field_X


    rad_lat_del  = 0
    rad_lat_diff = 0
    do l=2, (size(rad_lat,2)-1)
       rad_lat_del(:,l)   = (rad_lat(:,l+1)-rad_lat(:,l))/(rad_lat(:,l)-rad_lat(:,l-1))
       rad_lat_diff(:,l)  = (rad_lat(:,l+1)-rad_lat(:,l))
    enddo
    rad_lat_del(:,1)                = (rad_lat(:,2)-rad_lat(:,1))/(rad_lat(:,1)-rad_lat(:,size(rad_lat,2)))
    rad_lat_del(:,size(rad_lat,2))  = (rad_lat(:,1)-rad_lat(:,size(rad_lat,2)))/(rad_lat(:,size(rad_lat,2))-rad_lat(:,size(rad_lat,2)-1))
    rad_lat_diff(:,1)               = rad_lat(:,2) - rad_lat(:,1)
    rad_lat_diff(:,size(rad_lat,2)) = rad_lat(:,1) - rad_lat(:,size(rad_lat,2))


    d_field_X_d_lat = 0
    do k=2, (size(field_X,2)-1)
       d_field_X_d_lat(:,k) = (field_X(:,k+1) - field_X(:,k)*(1-rad_lat_del(:,k)**2) - field_X(:,k-1)*(rad_lat_del(:,k)**2))/((1+rad_lat_del(:,k))*rad_lat_diff(:,k))
    enddo
end subroutine mrdnl_gradient

subroutine gauss_smooth_field(field_X, rad_lat, stand_dev, gauss_smooth_field_X)

    integer :: i, j

    real, intent(in) ::                                                    &
          stand_dev

    real, dimension(:,:), intent(in) ::                                    &
         field_X             ! 2d input field

    real, dimension(:,:), intent(in) ::                                    &
         rad_lat             ! 2d latitude in radian

    real ::                                                                &
         coeff_max

    real, dimension(:,:), intent(out) ::                                   &
         gauss_smooth_field_X

   coeff_max             = 0
   smooth_factor         = 0
   smooth_field_X        = 0
   gauss_smooth_field_X  = 0

   coeff_max = 1/(stand_dev*sqrt(2*pi))
   do j = 2, (size(field_X,2)-1)
      do i = 2, (size(field_X,2)-1)
         smooth_factor(:,i,j) = coeff_max*exp(-((rad_lat(:,i)-rad_lat(:,j))**2)/(2*(stand_dev**2))) &
                   *(rad_lat(:,i+1)-rad_lat(:,i-1))/2
         smooth_field_X(:,i,j) = field_X(:,i)*smooth_factor(:,i,j)
      enddo
      gauss_smooth_field_X(:,j) = sum(smooth_field_X(:,1:size(field_X,2),j),2)/sum(smooth_factor(:,1:size(field_X,2),j),2)
   enddo

end subroutine gauss_smooth_field


!=================================================================================================================================


!=================================================================================================================================

subroutine mixed_layer (                                               &
     is,                                                               &
     ie,                                                               &
     js,                                                               &
     je,                                                               &
     Time,                                                             &
     t_surf,                                                           &
     h_ice,                                                            &
     a_ice,                                                            &
     t_ml,                                                             &
     flux_t,                                                           &
     flux_q,                                                           &
     flux_r,                                                           &
     flux_u,                                                           &
     dt,                                                               &
     net_surf_sw_down,                                                 &
     surf_lw_down,                                                     &
     Tri_surf,                                                         &
     dhdt_surf,                                                        &
     dedt_surf,                                                        &
     dedq_surf,                                                        &
     drdt_surf,                                                        &
     dhdt_atm,                                                         &
     dedq_atm)




! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)       :: Time
real, intent(in),  dimension(:,:) :: &
     net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:) :: &
     flux_t,    flux_q,     flux_r,    flux_u
real, intent(inout), dimension(:,:) :: t_surf, h_ice, a_ice, t_ml
real, intent(in), dimension(:,:) :: &
   dhdt_surf, dedt_surf, dedq_surf, &
   drdt_surf, dhdt_atm, dedq_atm
real, intent(in) :: dt
type(surf_diff_type), intent(inout) :: Tri_surf
integer, intent(in) :: is, ie, js, je
integer :: lat_max, nh_lath, sh_lath, j
real :: rad_lat_std

! for sfc_melt_from_file, sea ice by Ian Eisenman, XZ 2/10/2018
integer :: i, j, n, seconds, days
real    :: days_in_year    = 360 ! how many model days in solar year
real    :: day = 0.0

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps

! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature

gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_q + dedq_atm))

fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_q + Tri_surf%dtmass * flux_q)

en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf

!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q + dedq_atm * fn_q
alpha_lw = flux_r

beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

!##########################################################################################
!                                START OCEAN QFLUX
!##########################################################################################
if (ekman_layer) then

! INITIALIZE
ocean_qflux = 0
flux_u_1d = 0
flux_m_1d = 0
t_surf_1d = 0
flux_u_replicate = 0
flux_m_replicate = 0
t_surf_replicate = 0
rad_lat_global            = 0
rad_lat_diff_global       = 0
flux_u_global             = 0
flux_m_global             = 0
t_surf_global             = 0
d_t_surf_d_lat_global     = 0
flux_h_integrand_global   = 0
flux_h_global             = 0
d_flux_h_cos_d_lat_global = 0
div_flux_h_global         = 0
smooth_div_flux_h_global  = 0

call get_lat_max(lat_max)

call mpp_global_field(grid_domain, rad_lat_2d, rad_lat_global)
do j=2,(lat_max-1)
   rad_lat_diff_global(:,j)  = (rad_lat_global(:,j+1)-rad_lat_global(:,j-1))/2
enddo
coriolis_global = 2.0*omega*sin(rad_lat_global)

! A. Fields are zonally averaged, zonally replicated, and carried globally
flux_u_1d  = sum(flux_u,1)/size(flux_u,1)
flux_m_1d  = flux_u_1d / coriolis_global(1,:)
t_surf_1d  = sum(t_surf,1)/size(t_surf,1)

call replicate(is, ie, flux_u_replicate, flux_u_1d)
call mpp_global_field(grid_domain, flux_u_replicate, flux_u_global)

call replicate(is, ie, flux_m_replicate, flux_m_1d)
call mpp_global_field(grid_domain, flux_m_replicate, flux_m_global)

call replicate(is, ie, t_surf_replicate, t_surf_1d)
call mpp_global_field(grid_domain, t_surf_replicate, t_surf_global)

call mrdnl_gradient(t_surf_global, rad_lat_global, d_t_surf_d_lat_global)

! B. Determined latitude where surface wind changes sign in:
! a. The Northern Hemisphere
nh_lath = lat_max/2+1
do j=(lat_max-1),(lat_max/2+1),-1
   if ((sign(1.0, flux_u_global(1,j+1))  .ne.      &
         sign(1.0, flux_u_global(1,j)))   .and.     &
             (flux_u_global(1,j)>0)         .and.     &
       (j .ge. ((lat_max/2+1) + lat_max/12)))  then
      nh_lath = j+1
   endif
enddo
!b. The Southern Hemisphere
sh_lath = lat_max/2
do j = 2,(lat_max/2),1
   if ((sign(1.0, flux_u_global(1,j-1)) .ne.      &
        sign(1.0, flux_u_global(1,j)))   .and.     &
            (flux_u_global(1,j)>0)         .and.     &
      (j .le. (lat_max/2 - lat_max/12)))      then
        sh_lath = j-1
   endif
enddo

! C. Compute integrand of enthalpy flux in the tropics
do j = (sh_lath+1), (nh_lath-1)
   flux_h_integrand_global(:,j) = cp_ocean*flux_u_global(:,j)/coriolis_global(:,j)*d_t_surf_d_lat_global(:,j)
enddo

! D. Compute enthalpy flux in the tropics
do j = (sh_lath), (lat_max/2)
    flux_h_global(:,j) =  sum(flux_h_integrand_global(:,sh_lath:j)*rad_lat_diff_global(:,sh_lath:j),2)
enddo
do j = (lat_max/2+1), (nh_lath)
    flux_h_global(:,j) = - sum(flux_h_integrand_global(:,j:nh_lath)*rad_lat_diff_global(:,j:nh_lath),2)
enddo

! E. Compute enthalpy flux divergence
call mrdnl_gradient(flux_h_global*cos(rad_lat_global), rad_lat_global, d_flux_h_cos_d_lat_global)
do j=(sh_lath), (nh_lath)
   div_flux_h_global(:,j) = 1/(radius*cos(rad_lat_global(:,j)))*d_flux_h_cos_d_lat_global(:,j)
enddo

! F. Gaussian smooting function
rad_lat_std = 7*PI/180
if (gauss_smooth) then
   call gauss_smooth_field(div_flux_h_global, rad_lat_global, rad_lat_std, smooth_div_flux_h_global)
else
   smooth_div_flux_h_global = div_flux_h_global
endif

! G. Set flux_oceanq
do j= js, je
   ocean_qflux(:,j) = smooth_div_flux_h_global(:,j)
enddo

endif
!write(*,*) 'smooth_div_flux_h', smooth_div_flux_h_global(1,:)
!write(*,*) 'sh_lath', sh_lath, 'nh_lath', nh_lath
!##########################################################################################
!                                 END OCEAN QFLUX
!##########################################################################################


!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw ! + ocean_qflux, XZ 02/2018
t_surf_dependence = beta_t * CP_AIR + beta_lw


if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

! for calculation of ice sfc temperature, sea ice by Ian Eisenman, XZ 02/2018
dFdt_surf = dhdt_surf + drdt_surf + HLV * (dedt_surf) ! d(corrected_flux)/d(T_surf)

! Sea ice by Ian Eisenman, XZ 02/2018
! === (3) surface temperature increment is calculated with
! explicit forward time step using flux corrected for implicit atm;
! next, (4) surface state is updated  ===

! ======================================================================
!
!          Ocean mixed layer and sea ice model equations
!
! ======================================================================
!
! Mixed layer temperature is evolved is t_ml, and sea ice thickness is
! h_ice. Atmosphere cares about surface temparature (t_surf). Where
! h_ice=0, t_surf=t_ml, but where h_ice>0, t_surf=t_ice which is the
! ice surface temperature (assumed to be in steady-state: "zero layer
! model"). So h_ice, t_ml, and t_surf are all saved at each time step
! (t_ice does not need to be saved).
!
! This model is equivalent to the Semtner (1976) "zero-layer model",
! with several differences (no snow, ice latent heat is same at base
! as surface, ice surface temperature calculated by including sensible
! and latent heat derivatives in dF/dT_surf, inclusion of frazil
! growth).

if ( sea_ice_in_mixed_layer .and. .not. ice_as_albedo_only ) then ! === use sea ice model ===

  ! calculate values of delta_h_ice, delta_t_ml, and delta_t_ice

  where ( h_ice .le. 0 ) ! = ice-free = [ should be equivalent to use .eq. instead of .le. ]
    delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
    delta_h_ice = 0
  elsewhere ! = ice-covered = ( h>0 )
    delta_t_ml = - ( ice_basal_flux_const * ( t_ml - t_ice_base ) + ocean_qflux ) * dt/(depth*RHO_CP)
    delta_h_ice = ( corrected_flux - ice_basal_flux_const * ( t_ml - t_ice_base ) ) * dt/L_ice
  endwhere
  where ( t_ml + delta_t_ml .lt. TFREEZE ) ! = frazil growth =
    delta_h_ice = delta_h_ice - ( t_ml + delta_t_ml - TFREEZE ) * (depth*RHO_CP)/L_ice
    delta_t_ml = TFREEZE - t_ml
  endwhere
  where ( ( h_ice .gt. 0 ) .and. ( h_ice + delta_h_ice .le. 0 ) ) ! = complete ablation =
    delta_t_ml = delta_t_ml - ( h_ice + delta_h_ice ) * L_ice/(depth*RHO_CP)
    delta_h_ice = - h_ice
  endwhere
  ! = update surface temperature =
  if ( sfc_melt_from_file ) then ! sfc melt from seasonally varying sfc temp specified from file
    !
    ! linearly interpolate from input file time to model time at each location
    ! interp1( ( (1:num_input_times)-0.5 )*days_in_year/num_input_times, albedo(x,y,:), day )
    ! = interp1( 1:num_input_times, albedo(x,y,:), day*num_input_times/days_in_year+0.5 )
    ! find model time
    call get_time(Time,seconds,days)
    day = days + seconds/86400
    ! make sure day is between 0 and days_in_year=360
    do while (day .lt. 0)
      day = day + days_in_year
    end do
    do while (day .ge. days_in_year)
      day = day - days_in_year
    end do
    ! find index of nearest input time below
    n=floor(day*num_input_times/days_in_year+1.5)
    do i = 1, size(input_t_sfc,1)
      do j = 1, size(input_t_sfc,2)
        t_surf_for_melt(i,j)=input_t_sfc(i,j,n) + (input_t_sfc(i,j,n+1)-input_t_sfc(i,j,n)) * &
           ( (day*num_input_times/days_in_year+1.5)-n )
      enddo
    enddo
    ! compute surface melt
    where ( h_ice + delta_h_ice .gt. 0 ) ! surface is ice-covered
      ! calculate increment in steady-state ice surface temperature
      delta_t_ice = ( - corrected_flux + k_ice / (h_ice + delta_h_ice) * ( t_ice_base - t_surf ) ) &
                    / ( k_ice / (h_ice + delta_h_ice) + dFdt_surf )
      ! in grid boxes with ice, wherever input t_surf=t_fr, make t_surf=t_fr;
      ! otherwise, let t_surf be whatever it wants (even t_surf>t_fr)
      where ( t_surf_for_melt .ge. TFREEZE ) ! surface ablation
        delta_t_ice = TFREEZE - t_surf
      endwhere
      ! surface is ice-covered, so update t_surf as ice surface temperature
      t_surf = t_surf + delta_t_ice
    elsewhere ! ice-free, so update t_surf as mixed layer temperature
      t_surf = t_ml + delta_t_ml
    endwhere
    !
  else ! no file specifying sfc melt (default)
    ! compute surface melt
    where ( h_ice + delta_h_ice .gt. 0 ) ! surface is ice-covered
      ! calculate increment in steady-state ice surface temperature
      delta_t_ice = ( - corrected_flux + k_ice / (h_ice + delta_h_ice) * ( t_ice_base - t_surf ) ) &
                    / ( k_ice / (h_ice + delta_h_ice) + dFdt_surf )
      where ( t_surf + delta_t_ice .gt. TFREEZE ) ! surface ablation
        delta_t_ice = TFREEZE - t_surf
      endwhere
      ! surface is ice-covered, so update t_surf as ice surface temperature
      t_surf = t_surf + delta_t_ice
    elsewhere ! ice-free, so update t_surf as mixed layer temperature
      t_surf = t_ml + delta_t_ml
    endwhere
  endif

else ! === do not use sea ice model: just evolve mixed layer with explicit step ===
  ! where ( land_mask .eq. 0 ) ! ocean
  !   delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
  ! elsewhere ! land
  !   delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth_land*RHO_CP)
  ! endwhere
  delta_t_ml = - ( corrected_flux + ocean_qflux) * dt/(depth*RHO_CP)
  delta_h_ice = 0 ! do not evolve sea ice (no ice model)
  ! t_surf and t_ml are equal (they differ only when mixed layer is ice-covered)
  t_surf = t_ml + delta_t_ml
endif

! = update state =
t_ml = t_ml + delta_t_ml
h_ice = h_ice + delta_h_ice

where ( h_ice .gt. 0 )
  a_ice = 1
elsewhere
  a_ice = 0
endwhere

if ( sea_ice_in_mixed_layer .and. ice_as_albedo_only ) then ! === model with only albedo changing when T < Tfr ===
  ! a_ice is passed to radiation module for albedo. just set a_ice=1 where T<Tfr.
  where ( t_ml .le. tfreeze)
    a_ice = 1
  elsewhere
    a_ice = 0
  endwhere
endif

! ======================================================================

! === (5) save increments in T and Q for lowest atm layer ===

Tri_surf%delta_t = fn_t
Tri_surf%delta_q = fn_q

! Original implicit stepping is commented out, XZ 02/2018
!
! ! Now update the mixed layer surface temperature using an implicit step
! !
! !eff_heat_capacity = depth * RHO_CP + t_surf_dependence * dt
! eff_heat_capacity = depth_map * RHO_CP + t_surf_dependence * dt
!
! if (any(eff_heat_capacity .eq. 0.0))  then
!   write(*,*) 'mixed_layer: error', eff_heat_capacity
!   call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
! end if
!
! delta_t_surf = - corrected_flux  * dt / eff_heat_capacity
!
! t_surf = t_surf + delta_t_surf

!
! Finally calculate the increments for the lowest atmospheric layer
!
! Tri_surf%delta_t = fn_t + en_t * delta_t_surf
! Tri_surf%delta_q = fn_q + en_q * delta_t_surf


!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
!
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q, Time)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time)
if(id_flux_oceanq_sym > 0)   used = send_data(id_flux_oceanq_sym, ocean_qflux_sym, Time)
if(id_flux_zonasym > 0)   used = send_data(id_flux_zonasym, zon_asym_flux, Time)
if(id_flux_walker > 0)    used = send_data(id_flux_walker,  walker_flux, Time)
! Sea ice by Ian Eisenman, XZ 2/10/2018
if(id_h_ice > 0) used = send_data(id_h_ice, h_ice, Time)
if(id_a_ice > 0) used = send_data(id_a_ice, a_ice, Time)
if(id_t_ml > 0) used = send_data(id_t_ml, t_ml, Time)

end subroutine mixed_layer

!=================================================================================================================================

subroutine mixed_layer_end(t_surf, h_ice, a_ice, t_ml, q_surf, u_surf, v_surf, gust, bucket_depth)

real, intent(inout), dimension(:,:) :: t_surf, q_surf, u_surf, v_surf, gust, h_ice, a_ice, t_ml ! Added by XZ 02/2018
real, intent(inout), dimension(:,:,:) :: bucket_depth
integer:: unit

if (ekman_layer) then

deallocate(flux_u_replicate)
deallocate(flux_m_replicate)
deallocate(t_surf_replicate)

deallocate(rad_lat_global)
deallocate(rad_lat_diff_global)
deallocate(coriolis_global)
deallocate(flux_u_global)
deallocate(flux_m_global)
deallocate(t_surf_global)

deallocate(d_t_surf_d_lat_global)

deallocate(flux_h_integrand_global)
deallocate(flux_h_global)
deallocate(d_flux_h_cos_d_lat_global)
deallocate(div_flux_h_global)
deallocate(smooth_div_flux_h_global)

endif


if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
! added ZTAN
call write_data(trim('RESTART/mixed_layer.res'), 'q_surf',   q_surf, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'u_surf',   u_surf, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'v_surf',   v_surf, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'gust',     gust,   grid_domain)
! end of addition
call write_data(trim('RESTART/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)
! Sea ice by Ian Eisenman, XZ 2/11/2018
call write_data(trim('RESTART/mixed_layer.res'), 'h_ice', h_ice, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'a_ice', a_ice, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 't_ml', t_ml, grid_domain)

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
