! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! Module modified for use with MSHM
! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use        types_mod, only : r8, i8, MISSING_R8,    &
                             obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date,&
                             set_calendar_type,increment_time 

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, &
                             get_location, &
                             VERTISSURFACE, VERTISHEIGHT

use    utilities_mod, only : register_module, error_handler, &
                             open_file, close_file, &
                             get_unit, E_ERR, E_MSG,  to_upper, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, file_exist, check_namelist_read

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 nc_add_attribute_to_variable, nc_define_dimension,   &
                                 nc_define_real_variable, nc_define_integer_variable, &
                                 nc_put_variable, nc_open_file_readonly,              &
                                 nc_create_file, nc_close_file, &
                                 nc_get_variable, nc_define_real_scalar,  &
                                 nc_define_integer_scalar, nc_get_dimension_size 

use state_structure_mod, only : add_domain, get_domain_size, get_kind_index, &
                                get_model_variable_indices,get_varid_from_kind,&
                                get_dart_vector_index

use     obs_kind_mod,      only : get_index_for_quantity

use ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies,    &   !CPS , read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step
use netcdf

implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time,       &
          get_model_time

!CPS specific to mshm but not required by DART
public :: initialize_mshm_geo,      &
          write_mshm_restart,       &
          write_dart_netcdf         

character(len=256), parameter :: source   = "model_mod.f90"
logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step 

character(len=256)               :: string1, string2, string3
integer                          :: model_size
type(time_type)                  :: model_time

integer(kind=4)                  :: &
            ny,                     &    ! Latitude dimension, 
            nx,                     &    ! Longitude dimension,
            nz,                     &    ! Maximum vertical levels
            npts                         ! Total number of grid points
real(r8)                         :: &
            dlat,                   &    ! Grid resolution (degrees)
            dlon                         ! Grid resolution (degrees)

real(r8), allocatable            :: &
           lat(:),lat2d(:,:),       &    ! Geo-coordinates
           lon(:),lon2d(:,:),       &    ! Geo-coordinates
           levels(:)

real(r8), parameter              :: missingValue = -9999._r8

! Codes for interpreting the columns of the state_variables
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not
integer, parameter :: MAX_STATE_VARIABLES = 40 
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5

! Namelist option for setup at runtime
character(len=256)               :: mshm_geo_file = 'latlon_mshm.ascii'
character(len=256)               :: mshm_restart_nc = 'dart_prior.nc'
character(len=32)                :: calendar = 'Gregorian'
integer                          :: time_step_days    = 0
integer                          :: mshm_sv_fields    = 10    !+3 added 
integer                          :: time_step_seconds = 3600
character(len=obstypelength)     :: mshm_variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' '

namelist /model_nml/           &
     mshm_geo_file,            &
     mshm_sv_fields,           &
     time_step_days,           &
     calendar,                 &
     mshm_restart_nc,          &     !File created by mshm to dart
     time_step_seconds,        &
     mshm_variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io
integer  :: num_vars
character(len=obstypelength) :: var_names(MAX_STATE_VARIABLES)
real(r8) :: var_ranges(MAX_STATE_VARIABLES,2)
logical  :: var_update(MAX_STATE_VARIABLES)
integer  :: var_qtys(  MAX_STATE_VARIABLES)

!
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source)

! Read the DART namelist
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Set the calendar and get time
call set_calendar_type( calendar ) 
model_time = read_model_time(mshm_restart_nc)

! Get geo location
 call  extract_nc_geo(mshm_restart_nc)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(time_step_seconds, &
                                  time_step_days)

! Define which variables to use from the model states
 call check_mvar_assim(mshm_variables, num_vars, &
                   var_names, var_qtys, var_ranges, var_update)

 dom_id = add_domain(mshm_restart_nc, num_vars, var_names, &
                kind_list=var_qtys, &
                clamp_vals=var_ranges(1:num_vars,:), &
                update_list=var_update)

! Set model size
 model_size = get_model_size()


end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

! We want to get only size of choosen state vectors
get_model_size = get_domain_size(dom_id)   !model_size

end function get_model_size


!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.
!
! For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)
!
character(len=256), parameter :: routine   = "model_interpolate"
real(r8) :: loc(3), loc_lon, loc_lat, loc_depth
integer  :: varid
integer  :: geo_inds(4)
real(r8) :: geo_wgts(2)
integer  :: domid
!
if ( .not. module_initialized ) call static_init_model

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
expected_obs(:) = MISSING_R8

! istatus for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

! Start the horizontal interpolation code here

loc       = get_location(location) ! loc is in DEGREES
loc_lon   = loc(1)
loc_lat   = loc(2)
loc_depth = loc(3)

varid = get_varid_from_kind(1, qty)
!print*, routine, varid, loc_lon, loc_lat, loc_depth

call get_corners(loc_lon, loc_lat, geo_inds, geo_wgts)

call horizontal_interpolate(state_handle, ens_size, geo_inds, geo_wgts, &
                            loc_depth, varid, expected_obs) 
istatus(:) = 0

end subroutine model_interpolate

!------------------------------------------------------------------
! Horizontal interpolation for surface variable
subroutine horizontal_interpolate(state_handle, ens_size, ij_inds, ij_wgts, zloc, vid, fld)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,           intent(in)   :: ij_inds(4)
real(r8),           intent(in)  :: ij_wgts(2), zloc
integer,             intent(in) :: vid
real(r8),           intent(out) :: fld(ens_size)
!
real(r8)     :: x_iul(ens_size) ,x_iur(ens_size)
real(r8)     :: x_ill(ens_size), x_ilr(ens_size)
real(r8)     :: ij_rwgts(2)
integer(i8)  :: ill, iul, ilr, iur
integer      :: domid

!
domid = 1

ij_rwgts(1) = 1.0_r8 - ij_wgts(1)  !ifrac
ij_rwgts(2) = 1.0_r8 - ij_wgts(2)  !jfrac
 
ill = get_dart_vector_index(ij_inds(1), ij_inds(3), int(zloc), domid, vid)
iul = get_dart_vector_index(ij_inds(1), ij_inds(4), int(zloc), domid, vid)
ilr = get_dart_vector_index(ij_inds(2), ij_inds(3), int(zloc), domid, vid)
iur = get_dart_vector_index(ij_inds(2), ij_inds(4), int(zloc), domid, vid)

x_ill = get_state(ill, state_handle)
x_ilr = get_state(ilr, state_handle)
x_iul = get_state(iul, state_handle)
x_iur = get_state(iur, state_handle)

fld   = ij_rwgts(2) * ( ij_rwgts(1) * x_ill + ij_wgts(1) * x_ilr) + &
        ij_wgts(2)  * ( ij_rwgts(1) * x_iul + ij_wgts(1) * x_iur)

end subroutine horizontal_interpolate

!------------------------------------------------------------------
! Get corers of the model grid for observed location
subroutine get_corners(ilon, ilat, ij_inds, ij_wgts, ijstatus)

real(r8), intent(in)  :: ilon, ilat
integer,  intent(out) :: ij_inds(4)   ! left, right, bottom, top
real(r8), intent(out) :: ij_wgts(2)   ! x- and y- fractions
integer,  intent(out), optional :: ijstatus      !

!
integer               :: indarr(1)
real(r8)              :: idlon, idlat, xlon(nx), ylat(ny)

! Initialize
  ij_inds = -1
  ij_wgts = 0.0_r8
  !ijstatus = 1 

  if (npts.gt.1) then     !CPS for single pt simulations

!    | ............ dlon  ............. |
!    | ..idlon .. |
!    |------------|---------------------|
! lonleft       ilon                 lonright
! ij_inds(1)                         ij_inds(2)

  xlon        = lon2d(:,1) 
  indarr      = minloc(abs(xlon-ilon))

  if (xlon(indarr(1)).gt.ilon) then
    ij_inds(1)  = indarr(1) - 1
  else
    ij_inds(1)  = indarr(1)
  endif
  ij_inds(2)  = ij_inds(1) + 1

  !For boundary
  ij_inds(1) = max(ij_inds(1),1)
  ij_inds(2) = min(ij_inds(2),nx)

  idlon       = ilon - xlon(ij_inds(1))
  ij_wgts(1)   = idlon/dlon


!
!    | ............ dlat  ............. |
!    | ..idlat .. |
!    |------------|---------------------|
! latbot         ilat                 lattop
! ij_inds(3)                         ij_inds(4)

  ylat        = lat2d(1,:)
  indarr      = minloc(abs(ylat-ilat))
  if (ylat(indarr(1)).gt.ilat) then
    ij_inds(3) = indarr(1) - 1
  else
    ij_inds(3) = indarr(1)
  endif
  ij_inds(4)   = ij_inds(3) + 1

  !For boundary
  ij_inds(3) = max(ij_inds(3),1)
  ij_inds(4) = min(ij_inds(4),ny)

  idlat        = ilat - ylat(ij_inds(3))
  ij_wgts(2)   = idlat/dlat

  else

  xlon        = lon2d(:,1)
  ylat        = lat2d(1,:)
  ij_wgts  = 0.5_r8
  ij_inds = 1

  end if

  !ijstatus = 0 
  write(*,*) 'longitude index ', ij_inds(1:2)
  write(*,*) 'latitude index ',  ij_inds(3:4)
  write(*,*) 'wgts i j ', ij_wgts(1:2)
  write(*,*) 'lonW, lon, lonE ', xlon(ij_inds(1)), ilon, xlon(ij_inds(2))
  write(*,*) 'latS, lat, latN ', ylat(ij_inds(3)), ilat, ylat(ij_inds(4))

end subroutine get_corners


!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations

function apply_clamping(ivar, posterior) result (slab)

integer,  intent(in) :: ivar
real(r8), intent(in) :: posterior(:)
real(r8)             :: slab(size(posterior))

slab = posterior

!if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
!    (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
!   where ((slab /= MISSING_R8) .and. &
!          (slab > progvar(ivar)%maxvalue)) &
!           slab = progvar(ivar)%maxvalue
!endif

!if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
!    (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
!   where ((slab /= MISSING_R8) .and. &
!          (slab < progvar(ivar)%minvalue)) &
!           slab = progvar(ivar)%minvalue
!endif

end function apply_clamping

!------------------------------------------------------------------
!Read and parse MSHM variables to be assimilated/updated

subroutine check_mvar_assim(state_variables, nfields, &
                      var_names, var_qtys, var_ranges, var_update)

character(len=*),                 intent(in)  :: state_variables(:,:)
integer,                          intent(out) :: nfields
character(len=*), intent(out) :: var_names(:)
real(r8),         intent(out) :: var_ranges(:,:)
logical,          intent(out) :: var_update(:)
integer ,         intent(out) :: var_qtys(:)

integer :: i, io, quantity
real(r8) :: minvalue, maxvalue

character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5

nfields=0
MyLoop : do i = 1, size(state_variables,2) 

   varname      = state_variables(VT_VARNAMEINDX,i)
   dartstr      = state_variables(VT_KINDINDX   ,i)
   minvalstring = state_variables(VT_MINVALINDX ,i)
   maxvalstring = state_variables(VT_MAXVALINDX ,i)
   state_or_aux = state_variables(VT_STATEINDX  ,i)

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml: variable list not fully specified'
      call error_handler(E_ERR,'check_mvar_assim', string1, source)
   endif

   ! Make sure DART kind is valid
   quantity = get_index_for_quantity(dartstr)
   if( quantity < 0 ) then
      write(string1,'(''there is no obs_kind "'',a,''" in obs_kind_mod.f90'')') &
                    trim(dartstr)
      call error_handler(E_ERR,'check_mvar_assim',string1,source)
   endif

   nfields = nfields + 1

   var_names(nfields)  = varname
   var_qtys(  nfields)   = quantity 
   var_ranges(nfields,:) = (/ MISSING_R8, MISSING_R8 /)
   var_update(nfields)   = .false.   ! at least initially

   ! convert the [min,max]valstrings to numeric values if possible
   read(minvalstring,*,iostat=io) minvalue
   if (io == 0) var_ranges(nfields,1) = minvalue

   read(maxvalstring,*,iostat=io) maxvalue
   if (io == 0) var_ranges(nfields,2) = maxvalue

   call to_upper(state_or_aux)
   if (state_or_aux == 'UPDATE') var_update(nfields) = .true.

enddo MyLoop

end subroutine check_mvar_assim

!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

!CPS
integer  :: iloc, jloc, kloc
integer  :: varid,domid,qtyid

if ( .not. module_initialized ) call static_init_model

!****** 2d and 3d variables indices differ, kloc, iloc, jloc for 3D
call get_model_variable_indices(index_in, iloc, jloc, kloc,varid,domid,qtyid)

!CPS model_mod_check shows the way to format indices
!print*, "index_in, iloc, jloc, kloc, varid,domid,qtyid"
!print*, index_in, iloc, jloc, kloc, varid,domid,qtyid
!print*, "lon-lat " , lon2d(iloc,jloc), lat2d(iloc,jloc)

location = set_location( lon2d(iloc,jloc), lat2d(iloc,jloc), 0.0_r8,VERTISSURFACE)

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(qty)) qty = qtyid  

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files
subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "MSHM")

!
call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!---------------------------------------------------------
!! Initialize mshm geolocation for dart to mshm and mshm to dart
subroutine initialize_mshm_geo()

integer        :: iunit, io
integer        :: jj, ii

! Read the DART namelist
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Set the calendar
call set_calendar_type( calendar )

! Read invariant grid information
  call read_mshm_geo(mshm_geo_file)
  allocate(lat2d(nx,ny))
  allocate(lon2d(nx,ny))

  do jj=1,ny
  do ii=1,nx
    lat2d(ii,jj) = lat(jj)
    lon2d(ii,jj) = lon(ii)
  end do
  end do

!CPS  lat2d = reshape(lat,(/ nx, ny /))
!CPS  lon2d = reshape(lon,(/ nx, ny /))

  write(*,*) ' MSHM dimensions (ny,nx,nz) are ...', ny, nx, nz
  write(*,*) ' MSHM grid resolutions are (dlat,dlon) are ...', dlat, dlon
  write(*,*) ' MSHM total grid points (npts) are ...', npts
  write(*,*) ' MSHM total number of fields  ...', mshm_sv_fields
end subroutine initialize_mshm_geo

!Read mshm geolocation file
subroutine read_mshm_geo(filename)
character(len=256),    intent(in):: filename

integer(kind=4)                :: nudat, izerr
integer(kind=4)                :: pp
character(len=256)             :: errmsg
!
  nudat   = get_unit()
  open(nudat,file=trim(filename))
  ! read the header
  read(nudat, *, iostat=izerr) ny, nx, nz, dlat, dlon
  if (izerr /=0) then
    write(errmsg,*) 'Unable to read ny nx nz dlat dlon from ',trim(filename)
    call error_handler(E_ERR,'read_mshm_geo',errmsg,source)
  end if

  ! allocate lat/lon
  npts = ny*nx
  !CPS allocate(lon(npts))
  !CPS allocate(lat(npts))
  allocate(lon(nx))
  allocate(lat(ny))
  allocate(levels(nz))

  read(nudat,*,iostat=izerr) lat(:)                   !Latitude
  if (izerr /=0) then
    errmsg = "Unable to read lat" 
    call error_handler(E_ERR,'read_mshm_geo',errmsg,source)
  end if
  read(nudat,*,iostat=izerr) lon(:)                   !Longitude
  if (izerr /=0) then
    errmsg = "Unable to read lon"
    call error_handler(E_ERR,'read_mshm_geo',errmsg,source)
  end if
  close(nudat)
end subroutine read_mshm_geo

!---------------------------------------------------------
subroutine constrain_snowdensity(srho,sswe,sdepth,csr)
! Update snow density after assimilation and constrain it to physical limits
! This means that the SWE will be updated again for mass conservation

! csr=1
! For snow addition, we update it within the lower limits of fresh snow density
! to maximum of ice density
! csr=2
! For snow removal, we retain the model snow density, so it means that the
! assimilated SWE is again corrected to the model value

real(r8),   intent(inout)  :: srho, sswe, sdepth
integer,    intent(in)     :: csr
real(r8)                   :: tsrho

! Based on fresh snow density parameterization for top layer
! at -20 C (graph plot of range for NASA SnowEX site 2017)
!
! First compute new snow density
 tsrho = sswe*1000._r8/sdepth

 if (csr.eq.1) then
  ! snow addition
  if (tsrho.lt. 68._r8) then
    srho = 68._r8
  else if (tsrho .gt. 917._r8) then
    srho = 917._r8
  else
    srho = tsrho !CPS 
  end if
 else if (csr.eq.2) then
  ! snow removal
  ! keep snow density same
 end if

 sswe = srho * sdepth /1000._r8

end subroutine constrain_snowdensity

!---------------------------------------------------------
!Read the mshm ascii restart file
subroutine read_mshm_restart(filename,stime,thrs,mvar1,mvar2,mvar3)

character(len=*),    intent(in)      :: filename
real(r8),            intent(out)     :: mvar3(mshm_sv_fields,nz,npts)
real(r8),            intent(out)     :: mvar2(3,npts)
integer,             intent(out)     :: mvar1(npts)
type(time_type),     intent(out)     :: stime
real(r8),            intent(out)     :: thrs

integer(kind=4)                :: nudat, i, pp, izerr
integer(kind=4)                :: yyyy,mm,dd,hh,mn,ss
!
  nudat   = get_unit()
  open(nudat,file=trim(filename))
  do pp = 1, npts
    if (pp.eq.1) then
      read (nudat,*,iostat=izerr) thrs,yyyy,mm,dd,hh,mn,ss !Timee 
    end if
    read(nudat,*,iostat=izerr) mvar1(pp)                  !Snow Layers [-]
    read(nudat,*,iostat=izerr) mvar3(1,:,pp)              !SWE [m]
    read(nudat,*,iostat=izerr) mvar3(2,:,pp)              !Snow Depth [m]
    read(nudat,*,iostat=izerr) mvar3(3,:,pp)              !Liquid Water [m]
    read(nudat,*,iostat=izerr) mvar3(4,:,pp)              !Density [kgm-3]
    read(nudat,*,iostat=izerr) mvar3(5,:,pp)              !Corr. Length [mm]
    read(nudat,*,iostat=izerr) mvar3(6,:,pp)        !Dendricity 
    read(nudat,*,iostat=izerr) mvar3(7,:,pp)        !Sphericity
    read(nudat,*,iostat=izerr) mvar3(8,:,pp)        !Grain Size [m]
    read(nudat,*,iostat=izerr) mvar3(9,:,pp)              !Snow Temp. [K]
    read(nudat,*,iostat=izerr) mvar3(10,:,pp)              !Porosity [-]
    read(nudat,*,iostat=izerr) mvar2(1,pp)                  !Soil Temp. [K]
    read(nudat,*,iostat=izerr) mvar2(2,pp)                  !Total SWE [m]
    read(nudat,*,iostat=izerr) mvar2(3,pp)                  !Total Snow Depth [m]
  end do
  close(nudat)
  !Convert time to DART units
  write(*,*) ' MSHM start date is  ...', yyyy, mm,dd,hh,mn,ss
  write(*,*) ' MSHM simulation hours is  ...', thrs 
  stime = set_date(yyyy,mm,dd,hh,mn,ss)

end subroutine read_mshm_restart

!---------------------------------------------------------
function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time


!---------------------------------------------------------
subroutine write_model_time(ncid, dart_time)
integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time
integer :: iyear, imonth, iday, ihour, imin, isec
integer :: iunit

call get_date(dart_time, iyear, imonth, iday, ihour, imin, isec)

write(*,*) "Model time : ", iyear, imonth, iday, ihour, imin, isec

!CPS need this with observation files 
iunit = open_file('dart_assim_time.txt',form='formatted')
write(iunit, '(i4.4,i2.2,i2.2,i2.2,i2.2,i2.2)') iyear,imonth,iday,ihour,imin,isec
call close_file(iunit)

end subroutine write_model_time

! Get Model time based on start time and hours of simulations
!---------------------------------------------------------
function read_model_time(filename)

character(len=*), intent(in) :: filename

real(r8)                            :: thrs
type(time_type)                     :: stime
type(time_type)                     :: read_model_time 
integer                             :: iyear,imonth,iday,ihour,imin,isec
integer(kind=4)                     :: tdays, tseconds
integer                             :: ncid

  if ( .not. file_exist(filename) ) then
    write(string1,*) 'cannot open file ', trim(filename),' for reading.'
    call error_handler(E_ERR,'read_model_time',string1,source)
  endif

  ncid = nc_open_file_readonly(filename,'read_model_time')
  ! Get Time variables
  call nc_get_variable(ncid, 'simulation_hours', thrs, source)
  call nc_get_variable(ncid, 'start_year',  iyear, source)
  call nc_get_variable(ncid, 'start_month',  imonth, source)
  call nc_get_variable(ncid, 'start_day',  iday, source)
  call nc_get_variable(ncid, 'start_hour',  ihour, source)
  call nc_get_variable(ncid, 'start_minute',  imin, source)
  call nc_get_variable(ncid, 'start_second',  isec, source)

  tdays = thrs/24
  tseconds = INT((thrs/24. - tdays)*86400)
  stime = set_date(iyear,imonth,iday,ihour,imin,isec)
  read_model_time = increment_time(stime, tseconds, tdays)
  
  call write_model_time(ncid, read_model_time)

  call nc_close_file(ncid)

end function read_model_time

! REad geolocation from netcdf file
!---------------------------------------------------------
subroutine extract_nc_geo(filename)

character(len=*), intent(in) :: filename

integer                      :: ncid
integer                      :: j, i
real(r8) , allocatable       :: rlon(:),rlat(:)

 if ( .not. file_exist(filename) ) then
    write(string1,*) 'cannot open file ', trim(filename),' for reading.'
    call error_handler(E_ERR,'extract_nc_geo',string1,source)
 endif

 ncid = nc_open_file_readonly(filename,'extract_nc_geo') 

 ! Get dimension size
 nx =  nc_get_dimension_size(ncid,'lon',source)
 ny =  nc_get_dimension_size(ncid,'lat',source)

 ! Set 1D vector size
 npts = ny*nx

 allocate(rlon(nx))
 allocate(rlat(ny))

 ! Get geo location
 call nc_get_variable(ncid, 'lon', rlon, source)
 call nc_get_variable(ncid, 'lat', rlat, source)

 allocate(lat2d(nx,ny))
 allocate(lon2d(nx,ny))

 do j=1,ny
   lat2d(:,j) = rlat(j)
 end do
 do i=1,nx
   lon2d(i,:) = rlon(i)
 end do

 !CPS This does not work for single point simulations and since wgts are
 ! pre-specificed to 0.5, not needed, 
 ! but we need this for 2D grids 
 !For regular grids, update dlat and dlon
 if (npts.gt.1) then
   dlat = rlat(2) - rlat(1)
   dlon = rlon(2) - rlon(1)
 end if

 call nc_close_file(ncid)
 
end subroutine extract_nc_geo

!---------------------------------------------------------
! Write out update mshm restart file
subroutine write_mshm_restart(dartfile,rdartfile,ufilename,fpartition)

character(len=256), intent(in)     :: dartfile,rdartfile,ufilename
integer(kind=4), intent(in)        :: fpartition

integer(kind=4)                :: nudat, i, pp, izerr
type(time_type)                 :: mtime
character*100           :: px, os
real(r8)                :: mvar3(mshm_sv_fields,nz,npts)
real(r8)                :: mvar2(3,npts)
real(r8)                :: pvar2(2,npts)
integer                 :: mvar1(npts)
real(r8)                :: thrs
integer                 :: iyear,imonth,iday,ihour,imin,isec
type(time_type)         :: stime
!
  if ( .not. file_exist(dartfile) ) then
    write(string1,*) 'cannot open file ', trim(dartfile),' for reading.'
    call error_handler(E_ERR,'restart_file_to_sv',string1,source)
  endif

  if ( .not. file_exist(rdartfile) ) then
    write(string1,*) 'cannot open file ', trim(rdartfile),' for reading.'
    call error_handler(E_ERR,'restart_file_to_sv',string1,source)
  endif

  ! Read dart updated netcdf file
  call read_dart_netcdf(dartfile, thrs, stime, mvar1,mvar2,mvar3)

  ! Read dart updated netcdf file
  call read_dartposterior_netcdf(rdartfile, pvar2)

  ! Repartition 2D integrated variable into model states
  call repartition_integratedvar(mvar1,mvar2,mvar3,pvar2,fpartition)

  !Get date for header of the restart file
  call get_date(stime, iyear,imonth,iday,ihour,imin,isec)

  ! Clamp variables if necessary
  !CPS if (desired) rbuf = apply_clamping(ivar, sv) 

  ! Write the ascii file 
  !npts should be nz but ok !  
  write(px,'(i5)') npts ! Transfer integer into character variable
  write(os,*) '(',trim(px),'(F15.6,1x))' ! output style, machine error,

  nudat = get_unit()
  open(nudat,file=ufilename)
  do pp = 1, npts
    if (pp.eq.1) then
      write(nudat,*)  thrs,iyear,imonth,iday,ihour,imin,isec        ! Write time variable here for future
    end if
    write(nudat,'(i2)') mvar1(pp)                   !Snow Layers [-]
    write(nudat,os) (mvar3(1,i,pp),i=1,nz)      !SWE [m]
    write(nudat,os) (mvar3(2,i,pp),i=1,nz)      !Snow Depth [m]
    write(nudat,os) (mvar3(3,i,pp),i=1,nz)      !Liquid Water [m]
    write(nudat,os) (mvar3(4,i,pp),i=1,nz)      !Density [kgm-3]
    write(nudat,os) (mvar3(5,i,pp),i=1,nz)      !Corr. Length [mm]
    write(nudat,os) (mvar3(6,i,pp),i=1,nz)   !Dendricity 
    write(nudat,os) (mvar3(7,i,pp),i=1,nz)   !Sphericity  
    write(nudat,os) (mvar3(8,i,pp),i=1,nz)   !Grain Size [m] 
    write(nudat,os) (mvar3(9,i,pp),i=1,nz)      !Snow Temp. [K]
    write(nudat,os) (mvar3(10,i,pp),i=1,nz)      !Porosity [-]
    write(nudat,os) mvar2(1,pp)                 !Soil Temp. [K]
    write(nudat,os) pvar2(1,pp)                 !Total SWE [m]
    write(nudat,os) pvar2(2,pp)                 !Total SD [m]
  end do
  close(nudat)

end subroutine write_mshm_restart

!---------------------------------------------------------
! This routine partitions updated increments in snow variables
! to the 3D variables, swe, snow_depth, density
subroutine repartition_integratedvar(mvar1,mvar2,mvar3,pvar2,fpartition)
integer, intent(inout)               :: mvar1(npts)
real(r8), intent(inout)              :: mvar2(3,npts)
real(r8), intent(inout)              :: mvar3(mshm_sv_fields,nz,npts)
real(r8), intent(in)                 :: pvar2(2,npts)
integer(kind=4), intent(in)          :: fpartition

integer               :: iz, ip
real(r8)              :: frac1, frac2
real(r8)              :: incr1, incr2, tmp1, tmp2

! SWE mvar3(1,:,:)
! Snow Depth mvar3(2,:,:)
! Liquid Water mvar3(3,:,:)
! Density mvar3(4,:,:)

 if (fpartition .eq. 1) then
 ! Use fraction profile to distribute increments vertically
   do ip = 1, npts
   !CPS write(*,*) 'repartition algorithm : DA updates: ', pvar2(1,ip),pvar2(2,ip) 
   if (pvar2(1,ip).gt.0._r8 .and. pvar2(2,ip).gt.0._r8) then       ! Use both 
     incr1 = pvar2(1,ip) - sum(mvar3(1,:,ip))       !SWE increment (+,-)
     incr2 = pvar2(2,ip) - sum(mvar3(2,:,ip))       !Snow Depth increment (+,-)
     !CPS write(*,*) 'repartition algorithm :Increments : ', incr1, incr2
     ! Scenarios:
     ! 1) DART adds snow
     if (incr1.gt.0._r8 .and. incr2.gt.0._r8) then
       ! to no snow pixel
       if (mvar1(ip).eq.0) then
          !CPS write(*,*) "repartition algorithm :1) DART add snow to no snow pixel"
          !CPS Need to initialize other varaiables as well
          mvar1(ip) = 1
          mvar3(1,1,ip) = incr1    !SWE
          mvar3(2,1,ip) = incr2    !Snow Depth
          mvar3(3,1,ip) = 0._r8    !Liquid water
          ! mvar3(4,1,ip) = mvar3(1,1,ip) * 1000._r8 / mvar3(2,1,ip) !Snow Density
          ! Update and constrain snow density 
          call constrain_snowdensity(mvar3(4,1,ip),mvar3(1,1,ip),mvar3(2,1,ip),1)

          mvar3(5,1,ip) = 0.07_r8  !Correlation length[mm]
          !CPSSM
          mvar3(6,1,ip) = 1.0_r8  !Dendricity
          mvar3(7,1,ip) = 0.5_r8  !Sphericity
          mvar3(8,1,ip) = 1.E-4_r8  !Grain Size [m]
          !CPSSM
          mvar3(9,1,ip) = 273.15_r8 ! Snow Temperature
          mvar3(10,1,ip) = 1.0_r8-(mvar3(4,1,ip)/917._r8)  !Porosity 
       else
         do iz = mvar1(ip),1,-1         !Top to Bottom Layer
           ! Fraction profile
           !CPS frac1 = mvar3(1,iz,ip)/mvar2(1,ip)     !SWE Fraction
           !CPS frac2 = mvar3(2,iz,ip)/mvar2(3,ip)     !Snow Depth Fraction
           if (iz.eq.mvar1(ip)) then
             frac1 = 1._r8
             frac2 = 1._r8
           else
             frac1 = 0._r8
             frac2 = 0._r8
           end if
           ! Update SWE and Snow Depth 
           mvar3(1,iz,ip) = mvar3(1,iz,ip) + frac1*incr1    !SWE(z)
           mvar3(2,iz,ip) = mvar3(2,iz,ip) + frac2*incr2    !Snow Depth (z)
           ! Update snow density for conservation
           ! mvar3(4,iz,ip) = mvar3(1,iz,ip) * 1000. / mvar3(2,iz,ip)   !
           ! Update and constrain snow density 
           call constrain_snowdensity(mvar3(4,iz,ip),mvar3(1,iz,ip),mvar3(2,iz,ip),1)

           mvar3(10,iz,ip) = 1.0_r8-(mvar3(4,iz,ip)/917._r8)  !Porosity CPSSM

           !CPS write(*,*) 'repartition algorithm :DART add snow to ',iz,frac1,frac2
         end do 
       end if  ! 1) 
       ! Update integrated var
       mvar2(2,ip) = pvar2(1,ip)
       mvar2(3,ip) = pvar2(2,ip)
     
     else if (incr1.le.0._r8 .and. incr2.le.0._r8) then
     ! 2) DART removes snow
       if (mvar1(ip).gt.0) then
          !CPS write(*,*) 'repartition algorithm :2) DART removes snow ...'
         do iz = mvar1(ip),1,-1         !Top to Bottom Layer
           ! Update SWE and Snow Depth
           tmp1 = mvar3(1,iz,ip)
           tmp2 = mvar3(2,iz,ip)
           mvar3(1,iz,ip) = max(tmp1 + incr1, 0._r8)    !SWE(z)
           mvar3(2,iz,ip) = max(tmp2 + incr2, 0._r8)    !Snow Depth (z)
           ! Update increments to remove snow from lower layers
           incr1 = min(tmp1 + incr1, 0._r8)
           incr2 = min(tmp2 + incr2, 0._r8) 
           ! Update snow density for conservation
           ! Use both as otherwise Infinite density possible if SWE only used
           if (mvar3(1,iz,ip).eq.0._r8 .or. mvar3(2,iz,ip).eq.0._r8) then      !Use both 
             !CPS write(*,*) 'repartition algorithm : Removing layer ',iz,incr1,incr2
             mvar3(:,iz,ip) = 0._r8    ! Set all state of layer to zero
             mvar1(ip) = mvar1(ip) -1  ! Update snow layer 
           else
             ! mvar3(4,iz,ip) = mvar3(1,iz,ip) * 1000. / mvar3(2,iz,ip)   !
             ! Update and constrain snow density 
             call constrain_snowdensity(mvar3(4,iz,ip),mvar3(1,iz,ip),mvar3(2,iz,ip),2) !-1 to 2
           end if
           mvar3(10,iz,ip) = 1.0_r8-(mvar3(4,iz,ip)/917._r8)  !Porosity CPSSM
           !CPS write(*,*)'repartition algorithm : iz',iz,incr1,incr2,mvar3(:,iz,ip)
         end do
       end if
       ! Update integrated var
       mvar2(2,ip) = pvar2(1,ip)
       mvar2(3,ip) = pvar2(2,ip)

     else
       !CPS if the signs of increments are opposite, we do not update
       !need to futher explore this situation, what it means
       !CPS write(*,*)'repartition algorithm : not updating snow variables'
     end if

   else          !CPS if (pvar2(1,ip) .le. 0._r8) then
   ! DART removes snow completely from the pixel
     !CPS write(*,*) 'repartition algorithm : Removing entire snow layer' 
     mvar1(ip)     = 0
     mvar2(2,ip)   = 0._r8
     mvar2(3,ip)   = 0._r8
     mvar3(:,:,ip) = 0._r8
   end if
   !CPS write(*,*) 'repartition algorithm : Updates ', mvar2(2,ip),mvar2(3,ip)
   end do
 
 else
   string1 = "Code not written for repartition_swe =/ 1" 
   call error_handler(E_ERR,'repartition_integratedvar', string1, source)
 end if   !repartition switch
 
end subroutine repartition_integratedvar

!!---------------------------------------------------------
!This routine reads posterior netcdf output from DART
subroutine read_dartposterior_netcdf(dfilename, pvar2)

character(len=256), intent(in)        :: dfilename
real(r8), intent(out)                 :: pvar2(2,npts)

integer                               :: ncid
real(r8)                              :: ncvar2(nx,ny)

!Read and reshape netcdf data 
  ncid = nc_open_file_readonly(dfilename)

  call nc_get_variable(ncid, 'H2OSNO', ncvar2)
  pvar2(1,:) = reshape(ncvar2,(/ npts /))

  call nc_get_variable(ncid, 'SNOW_DEPTH', ncvar2)
  pvar2(2,:) = reshape(ncvar2,(/ npts /))
!
  call nc_close_file(ncid)

end subroutine read_dartposterior_netcdf

!---------------------------------------------------------
!This routine reads prior netcdf output from DART
subroutine read_dart_netcdf(dfilename,thrs, stime, mvar1,mvar2,mvar3)

character(len=256), intent(in)        :: dfilename
real(r8), intent(out)                 :: mvar3(mshm_sv_fields,nz,npts)
real(r8), intent(out)                 :: mvar2(3,npts)
integer, intent(out)                  :: mvar1(npts)
real(r8), intent(out)                 :: thrs
type(time_type), intent(out)          :: stime

integer                               :: ncid, nslev, pp
real(r8)                              :: ncvar3(nz,nx,ny),ncvar2(nx,ny)
integer                               :: ncvar1(nx,ny)
integer                               :: iyear, imonth, iday, ihour, imin, isec


  !Read and reshape netcdf data 
  ncid = nc_open_file_readonly(dfilename)

  ! First get Time variables
  call nc_get_variable(ncid, 'simulation_hours', thrs, source)
  call nc_get_variable(ncid, 'start_year',  iyear, source)
  call nc_get_variable(ncid, 'start_month',  imonth, source)
  call nc_get_variable(ncid, 'start_day',  iday, source)
  call nc_get_variable(ncid, 'start_hour',  ihour, source)
  call nc_get_variable(ncid, 'start_minute',  imin, source)
  call nc_get_variable(ncid, 'start_second',  isec, source)
  !
  stime =  set_date(iyear, imonth, iday, ihour, imin, isec) 
 
  call nc_get_variable(ncid, 'SNOLEV', ncvar1)
  mvar1        = reshape(ncvar1, (/ npts /))

  call nc_get_variable(ncid, 'SWE', ncvar3)
  mvar3(1,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'SNODEP', ncvar3)
  mvar3(2,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'H2OLIQ', ncvar3)
  mvar3(3,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'RHOSNO', ncvar3)
  mvar3(4,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'SNOCL', ncvar3)
  mvar3(5,:,:) = reshape(ncvar3,(/ nz, npts /))
!CPSSM
  call nc_get_variable(ncid, 'SNODENDR', ncvar3)
  mvar3(6,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'SNOSPHER', ncvar3)
  mvar3(7,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'SNOGRAIN', ncvar3)
  mvar3(8,:,:) = reshape(ncvar3,(/ nz, npts /))
!CPSSM
  call nc_get_variable(ncid, 'TSNO', ncvar3)
  mvar3(9,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'PHI', ncvar3)
  mvar3(10,:,:) = reshape(ncvar3,(/ nz, npts /))

  call nc_get_variable(ncid, 'TSOI', ncvar2)
  mvar2(1,:) = reshape(ncvar2,(/ npts /))

  call nc_get_variable(ncid, 'H2OSNO', ncvar2)
  mvar2(2,:) = reshape(ncvar2,(/ npts /))

  call nc_get_variable(ncid, 'SNOW_DEPTH', ncvar2)
  mvar2(3,:) = reshape(ncvar2,(/ npts /))
!
  call nc_close_file(ncid)

! Remove missing_value and fill with zeros
! for no snow layers
  do pp =1, npts
    nslev = mvar1(pp)      !Number of snow layers
    if (nslev.eq.0) then
      mvar3(:,:,pp) = 0._r8 
    else
      mvar3(:,nslev+1:nz,pp) = 0._r8
    end if
  end do

end subroutine read_dart_netcdf

!---------------------------------------------------------
! This routine writes netcdf data for DART
subroutine write_dart_netcdf(rfilename, dfilename)

character(len=256), intent(in)        :: rfilename, dfilename
!
real(r8)                             :: mvar3(mshm_sv_fields,nz,npts)
real(r8)                             :: mvar2(3,npts)
integer                              :: mvar1(npts)
type(time_type)                      :: stime
real(r8)                             :: thrs

character(len=*), parameter :: source = 'write_dart_netcdf'
integer                               :: ncid, pp, nslev
integer                               :: snl(nx,ny)
real(r8)                              :: tss(nx,ny)
real(r8)                              :: prg(nz,nx,ny)
integer                               :: iyear,imonth,iday,ihour,imin,isec

!CPS if ( .not. module_initialized ) call static_init_model

!First read mshm restart file
call read_mshm_restart(rfilename,stime,thrs,mvar1,mvar2,mvar3)

!Get dates, we need to store thrs, and stime dates 
call get_date(stime, iyear,imonth,iday,ihour,imin,isec)

! Mask 3D data for no snow layers
do pp =1, npts
  nslev = mvar1(pp)      !Number of snow layers
  if (nslev.eq.0) then
     mvar3(:,:,pp) = missingValue
  else
     mvar3(:,nslev+1:nz,pp) = missingValue
  end if
end do

!Write dart netcdf file
ncid = nc_create_file(dfilename,source)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "MSHM")

!Dimensions
call nc_define_dimension(ncid, 'lon', nx, source)   ! Longitude 
call nc_define_dimension(ncid, 'lat', ny, source)   ! Latitude 
call nc_define_dimension(ncid, 'lev', nz, source)     ! Layers
!
!Geo-location (regular co-ordinate)
call nc_define_real_variable(ncid, 'lon', (/ 'lon' /), source)
call nc_add_attribute_to_variable(ncid, 'lon', 'long_name', 'Longitude',source)
call nc_add_attribute_to_variable(ncid, 'lon', 'units', 'degrees_east',source)
call nc_add_attribute_to_variable(ncid,'lon','valid_range',(/-180.0_r8,180.0_r8/),source)
!
call nc_define_real_variable(ncid, 'lat', (/ 'lat' /), source)
call nc_add_attribute_to_variable(ncid, 'lat', 'long_name', 'Longitude',source)
call nc_add_attribute_to_variable(ncid, 'lat', 'units', 'degrees_north',source)
call nc_add_attribute_to_variable(ncid,'lat','valid_range',(/-90.0_r8,90.0_r8/),source)
!

! Model State Variables
call nc_define_integer_variable(ncid, 'SNOLEV', (/ 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNOLEV','long_name','Snow Layers',source)
call nc_add_attribute_to_variable(ncid, 'SNOLEV', 'units', ' ', source)
!
call nc_define_real_variable(ncid, 'TSOI', (/ 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'TSOI','long_name','Soil Temperature',source)
call nc_add_attribute_to_variable(ncid, 'TSOI', 'units', 'K', source)
call nc_add_attribute_to_variable(ncid, 'TSOI','missing_value',missingValue,source)
!
call nc_define_real_variable(ncid, 'H2OSNO', (/ 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'H2OSNO','long_name','Snow Water Equivalent',source)
call nc_add_attribute_to_variable(ncid, 'H2OSNO', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid, 'H2OSNO','missing_value',missingValue,source)
!
call nc_define_real_variable(ncid, 'SNOW_DEPTH', (/ 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNOW_DEPTH','long_name','Snow Depth',source) 
call nc_add_attribute_to_variable(ncid, 'SNOW_DEPTH', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid, 'SNOW_DEPTH','missing_value',missingValue,source)
!
call nc_define_real_variable(ncid, 'SWE', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SWE','long_name','Snow Water Equivalent',source)
call nc_add_attribute_to_variable(ncid, 'SWE', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid, 'SWE','missing_value',missingValue,source)
!
call nc_define_real_variable(ncid, 'SNODEP', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNODEP','long_name','Snow Depth',source)
call nc_add_attribute_to_variable(ncid, 'SNODEP', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid, 'SNODEP','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'H2OLIQ', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'H2OLIQ','long_name','Liquid Water')
call nc_add_attribute_to_variable(ncid, 'H2OLIQ', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid, 'H2OLIQ','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'RHOSNO', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'RHOSNO','long_name','Snow Density')
call nc_add_attribute_to_variable(ncid, 'RHOSNO', 'units', 'kgm-3', source)
call nc_add_attribute_to_variable(ncid, 'RHOSNO','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'SNOCL', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNOCL','long_name','Snow Correlation Length')
call nc_add_attribute_to_variable(ncid, 'SNOCL', 'units', 'mm', source)
call nc_add_attribute_to_variable(ncid, 'SNOCL','missing_value',missingValue,source)
!CPSSM
call nc_define_real_variable(ncid, 'SNODENDR', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNODENDR','long_name','Snow Dendricity')
call nc_add_attribute_to_variable(ncid, 'SNODENDR', 'units', ' ', source)
call nc_add_attribute_to_variable(ncid,'SNODENDR','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'SNOSPHER', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNOSPHER','long_name','Snow Sphericity')
call nc_add_attribute_to_variable(ncid, 'SNOSPHER', 'units', ' ', source)
call nc_add_attribute_to_variable(ncid,'SNOSPHER','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'SNOGRAIN', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'SNOGRAIN','long_name','Snow Grain Size')
call nc_add_attribute_to_variable(ncid, 'SNOGRAIN', 'units', 'm', source)
call nc_add_attribute_to_variable(ncid,'SNOGRAIN','missing_value',missingValue,source)
!CPSSM
call nc_define_real_variable(ncid, 'TSNO', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'TSNO','long_name','Snow Temperature')
call nc_add_attribute_to_variable(ncid, 'TSNO', 'units', 'k', source)
call nc_add_attribute_to_variable(ncid, 'TSNO','missing_value',missingValue,source)

call nc_define_real_variable(ncid, 'PHI', (/ 'lev', 'lon', 'lat' /), source)
call nc_add_attribute_to_variable(ncid,'PHI','long_name','Porosity')
call nc_add_attribute_to_variable(ncid, 'PHI', 'units', 'm-3/m-3', source)
call nc_add_attribute_to_variable(ncid, 'PHI','missing_value',missingValue,source)

call nc_define_real_scalar(ncid, 'simulation_hours',source)
call nc_define_integer_scalar(ncid, 'start_year',source)
call nc_define_integer_scalar(ncid, 'start_month',source)
call nc_define_integer_scalar(ncid, 'start_day',source)
call nc_define_integer_scalar(ncid, 'start_hour',source)
call nc_define_integer_scalar(ncid, 'start_minute',source)
call nc_define_integer_scalar(ncid, 'start_second',source)

! Define model complete
call nc_end_define_mode(ncid)

! First put Time variables
call nc_put_variable(ncid, 'simulation_hours', thrs, source)
call nc_put_variable(ncid, 'start_year', 1, iyear, source)
call nc_put_variable(ncid, 'start_month',1,  imonth, source)
call nc_put_variable(ncid, 'start_day', 1, iday, source)
call nc_put_variable(ncid, 'start_hour', 1, ihour, source)
call nc_put_variable(ncid, 'start_minute', 1, imin, source)
call nc_put_variable(ncid, 'start_second', 1, isec, source)


call nc_put_variable(ncid, 'lon', lon2d(:,1), source)   !Regular grid
call nc_put_variable(ncid, 'lat', lat2d(1,:), source)   !Regular grid

! Need to reshape variables
snl = reshape(mvar1, (/ nx, ny /))
call nc_put_variable(ncid, 'SNOLEV', snl, source)

prg = reshape(mvar3(1,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SWE', prg, source)

prg = reshape(mvar3(2,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SNODEP', prg, source)

prg = reshape(mvar3(3,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'H2OLIQ', prg, source)

prg = reshape(mvar3(4,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'RHOSNO', prg, source)

prg = reshape(mvar3(5,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SNOCL', prg, source)
!CPSSM
prg = reshape(mvar3(6,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SNODENDR', prg, source)

prg = reshape(mvar3(7,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SNOSPHER', prg, source)

prg = reshape(mvar3(8,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'SNOGRAIN', prg, source)
!CPSSM
prg = reshape(mvar3(9,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'TSNO', prg, source)

prg = reshape(mvar3(10,:,:), (/ nz, nx, ny /))
call nc_put_variable(ncid, 'PHI', prg, source)

tss = reshape(mvar2(1,:), (/ nx, ny /))
call nc_put_variable(ncid, 'TSOI', tss, source)

tss = reshape(mvar2(2,:), (/ nx, ny /))
call nc_put_variable(ncid, 'H2OSNO', tss, source)

tss = reshape(mvar2(3,:), (/ nx, ny /))
call nc_put_variable(ncid, 'SNOW_DEPTH', tss, source)

call nc_close_file(ncid)

!CPS stime = set_date(iyear,imonth,iday,ihour,imin,isec)
!CPS new_Time = increment_time(stime, tseconds, tdays) 
!CPS call get_date(dart_time, iyear, imonth, iday, ihour, imin, isec)

end subroutine write_dart_netcdf

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

