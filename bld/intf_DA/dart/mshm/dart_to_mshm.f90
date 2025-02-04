! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: dart_to_model.f90 6311 2013-07-17 22:04:54Z thoar $

program dart_to_mshm

!----------------------------------------------------------------------
! purpose: interface between DART and the MSHM model
!
! method: Read DART state vector and overwrite values in a model restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called model_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance model to the requested time.
!
!         The dart_to_mshm_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Prabhakar Shrestha 05/23/2023
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file

use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date

use        model_mod, only : initialize_mshm_geo, write_mshm_restart, &
                             read_model_time

implicit none

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 256) :: dart_to_mshm_input_file     = 'dart_prior.nc'
logical               :: advance_time_present        = .false.
character(len=256)    :: mshm_up_restart_file        = 'mshm_up_restart'
character(len=256)    :: repartition_analysis_file   = 'dart_posterior.nc'
integer               :: repartition_swe             = 1

namelist /dart_to_mshm_nml/  dart_to_mshm_input_file, &
                            advance_time_present,     &
                            mshm_up_restart_file,     &
                            repartition_analysis_file, &
                            repartition_swe

!----------------------------------------------------------------------

integer               :: iunit, io, x_size, diff1, diff2
type(time_type)       :: model_time, adv_to_time, base_time
logical               :: verbose              = .FALSE.
!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_mshm', output_flag=verbose)

!----------------------------------------------------------------------
! Read model namelist and set grid sizes
!----------------------------------------------------------------------

call initialize_mshm_geo

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_mshm_nml", iunit)
read(iunit, nml = dart_to_mshm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_mshm_nml")

write(*,*)
write(*,*) 'dart_to_mshm: converting DART file ', "'"//trim(repartition_analysis_file)//"'"
write(*,*) 'to model restart files named        ', "'"//trim(mshm_up_restart_file)//"'" 

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

!if ( advance_time_present ) then
!else
!endif
!CPS call close_restart(iunit)

!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

write(*,*) 'Obtaining mshm state vector'
write(*,*) 'Calling write_mshm_file to restart file'
call write_mshm_restart(dart_to_mshm_input_file, repartition_analysis_file, mshm_up_restart_file, repartition_swe)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

model_time = read_model_time(dart_to_mshm_input_file)

call print_date( model_time,'dart_to_mshm:MSHM  model date')
call print_time( model_time,'dart_to_mshm:DART model time')
call print_date( model_time,'dart_to_mshm:MSHM model date',logfileunit)
call print_time( model_time,'dart_to_mshm:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_mshm:advance_to time')
call print_date(adv_to_time,'dart_to_mshm:advance_to date')
call print_time(adv_to_time,'dart_to_mshm:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_mshm:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_mshm
