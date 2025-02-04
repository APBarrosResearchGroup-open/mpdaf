! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_to_dart.f90 6311 2013-07-17 22:04:54Z thoar $

program mshm_to_dart

!----------------------------------------------------------------------
! purpose: interface between mshm and DART
!
! method: Read restart files of MSHM model state.
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
!
! 
! USAGE:  The mshm ascii filename is read from the input.nml
!         <edit mshm_to_dart_output_file in input.nml:model_to_dart_nml>
!         model_to_dart
!
! author: Prabhakar Shrestha 05/23/2023
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file,E_ERR, E_MSG, error_handler

use        model_mod, only : initialize_mshm_geo, write_dart_netcdf, &
                             read_model_time

use time_manager_mod, only : time_type, print_time, print_date

implicit none

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: mshm_to_dart_output_file  = 'dart_prior.nc'
character(len=256) :: mshm_restart_file         = 'restart_mshm'
namelist /mshm_to_dart_nml/    &
     mshm_to_dart_output_file,  &
     mshm_restart_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical                       :: verbose = .TRUE.
integer                       :: io, iunit, x_size
type(time_type)               :: model_time

character(len=512)   :: string1, string2

!======================================================================

call initialize_utilities(progname='mshm_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "mshm_to_dart_nml", iunit)
read(iunit, nml = mshm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "mshm_to_dart_nml") ! closes, too.

write(string1,*) 'converting mshm file "'//trim(mshm_restart_file)//'"'
write(string2,*) ' to DART file "'//trim(mshm_to_dart_output_file)//'"'
call error_handler(E_MSG,'mshm_to_dart',string1,text2=string2)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

call initialize_mshm_geo

call write_dart_netcdf(mshm_restart_file,mshm_to_dart_output_file)

model_time = read_model_time(mshm_to_dart_output_file)

!deallocate(x)

!----------------------------------------------------------------------
! finish up
!----------------------------------------------------------------------

call print_time(model_time, str='mshm_to_dart:DART model time')
call finalize_utilities()

end program mshm_to_dart
