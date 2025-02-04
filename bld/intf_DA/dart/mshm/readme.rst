mshm
==============

.. attention::
    Add your model documentation here.
1) cd build_templates
   cp  mkmf.template.gfortran mkmf.template  !Copy closest template for DART compilation 
2) Create DART template for porting mshm to DART
  cd DART/models
 ./new_model.sh mshm threed_sphere
3) Add function and  subroutines to read/write mshsm restart ascii files in model_mod.f90
   ++ read_model_time -reads simulation-hours and yyyy mm dd hh mm ss from netcdf file and converts to model_time
   ++ extract_nc_geo -reads geo location from netcdf file
   ++ check_mvar_assim - read and parse MSHM variables to be assimilated
   ++ get_state_meta_data - given index, obtain location and qty
   ++ mshm_to_dart   - converts mshm restart file to DART netcdf and updates missingValue for empty snow layers
         |
         initialize_mshm_geo  - reads grid sizes and geolocation from mshm invariant file
         write_dart_netcdf    - writes dart netcdf file
             |
             read_mshm_restart  - reads mshm ascii restart file
   ++ dart_to_mshm  - converts DART posterior netcdf to mshm ascii restart file, missingValue replace by zeros    
         |
         initialize_mshm_geo  - ''
         write_mshm_restart   - writes mshm restart file
             |
             read_dart_posterior_netcdf  - reads dart posterior netcdf
             read_dart_netcdf - reads dart prior netcdf
   Note: initialize_mshm_geo was added as an alternative to static_init, just to transfer ascii to netcdf and vice versa
         static_init directl reads time and model size from netcdf file
         Also check the namelist update in the input_nml or in the model_mod.f90

   ++ model_interpolate      - perform horizontal interpolation for given qty and location
   ++ horizontal_interpolate  - use corners and wgts for interpolation
   ++ get_corners       - get corners and wgts for location provided
   ++ check_mvar_assim  - read and parse MSHM variables to be assimilated/updated
   
4) input_nml
  ++ obs_type_files is undefined in preprocess_nml to use the state variables as direct observation type and keep 
     the FO as unity. DART differentiates between observation types and physical quantities. 
    
