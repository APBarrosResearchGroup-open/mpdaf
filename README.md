### MPDAF:
Multi-physics data-assimilation framework to assimilate measurements from
existing and upcoming remote sensing observations at high-spatial resolutions.

Team:

MPDAF was developed at UIUC within the project "Satellite Remote Sensing of Snow Water Resources 
and Hydrologic Prediction" funded  by the National Oceanic and Atmospheric Administration (NOAA), awarded to the 
Cooperative Institute for Research on Hydrology (CIROH) through the NOAA Cooperative Agreement with The University of Alabama, NA22NWS4320003. 

Prabhakar Shrestha, Ana P. Barros 
Department of Civil and Environmental Engineering, UIUC
email:  barros@illinois.edu

Citation:

Shrestha, P., & Barros, A. P. (2025). Multi‐ physics data assimilation framework for remotely sensed Snowpacks to improve water prediction. Water Resources Research, 61, e2024WR037885. https://doi.org/10.1029/2024WR037885

#### Usage:

MPDAF compiles and setups the run using the following two steps:

* ./build.ksh 0/1        
  >0 - Clean Compilation (For first time always)

  >1 - Complie by changing source code in tmp_XXXX directory
* ./setup.ksh            
  >Setup the defined experiment case (e.g., SnowEx in CCI cluster)
                       
#### Descriptions:
* setups/               contains different test cases for MSHM
* machines/             setups for different machines
* intf_DA/              main interface for data assimilation with included models 
* docs/                 tutorial to conduct DA using MSHM and DART
