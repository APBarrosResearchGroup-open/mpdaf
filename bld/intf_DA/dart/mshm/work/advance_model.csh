#!/bin/tcsh

#for perfect_model run
set fil="filter_control00000"
if ( -f ${fil} ) then
  rm ${fil}
endif

exit 0
