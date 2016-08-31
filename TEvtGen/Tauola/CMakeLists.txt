# **************************************************************************
# * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
# *                                                                        *
# * Author: The ALICE Off-line Project.                                    *
# * Contributors are mentioned in the code where appropriate.              *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

#Module
set(MODULE Tauola)
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/TEvtGen/${MODULE})

# Additional includes - alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/TEvtGen
)

# Sources - alphabetical order
set(SRCS
  DecayList.cxx
  TauolaLog.cxx
  Plots.cxx
  Tauola.cxx
  TauolaEvent.cxx
  TauolaHEPEVTEvent.cxx
  TauolaHEPEVTParticle.cxx
  TauolaHepMCEvent.cxx
  TauolaHepMCParticle.cxx
  TauolaParticle.cxx
  TauolaParticlePair.cxx
  f_Decay.cxx
  f_FilHep.cxx
  f_InitTauola.cxx
)

set ( FSRCS
  tauola-fortran/curr_cleo.f
  tauola-fortran/f3pi.f
  tauola-fortran/formf.f
  tauola-fortran/pkorb.f
  tauola-fortran/tauola.f
  tauola-fortran/new-currents/other-currents/frho_pi_belle.f
  tauola-fortran/new-currents/RChL-currents/rcht_3pi/f3pi_rcht.f
  tauola-fortran/new-currents/RChL-currents/rcht_3pi/funct_3pi.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/FA1RCHL.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/ffwid3pi.f 
  tauola-fortran/new-currents/RChL-currents/rcht_common/funct_rpt.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/gaus_integr.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/gfact.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/initA1Tab.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/initA1TabKKpi.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/value_parameter.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/wid_a1_fit.f
  tauola-fortran/new-currents/RChL-currents/rcht_common/wid_a1_fitKKpi.f
  tauola_extras.f
)

# No dictionary

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} ${FSRCS})

# Do not generate the ROOT map
# Dependencies
set(LIBDEPS HepMC )

# Additional compilation flags
if(${CMAKE_Fortran_COMPILER} MATCHES gfortran)
set( CMAKE_Fortran_FLAGS "-fno-automatic -fno-backslash -ffixed-line-length-132")
endif(${CMAKE_Fortran_COMPILER} MATCHES gfortran)
if(${CMAKE_Fortran_COMPILER} MATCHES ifort)
set ( CMAKE_Fortran_FLAGS "-extend-source")
endif(${CMAKE_Fortran_COMPILER} MATCHES ifort)
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE} 
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
