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
# **************************************************************************

# HMPID libraries
add_subdirectory(PHOSbase)
add_subdirectory(PHOSUtils)
add_subdirectory(PHOSrec)
add_subdirectory(PHOSsim)
add_subdirectory(PHOSpi0Calib)
add_subdirectory(PHOSshuttle)

if(DA)
    add_subdirectory(DA)
endif(DA)

# Install the macros, mapping, .root, .dat
install(DIRECTORY macros DESTINATION PHOS)
install(DIRECTORY mapping DESTINATION PHOS)

set(ROOTFILES
    PCA8pa15_0.5-100.root
    PCA8pa15_0.5-5.root
    PCA_pi0_40-120.root
   )

set(DATFILES
    Parameters_0.5_100.dat
    Parameters_0.5_5.dat
    Parameters.dat
   )

if(EXISTS PHOS/PCA8pa15_0.5-100.root)
    install(FILES ${ROOTFILES} ${DATFILES} DESTINATION PHOS)
endif()

message(STATUS "PHOS enabled")
