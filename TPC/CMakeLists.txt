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

add_subdirectory(TPCbase)
add_subdirectory(TPCutil)
add_subdirectory(TPCcalib)
add_subdirectory(TPCrec)
add_subdirectory(TPCsim)
add_subdirectory(TPCmon)
add_subdirectory(TPCupgrade)

if(DA)
    add_subdirectory(DA)
endif(DA)

# Installing the mapping
install(DIRECTORY mapping DESTINATION TPC)

# Installing the macros
install(DIRECTORY macros DESTINATION TPC)

# Install root & dat files

# The following two data files are required (ALIROOT-7111)
install(FILES AliTPCprf2dGEM.root AliTPCprf2d.root DESTINATION TPC)

# The following data files are optional
foreach(DATA_FILE CovMatrixDB_PbPb6000_B0.4T.root CovMatrixDB_pp_B0.4T.root conn_iroc.dat conn_oroc.dat)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${DATA_FILE})
    install(FILES ${DATA_FILE} DESTINATION TPC)
  elseif(FATAL_IF_MISSING)
    message(FATAL_ERROR "${DATA_FILE} is missing")
  endif()
endforeach()
