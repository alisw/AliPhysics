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

# Include the libraries
add_subdirectory(blastwave)
add_subdirectory(CME)
add_subdirectory(CME/TwoPwrtRP)
add_subdirectory(Forward)

# Install the macros
install(DIRECTORY Forward/macros DESTINATION PWGCF/FLOW)
install(DIRECTORY macros DESTINATION PWGCF/FLOW)

# Install calibration data
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/database)
  install(DIRECTORY database DESTINATION PWGCF/FLOW)
endif()
