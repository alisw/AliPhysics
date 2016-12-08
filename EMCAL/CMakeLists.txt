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

# EMCAL libraries
add_subdirectory(EMCALbase)
add_subdirectory(EMCALUtils)
add_subdirectory(EMCALraw)
add_subdirectory(EMCALrec)
add_subdirectory(EMCALsim)
add_subdirectory(EMCALTriggerBase)

if(DA)
    add_subdirectory(DA)
endif(DA)

# Install the macros, docs
install(DIRECTORY macros DESTINATION EMCAL)
install(DIRECTORY mapping DESTINATION EMCAL)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc)
    install(DIRECTORY doc DESTINATION EMCAL)
elseif(FATAL_IF_MISSING)
    message(FATAL_ERROR "doc is missing")
endif()

message(STATUS "EMCAL enabled")
