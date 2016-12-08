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

# Libraries
add_subdirectory(ITSUpgradeBase)
add_subdirectory(ITSUpgradeRec)
add_subdirectory(ITSUpgradeSim)

# Macros
install(DIRECTORY macros FT1 FT2 itsuTestBench DESTINATION ITSMFT/ITS)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/misc)
    install(DIRECTORY misc DESTINATION ITSMFT/ITS)
elseif(FATAL_IF_MISSING)
    message(FATAL_ERROR "misc is missing")
endif()

message(STATUS "ITS UPGRADE from ITSMFT enabled")
