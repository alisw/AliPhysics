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
# Do not forget the ending / from macros/
# The / will point to install the content of the folder to 
# macros/ITS/UPGRADE Without the / it will create a new folder 
# macros/ITS/UPGRADE/macros
install(DIRECTORY macros/ DESTINATION macros/ITS/UPGRADE)

message(STATUS "ITS UPGRADE enabled")
