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
set(MODULE PWGMCLEGO)
add_definitions(-D_MODULE_="${MODULE}")

# install all executable scripts
install(DIRECTORY .
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.sh"
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
  )

# install ThePEG configurations
install(DIRECTORY ThePEG
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.in"
  )

# install JEWEL configurations
install(DIRECTORY JEWEL
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.dat"
  )

# install CRMC configurations
install(DIRECTORY CRMC
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.param"
  )

# install EPOS configurations
install(DIRECTORY EPOS
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.optns"
  )

# install TunedOnPbPb configurations
install(DIRECTORY TunedOnPbPb
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.C"
  )

# install Herwig configurations
install(DIRECTORY Herwig
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.in"
  )

# install AMPT configurations
install(DIRECTORY AMPT
  DESTINATION PWG/MCLEGO
  FILES_MATCHING PATTERN "*.ampt"
  )
