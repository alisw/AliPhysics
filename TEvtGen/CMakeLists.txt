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
add_subdirectory (EvtGenExternal)
add_subdirectory (EvtGen)
add_subdirectory (THepMCParser)
add_subdirectory (HepMC)
add_subdirectory (Photos)
add_subdirectory (Tauola)
add_subdirectory (TEvtGen)

message(STATUS "TEvtGen enabled")
