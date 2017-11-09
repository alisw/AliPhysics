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
add_subdirectory(CaloTrackCorrBase)
add_subdirectory(DevNanoAOD)
add_subdirectory(EMCAL)
add_subdirectory(FLOW)
add_subdirectory(Glauber)
add_subdirectory(mftmuondep)
add_subdirectory(muon)
add_subdirectory(muondep)
add_subdirectory(Tools)
add_subdirectory(CALO)
add_subdirectory(TRD)
add_subdirectory(HMTF)
add_subdirectory(MCLEGO)
add_subdirectory(JETFW)
add_subdirectory(Cocktail)

message(STATUS "PWG enabled")
