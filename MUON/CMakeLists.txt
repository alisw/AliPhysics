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

# MUON libraries
add_subdirectory(MUONcore)
add_subdirectory(MUONbase)
add_subdirectory(MUONcalib)
add_subdirectory(MUONcalign)
add_subdirectory(MUONevaluation)
add_subdirectory(MUONgeometry)
add_subdirectory(MUONgraphics)
add_subdirectory(MUONmapping)
add_subdirectory(MUONmpgraphics)
add_subdirectory(MUONraw)
add_subdirectory(MUONrec)
add_subdirectory(MUONshuttle)
add_subdirectory(MUONsim)
add_subdirectory(MUONtrigger)
add_subdirectory(mchview)

if(DA)
    add_subdirectory(DA)
endif(DA)

# Install data folder
install(DIRECTORY data macros mapping DESTINATION MUON)

message(STATUS "MUON enabled")
