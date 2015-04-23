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

# PMD libraries
add_subdirectory(FMDbase)
add_subdirectory(FMDutil)
add_subdirectory(FMDrec)
add_subdirectory(FMDsim)

if(DA)
    add_subdirectory(DA)
endif()

# Install files to destination
install(DIRECTORY scripts DESTINATION FMD)
install(DIRECTORY macros DESTINATION FMD)
install(FILES 
            Config.C
            MakeFMDFullMisAlignment.C
            MakeFMDResMisAlignment.C
            MakeFMDZeroMisAlignment.C
            Reconstruct.C
            Simulate.C
            ViewFMD.C
        DESTINATION FMD)

message(STATUS "FMD enabled")
