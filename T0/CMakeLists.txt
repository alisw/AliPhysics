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

add_subdirectory(T0base)
add_subdirectory(T0calib)
add_subdirectory(T0rec)
add_subdirectory(T0sim)
add_subdirectory(T0shuttle)

if(DA)
    add_subdirectory(DA)
endif(DA)

# Install macros
install(FILES AddTaskT0Analysis.C DESTINATION T0)
install(DIRECTORY macros data DESTINATION T0)

message(STATUS "T0 enabled")
