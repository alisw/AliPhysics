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

# include the libraries
add_subdirectory (STEERBase)
add_subdirectory (ESD)
add_subdirectory (AOD)
add_subdirectory (CDB)
add_subdirectory (STEER)

# Macros
# Do not forget the ending / from macros/
# The / will point to install the content of the folder to 
# macros/STEER Without the / it will create a new folder 
# macros/STEER/macros
install(DIRECTORY macros/ DESTINATION macros/STEER)

message(STATUS "STEER enabled")
