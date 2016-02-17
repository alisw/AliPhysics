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

# Module
set(MODULE Utilities)
add_definitions(-D_MODULE_="${MODULE}")
# Aditional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIR} 
                   )

# Create an object to be reused in case of static libraries 
# Otherwise the sources will be compiled twice
# Add a library to the project using the object
add_executable(alihadd alihadd.cxx)
add_executable(alisync alisync.cxx)
target_link_libraries(alihadd RIO Hist Core)
target_link_libraries(alisync RIO Hist Core Net Tree)
 
# Installation
install(TARGETS alihadd RUNTIME DESTINATION bin)
install(TARGETS alisync RUNTIME DESTINATION bin)
