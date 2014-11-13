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

# MONITOR libraries
add_subdirectory(MONITOR)

# MONITOR executables
add_subdirectory(alitestproc)
add_subdirectory(deroot)
add_subdirectory(monitorCheck)
add_subdirectory(monitorGDC)

#if(ZEROMQ_FOUND)
    add_subdirectory(MONITORzmq)
#  alistoragemanager/CMakebinalistorage.pkg
#  alistoragemanager/CMakebinalistorageadmin.pkg
#  alistoragemanager/CMakebinalifakedisplay.pkg
#  alistoragemanager/CMakebinalistorageminimalexample.pkg
#  alistoragemanager/CMakebinalistoragegetpoints.pkg
#  alieventserver/CMakebinalieventserver.pkg
#endif(ZEROMQ_FOUND)


message(STATUS "MONITOR enabled")
