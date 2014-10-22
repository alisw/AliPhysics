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

option(HLT_CUDA   "Build HLT GPU tracker using CUDA"   OFF)
option(HLT_OPENCL "Build HLT GPU tracker using OpenCL" OFF)

# AliRoot debug mode
if(ALIDEBUG)
    add_definitions(-DDEBUG -D__DEBUG)
else()
    add_definitions(-DNDEBUG)
endif (ALIDEBUG)

# MLUC
if(ALIHLT_MLUCDIR)
    add_definitions(-Duse_logging)
    include_directories(${ALIHLT_MLUCDIR}/include)
    link_directories(${ALIHLT_MLUCDIR}/lib/tgt_${ALICE_TARGET})
endif (ALIHLT_MLUCDIR) 

add_subdirectory(CALO)
add_subdirectory(comp)

#file(GLOB PACKAGES CMake*.pkg)

#if(HLT_CUDA)
#  list(APPEND PACKAGES TPCLib/tracking-ca/cagpu/CMakelibAliHLTTPCCAGPU.pkg)
#endif(HLT_CUDA)

#if(HLT_OPENCL)
#  list(APPEND PACKAGES TPCLib/tracking-ca/cagpu/CMakelibAliHLTTPCCAGPUOpenCL.pkg)
#endif()

message(STATUS "HLT enabled")
