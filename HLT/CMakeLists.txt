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

add_definitions(-Duse_aliroot -DUSE_ROOT -DROWHOUGHPARAMS)

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

add_subdirectory(BASE)
add_subdirectory(CALO)
add_subdirectory(comp)
add_subdirectory(global)
add_subdirectory(EMCAL)
add_subdirectory(FMD)
add_subdirectory(ITS)
#add_subdirectory(FJWrapper) - not needed
add_subdirectory(MUON)
add_subdirectory(pendolino)
add_subdirectory(PHOS)
add_subdirectory(RCU)
add_subdirectory(SampleLib)
add_subdirectory(TPCLib)
add_subdirectory(TRD)
add_subdirectory(trigger)
add_subdirectory(VZERO)
add_subdirectory(ZDC)
# Enable EVE if ROOT has OpenGL support
if(ROOT_HASOPENGL)
    add_subdirectory(EVE)
endif()
add_subdirectory(QA)
add_subdirectory(rec)
add_subdirectory(shuttle)
add_subdirectory(sim)

include("CMakehltdoc.conf")
install(DIRECTORY macros DESTINATION HLT)

#install macros
install(DIRECTORY exa DESTINATION HLT)

message(STATUS "HLT enabled")
