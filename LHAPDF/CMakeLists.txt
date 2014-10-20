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
set(MODULE lhapdf)

# Module include folder
include_directories(${CMAKE_SOURCE_DIR}/LHAPDF)

# Additional include folders in alphabetical order except ROOT
include_directories(${ROOT_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/LHAPDF/lhapdf5.5.1
                    ${CMAKE_SOURCE_DIR}/LHAPDF/lhapdf5.5.1/include
                    ${CMAKE_SOURCE_DIR}/LHAPDF/lhapdf5.5.1/src
                   )

# Sources in alphabetical order
set(SRCS
    AliStructFuncType.cxx 
    lhapdf5.5.1/src/version.cxx  
    lhapdf5.5.1/src/getdatapath.cxx
   )

# Headers
set(HDRS
    AliStructFuncType.h
   )

# C sources
set(CSRCS
    lhapdf5.5.1/src/binreloc.c
   )

# Fortran sources
set(FSRCS
    lhapdf5.5.1/src/alphas.f
    lhapdf5.5.1/src/commoninit.f
    lhapdf5.5.1/src/description.f
    lhapdf5.5.1/src/eks98.f
    lhapdf5.5.1/src/eksarp.f
    lhapdf5.5.1/src/eps08.f
    lhapdf5.5.1/src/eps09.f
    lhapdf5.5.1/src/EVLCTEQ.f
    lhapdf5.5.1/src/evolution.f
    lhapdf5.5.1/src/inputPDF.f
    lhapdf5.5.1/src/LHpdflib.F
    lhapdf5.5.1/src/parameter.f
    lhapdf5.5.1/src/QCDNUM.f
    lhapdf5.5.1/src/QCDparams.f
    lhapdf5.5.1/src/Sqcdnum.f
    lhapdf5.5.1/src/wrapa02.f
    lhapdf5.5.1/src/wrapa02m.f
    lhapdf5.5.1/src/wrapabfkwpi.f
    lhapdf5.5.1/src/wrapacfgpg.f
    lhapdf5.5.1/src/wrapcteq5.f
    lhapdf5.5.1/src/wrapcteq6.f
    lhapdf5.5.1/src/wrapcteq6lg.f
    lhapdf5.5.1/src/wrapdgg.f
    lhapdf5.5.1/src/wrapdog.f
    lhapdf5.5.1/src/wrapEVLCTEQ.f
    lhapdf5.5.1/src/wrapevolve.f
    lhapdf5.5.1/src/wrapgrv.f
    lhapdf5.5.1/src/wrapgrvg.f
    lhapdf5.5.1/src/wrapgrvpi.f
    lhapdf5.5.1/src/wrapgsg96.f
    lhapdf5.5.1/src/wrapgsg.f
    lhapdf5.5.1/src/wraph1.f
    lhapdf5.5.1/src/wraplacg.f
    lhapdf5.5.1/src/wrapowpi.f
    lhapdf5.5.1/src/wrapmrst.f
    lhapdf5.5.1/src/wrapQCDNUM3.f
    lhapdf5.5.1/src/wrapQCDNUM4.f
    lhapdf5.5.1/src/wrapQCDNUM.f
    lhapdf5.5.1/src/wrapsasg.f
    lhapdf5.5.1/src/wrapsmrspi.f
    lhapdf5.5.1/src/wrapwhitg.f
    lhapdf5.5.1/src/wrapzeus.f
    lhapdf5.5.1/src/lhaglue.f
    lhapdf5.5.1/src/wrapUSER.f
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS Core)
generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library(${MODULE} SHARED ${SRCS} G__${MODULE}.cxx ${CSRCS} ${FSRCS})

# Additional compilation flags
set(CMAKE_C_FLAGS  ${CMAKE_C_FLAGS} "-ansi")
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "-O -g" Fortran_FORMAT FREE )
add_definitions(-DDEFAULTPREFIXPATH="NONE"  -DDEFAULTLHAPATH="${CMAKE_SOURCE_DIR}/LHAPDF/PDFsets")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
install(DIRECTORY PDFsets DESTINATION LHAPDF)
