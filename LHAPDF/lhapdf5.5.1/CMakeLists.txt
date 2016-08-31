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
add_definitions(-D_MODULE_="${MODULE}")

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/LHAPDF/base
		    ${AliRoot_SOURCE_DIR}/LHAPDF/lhapdf5.5.1
                    ${AliRoot_SOURCE_DIR}/LHAPDF/lhapdf5.5.1/include
                    ${AliRoot_SOURCE_DIR}/LHAPDF/lhapdf5.5.1/src
                   )

# Additional include folders in alphabetical order
include_directories(${ROOT_INCLUDE_DIR}
                   )

# Sources in alphabetical order
set(SRCS
    src/version.cxx  
    src/getdatapath.cxx
   )

# Headers
set(HDRS
   )

# C sources
set(CSRCS
    src/binreloc.c
   )

# Fortran sources
set(FSRCS
    src/alphas.f
    src/commoninit.f
    src/description.f
    src/eks98.f
    src/eksarp.f
    src/eps08.f
    src/eps09.f
    src/EVLCTEQ.f
    src/evolution.f
    src/inputPDF.f
    src/LHpdflib.F
    src/parameter.f
    src/QCDNUM.f
    src/QCDparams.f
    src/Sqcdnum.f
    src/wrapa02.f
    src/wrapa02m.f
    src/wrapabfkwpi.f
    src/wrapacfgpg.f
    src/wrapcteq5.f
    src/wrapcteq6.f
    src/wrapcteq6lg.f
    src/wrapdgg.f
    src/wrapdog.f
    src/wrapEVLCTEQ.f
    src/wrapevolve.f
    src/wrapgrv.f
    src/wrapgrvg.f
    src/wrapgrvpi.f
    src/wrapgsg96.f
    src/wrapgsg.f
    src/wraph1.f
    src/wraplacg.f
    src/wrapowpi.f
    src/wrapmrst.f
    src/wrapQCDNUM3.f
    src/wrapQCDNUM4.f
    src/wrapQCDNUM.f
    src/wrapsasg.f
    src/wrapsmrspi.f
    src/wrapwhitg.f
    src/wrapzeus.f
    src/lhaglue.f
    src/wrapUSER.f
   )

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
# get_directory_property(incdirs INCLUDE_DIRECTORIES)
# generate_dictionary("${MODULE}" "${MODULE}LinkDef.h" "${HDRS}" "${incdirs}")

# Generate the ROOT map
# Dependecies
set(LIBDEPS Core lhapdfbase)
#generate_rootmap("${MODULE}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}LinkDef.h")

# Add a library to the project using the specified source files
add_library_tested(${MODULE} SHARED ${SRCS} ${CSRCS} ${FSRCS})
add_definitions(-DDEFAULTPREFIXPATH="NONE"  -DDEFAULTLHAPATH="${AliRoot_SOURCE_DIR}/LHAPDF/PDFsets")

# Additional compilation and linking flags
set(CMAKE_C_FLAGS "-ansi ${CMAKE_C_FLAGS}")
set(MODULE_COMPILE_FLAGS "")
set(MODULE_LINK_FLAGS)

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

# Target properties
set_target_properties(${MODULE} PROPERTIES FORTRAN_FLAGS "")
set_target_properties(${MODULE} PROPERTIES COMPILE_FLAGS "${MODULE_COMPILE_FLAGS}")
set_target_properties(${MODULE} PROPERTIES Fortran_FORMAT FREE)
set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "${MODULE_LINK_FLAGS}")

# Linking the library
target_link_libraries(${MODULE} ${LIBDEPS})

# Installation
install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)
