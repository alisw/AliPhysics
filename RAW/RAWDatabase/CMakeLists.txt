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

# Mandatory fields, set empty if not needed
# MODULE - target name
# MODULE_SRCS - list of sources that belong to the library
# MODULE_HDRS - list of headers that belong to the library
# MODULE_HDRS_INSTALL - list of headers that will be installed
# MODULE_INCLUDES - list of include directories
#                   ROOT and current folder are taken by default
#                   To add extra include folders use either 
#                       - include_directories(${MODULE} FOLDER)
#                       - target_include_directories(${MODULE}-object FOLDER)
# MODULE_ROOT_DEPENDENCIES
# MODULE_ALIROOT_DEPENDENCIES
# MODULE_COMPILE_FLAGS - compilation flags for the ${MODULE}-object target
# MODULE_LINK_FLAGS - linking flags for ${MODULE} target
#                   - other flags accesible thorough standard cmake procedures
#
# Target/directory properties are accesible thorough cmake standard procedures
# Target object name = ${MODULE}-object
# Target library name = ${MODULE}
# Shared and static library targets are using the same object target

###########################################################################
# Library description
###########################################################################

# Module name
set(MODULE RAWDatabase)

# Sources in alphabetical order
set(MODULE_SRCS
    AliAltroMapping.cxx
    AliCaloAltroMapping.cxx
    AliRawDataArray.cxx
    AliRawDataArrayV2.cxx
    AliRawData.cxx
    AliRawEquipment.cxx
    AliRawEquipmentHeader.cxx
    AliRawEquipmentV2.cxx
    AliRawEvent.cxx
    AliRawEventHeaderBase.cxx
    AliRawEventTag.cxx
    AliRawEventV2.cxx
    AliRawVEquipment.cxx
    AliRawVEvent.cxx
    AliStats.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" MODULE_HDRS "${MODULE_SRCS}")

# Extra headers in alphabetical order
set(MODULE_HDRS ${MODULE_HDRS}
    AliRawDataHeader.h
    AliRawDataHeaderV3.h
    AliRawEventHeaderVersions.h
   )

# List of headers to be installed
set(MODULE_HDRS_INSTALL ${MODULE_HDRS})

# List of include folders
set(MODULE_INCLUDES
        STEER/ESD
        STEER/STEERBase
   )

# Root dependecies - please solve all symbols by searching them inside the *.rootmap files
set(MODULE_ROOT_DEPENDENCIES Core Hist RIO)
# AliRoot dependencies
set(MODULE_ALIROOT_DEPENDENCIES ESD STEERBase)

# Compilation and linking flags
set(MODULE_COMPILE_FLAGS)
set(MODULE_LINK_FLAGS)

if(DATE_FOUND)
    set(MODULE_COMPILE_FLAGS "${DATE_CFLAGS}")
    set(MODULE_LINK_FLAGS "${DATE_LDFLAGS} ${DATE_LIBS}")
endif(DATE_FOUND)

# Additional compilation and linking flags
set(MODULE_COMPILE_FLAGS " ${MODULE_COMPILE_FLAGS}")

# System dependent: Modify the way the library is build
if(${CMAKE_SYSTEM} MATCHES Darwin)
    set(MODULE_LINK_FLAGS "-undefined dynamic_lookup ${MODULE_LINK_FLAGS}")
endif(${CMAKE_SYSTEM} MATCHES Darwin)

###########################################################################
# Generating the libraries
###########################################################################
# Shared librarys
generate_shared_library()

# Static library
if(ALIROOT_STATIC)
    generate_static_library()
endif(ALIROOT_STATIC)
