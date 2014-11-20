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

# AMORE - used to generate the DAs
# Flags are filled using amore-config
#               - AMORE_VERSION
#               - AMORE_STATICLIBS - libraries and linking folders for static build
#               - AMORE_AUXLIBS - auxiliary libraries for DA 
#               - AMORE_DEFINITIONS

set(AMORE_FOUND FALSE)

if(AMORE_CONFIG)
    # Checking AMORE version
    execute_process(COMMAND ${AMORE_CONFIG} --version OUTPUT_VARIABLE AMORE_VERSION ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving AMORE version : ${error}")
    endif(error)
    string(STRIP ${AMORE_VERSION} AMORE_VERSION)

    # Extract major, minor, and patch versions from
#    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+" "\\1" DATE_VERSION_MAJOR "${DATE_VERSION}")
#    string(REGEX REPLACE "^[0-9]+\\.([0-9]+)" "\\1" DATE_VERSION_MINOR "${DATE_VERSION}")
    message(STATUS "AMORE version ${AMORE_VERSION} found.")

    # Checking AMORE static libraries
    execute_process(COMMAND ${AMORE_CONFIG} --ldflags-da-static OUTPUT_VARIABLE AMORE_STATICLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving AMORE static libraries : ${error}")
    endif(error)
    string(STRIP ${AMORE_STATICLIBS} AMORE_STATICLIBS)
    string(REPLACE "\n" " " AMORE_STATICLIBS ${AMORE_STATICLIBS})

    # Checking AMORE auxiliary libraries
    execute_process(COMMAND ${AMORE_CONFIG} --auxlibs OUTPUT_VARIABLE AMORE_AUXLIBS ERROR_VARIABLE error OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(error)
        message(FATAL_ERROR "Error retrieving AMORE auxiliary libraries : ${error}")
    endif(error)
    string(STRIP ${AMORE_AUXLIBS} AMORE_AUXLIBS)

    set(AMORE_DEFINITIONS "-DALI_AMORE")
    set(AMORE_FOUND TRUE)

endif(AMORE_CONFIG)
