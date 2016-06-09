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
#               - AMORE_CFLAGS
#               - AMORE_INCLUDE_DIR

set(AMORE_FOUND FALSE)

if(AMORE_CONFIG)

    # Helper to execute amore-config with the given options in the DATE environment.
    # Output is sent to _OUTVAR.
    # Errors are fatal.
    macro(amore_config _OPTS _OUTVAR)
        execute_process(COMMAND sh -c "source ${DATE_ENV} && export ROOTSYS=${ROOTSYS} && ${AMORE_CONFIG} ${_OPTS}"
                        OUTPUT_VARIABLE _OUTVAR_RAW
                        ERROR_VARIABLE _ERR
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(_ERR)
            message(FATAL_ERROR "Error executing ${DATE_CONFIG} ${_OPTS}")
        endif()
        string(STRIP "${_OUTVAR_RAW}" _OUTVAR_STRIPPED)
        string(REPLACE "\n" " " ${_OUTVAR} "${_OUTVAR_STRIPPED}")
    endmacro()

    # Checking AMORE version
    amore_config(--version AMORE_VERSION)
    message(STATUS "AMORE version ${AMORE_VERSION} found")

    # Checking AMORE static libraries
    amore_config(--ldflags-da-static AMORE_STATICLIBS)
    message(STATUS "AMORE static libraries: ${AMORE_STATICLIBS}")

    # Checking AMORE auxiliary libraries
    amore_config(--auxlibs-list AMORE_AUXLIBS)
    message(STATUS "AMORE auxiliary libraries: ${AMORE_AUXLIBS}")

    # Checking AMORE cflags
    amore_config(--cflags AMORE_CFLAGS)
    message(STATUS "AMORE CFLAGS: ${AMORE_CFLAGS}")

    # Checking AMORE cflags
    amore_config(--include-dir AMORE_INCLUDE_DIR)
    message(STATUS "AMORE include directory: ${AMORE_INCLUDE_DIR}")

    set(AMORE_DEFINITIONS "-DALI_AMORE")
    set(AMORE_FOUND TRUE)

endif(AMORE_CONFIG)
