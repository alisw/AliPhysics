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

# base subdirectory
add_subdirectory(base)
# liblhapdf.so
add_subdirectory(lhapdf5.5.1)
# liblhapdf-5.9.1
add_subdirectory(lhapdf-5.9.1)

if(EXISTS "LHAPDF/PDFsets")
    install(DIRECTORY PDFsets DESTINATION LHAPDF)
endif()

message(STATUS "LHAPDF enabled. Default liblhapdf.so -> 5.5.1. liblhapdf-5.9.1.so available.")
