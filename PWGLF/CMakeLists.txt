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

# Include the libraries
add_subdirectory (FORWARD)
add_subdirectory (QATasks)
add_subdirectory (RESONANCES)
add_subdirectory (RESONANCES/extra)
add_subdirectory (SPECTRA)
add_subdirectory (spectraJET)
add_subdirectory (STRANGENESS)
add_subdirectory (ThermalFits)
add_subdirectory (totEt)
add_subdirectory (ppVsMult)
add_subdirectory (NUCLEX)
add_subdirectory (CommonUtils)


message(STATUS "PWGLF enabled")
