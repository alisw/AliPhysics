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

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/AD/DA)

# Additional include folders in alphabetical order
include_directories(
		)
		   
set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static ADbase-static)
# Generate the DA using DetectorName, Algortihm
# The DA file has to be named DetectorNameAlgortihmda.cxx
# generateDA("DetectorName" "Algorithm" "${STATIC_DEPENDENCIES}")
generateDA("AD0" "PEDESTAL" "${STATIC_DEPENDENCIES}")
