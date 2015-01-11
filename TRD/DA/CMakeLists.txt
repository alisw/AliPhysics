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
include_directories(${AliRoot_SOURCE_DIR}/TRD/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                    ${AliRoot_SOURCE_DIR}/TPC/TPCbase
                   )

set(STATIC_DEPENDENCIES CDB-static RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "HALFCHAMBERSTATUS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES CDB-static RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "PEDESTAL" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatarec-static STEERBase-static TRDbase-static)
generateDA("TRD" "VDRIFT" "${STATIC_DEPENDENCIES}")
