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
include_directories(${AliRoot_SOURCE_DIR}/FMD/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/FMD/
                    ${AliRoot_SOURCE_DIR}/FMD/FMDbase
                    ${AliRoot_SOURCE_DIR}/FMD/FMDutil
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

set(STATIC_DEPENDENCIES FMDbase-static FMDutil-static RAWDatarec-static STEERBase-static)
generateDA("FMD" "Base" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES FMDbase-static FMDutil-static)
generateDA("FMD" "Gain" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES  FMDbase-static FMDutil-static)
generateDA("FMD" "Pedestal" "${STATIC_DEPENDENCIES}")
