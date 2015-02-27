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
include_directories(${AliRoot_SOURCE_DIR}/PHOS/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSbase
                    ${AliRoot_SOURCE_DIR}/PHOS/PHOSshuttle
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatasim
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "BCM" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "GAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static PHOSshuttle-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "LED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "PED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "CPVPED" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "CPVBCM" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES PHOSbase-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("PHOS" "CPVGAIN" "${STATIC_DEPENDENCIES}")
