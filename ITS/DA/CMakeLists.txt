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
include_directories(${AliRoot_SOURCE_DIR}/HMPID/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/ITS/ITSbase
                    ${AliRoot_SOURCE_DIR}/ITS/ITSrec
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )



set(STATIC_DEPENDENCIES ITSbase ITSrec RAWDatarec)
generateDA("ITS" "SDDBAS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static)
generateDA("ITS" "SDDGAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static)
generateDA("ITS" "SDDINJ" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("ITS" "SPDFO" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("ITS" "SPDPHYS" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static ITSbase-static RAWDatarec-static STEERBase-static)
generateDA("ITS" "SPDSCAN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSrec-static RAWDatarec-static STEERBase-static CDB-static)
generateDA("ITS" "SPDVertexDiamond" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ITSbase-static)
generateDA("ITS" "SSDPED" "${STATIC_DEPENDENCIES}")

