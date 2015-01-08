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
include_directories(${AliRoot_SOURCE_DIR}/MUON/DA)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcalib
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/MUON/MUONraw
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )


set(STATIC_DEPENDENCIES MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatarec-static)
generateDA("MUON" "TRG" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ESD-static MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatarec-static STEERBase-static)
generateDA("MUON" "TRKGAIN" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatabase-static RAWDatarec-static STEERBase-static)
generateDA("MUON" "TRKOCC" "${STATIC_DEPENDENCIES}")

set(STATIC_DEPENDENCIES ESD-static MUONcalib-static MUONcore-static MUONmapping-static MUONraw-static RAWDatarec-static STEERBase-static)
generateDA("MUON" "TRKPED" "${STATIC_DEPENDENCIES}")