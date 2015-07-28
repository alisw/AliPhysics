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

# alieve executable

# Module include folder
include_directories(${AliRoot_SOURCE_DIR}/MUON/mchview)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_SOURCE_DIR}/MUON/MUONcore
                    ${AliRoot_SOURCE_DIR}/MUON/MUONgraphics
                    ${AliRoot_SOURCE_DIR}/MUON/MUONmapping
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# aliroot executable
# Sources in alphabetical order
set(SRCS
  mchview.cxx
)

add_executable(mchview ${SRCS})
target_link_libraries(mchview MUONmapping RAWDatabase MUONshuttle MUONgraphics MUONrec MUONgeometry MUONcalib MUONbase MUONraw RAWDatasim MUONtrigger MUONcore STEER CDB AOD RAWDatarec ESD RAWDatabase STEERBase ANALYSIS Physics Hist Geom MathCore Matrix Minuit VMC Gpad Graf HLTbase EG GenVector GeomPainter Proof Gui XMLParser Graf3d Thread Rint)

# Installation
install(TARGETS mchview RUNTIME DESTINATION bin)
