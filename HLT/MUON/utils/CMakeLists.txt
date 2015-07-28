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
include_directories(${AliRoot_SOURCE_DIR}/HLT/MUON/utils)

# dHLTdumpraw executable
set(SRCS1
    dHLTdumpraw.cxx
   )

add_executable(dHLTdumpraw ${SRCS1})

target_link_libraries(dHLTdumpraw AliHLTMUON HLTbase MUONmapping MUONevaluation MUONsim MUONrec MUONgeometry MUONcalib MUONbase MUONraw RAWDatasim MUONtrigger MUONcore STEER CDB AOD RAWDatarec ESD RAWDatabase STEERBase ANALYSIS AliHLTHOMER Geom Graf Hist MathCore XMLParser Physics Gui EG VMC Matrix Minuit Gpad GenVector GeomPainter Proof Graf3d Thread)

# dHLTrootify executable
set(SRCS2
    dHLTrootify.cxx
)

add_executable(dHLTrootify ${SRCS2})
target_link_libraries(dHLTrootify AliHLTMUON HLTbase MUONmapping MUONevaluation MUONsim MUONrec MUONgeometry MUONcalib MUONbase MUONraw RAWDatasim MUONtrigger MUONcore STEER CDB AOD RAWDatarec ESD RAWDatabase STEERBase ANALYSIS AliHLTHOMER Geom Graf Hist MathCore XMLParser Physics Gui EG VMC Matrix Minuit Gpad GenVector GeomPainter Proof Graf3d Thread)

# Installation
install(TARGETS dHLTdumpraw dHLTrootify RUNTIME DESTINATION bin)
