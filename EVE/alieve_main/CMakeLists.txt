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
include_directories(${AliRoot_SOURCE_DIR}/EVE/alieve_main)

# Additional include folders in alphabetical order
include_directories(
                    ${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/EVE/EveBase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# alieve executable
# Sources in alphabetical order
set(SRCS
    alieve_main.cxx
   )

# Headers from sources
string(REPLACE ".cxx" ".h" HDRS "${SRCS}")

add_executable(alieve ${SRCS})
target_link_libraries(alieve EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC STEER ESD CDB AOD RAWDatabase RAWDatarec Thread ANALYSIS HLTbase GenVector GeomPainter Proof Gui XMLParser ANALYSISalice ITSbase ITSrec TPCbase TPCrec TPCsim Ged RGL AliHLTHOMER ANALYSISaliceBase XMLIO RAWDatasim HistPainter ITSsim Spectrum ${LIB_FTGL} GLEW  STEERBase EveBase Eve Core Rint ASImage)

if(ZEROMQ_FOUND)
    target_link_libraries(alieve MONITORzmq)
    add_definitions(-DZMQ)
endif(ZEROMQ_FOUND)

# alihlteve executable
set(SRCSHLT
    alihlteve_main.cxx
)

add_executable(alihlteve ${SRCSHLT})
target_link_libraries(alihlteve EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC STEER ESD CDB AOD RAWDatabase RAWDatarec Thread ANALYSIS HLTbase GenVector GeomPainter Proof Gui XMLParser ANALYSISalice ITSbase ITSrec TPCbase TPCrec TPCsim Ged RGL AliHLTHOMER ANALYSISaliceBase XMLIO RAWDatasim HistPainter ITSsim Spectrum ${LIB_FTGL} GLEW STEERBase EveBase Eve Core Rint ASImage)

if(ZEROMQ_FOUND)
    target_link_libraries(alihlteve EveBase Eve Rint MONITORzmq)
endif(ZEROMQ_FOUND)

# Installation
install(TARGETS alieve alihlteve RUNTIME DESTINATION bin)
