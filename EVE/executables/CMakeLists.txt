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

# EVE executables

# Module include folder
include_directories(
                    ${AliRoot_SOURCE_DIR}/EVE/alieve_main
                   )

# Additional include folders in alphabetical order
include_directories(
                    ${ROOT_INCLUDE_DIR}
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSIS
                    ${AliRoot_SOURCE_DIR}/ANALYSIS/ANALYSISalice
                    ${AliRoot_SOURCE_DIR}/EVE/EveBase
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatarec
                    ${AliRoot_SOURCE_DIR}/RAW/RAWDatabase
                    ${AliRoot_SOURCE_DIR}/STEER/AOD
                    ${AliRoot_SOURCE_DIR}/STEER/CDB
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# aliceToJson executable
if(ROOT_HTTP_FOUND)
    include_directories(
        ${AliRoot_SOURCE_DIR}/EVE/aliceToJson
        ${AliRoot_SOURCE_DIR}/EVE/converters/ConverterJSON
        )
    add_executable(aliceToJson aliceToJson.cxx)
    target_link_libraries(aliceToJson ConverterJSON STEERBase STEER ESD)
endif(ROOT_HTTP_FOUND)

# alieve executable
add_executable(alieve alieve_main.cxx)
target_link_libraries(alieve EG Geom Gpad Graf3d Graf Hist MathCore Matrix Minuit Net Physics RIO Tree VMC STEER ESD CDB AOD RAWDatabase RAWDatarec Thread ANALYSIS HLTbase GenVector GeomPainter Proof Gui XMLParser ANALYSISalice ITSbase ITSrec TPCbase TPCrec TPCsim Ged RGL AliHLTHOMER ANALYSISaliceBase XMLIO RAWDatasim HistPainter ITSsim Spectrum ${LIB_FTGL} GLEW  STEERBase EveBase Eve Core Rint ASImage)
if(ZEROMQ_FOUND)
    target_link_libraries(alieve MONITORzmq)
    add_definitions(-DZMQ)
endif(ZEROMQ_FOUND)

# alimyo executable
#add_executable(alimyo hello-myo.cxx)
#set(MYO_LIBRARY ${AliRoot_SOURCE_DIR}/EVE/resources/myo/myo.framework)
#target_link_libraries(alimyo ${MYO_LIBRARY})
#SET_TARGET_PROPERTIES(alimyo PROPERTIES INSTALL_RPATH "${AliRoot_SOURCE_DIR}/EVE/resources/myo/")


# Installation
install(TARGETS alieve RUNTIME DESTINATION bin)

if(ROOT_HTTP_FOUND)
    install(TARGETS aliceToJson RUNTIME DESTINATION bin)
endif(ROOT_HTTP_FOUND)