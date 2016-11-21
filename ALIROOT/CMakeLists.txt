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
include_directories(${AliRoot_SOURCE_DIR}/ALIROOT)

# Additional include folders in alphabetical order
include_directories(
                    ${AliRoot_BINARY_DIR}/version
                    ${AliRoot_SOURCE_DIR}/STEER/STEER
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

# aliroot executable
# Sources in alphabetical order
set(SRCS
  aliroot.cxx
)
  
set(CSRCS
  fpu.c
)

add_executable(aliroot ${SRCS} ${CSRCS})
target_link_libraries(aliroot MUONevaluation MUONmapping MUONshuttle MUONgraphics MUONsim MUONrec MUONgeometry MUONcalib MUONbase MUONraw MUONtrigger MUONcore TPCbase TPCsim TPCrec TPCutil  ITSbase ITSsim ITSrec PMDbase PMDsim PMDrec TRDbase TRDsim TRDrec FMDbase FMDsim FMDrec TOFbase TOFrec TOFsim PHOSUtils PHOSbase PHOSsim PHOSrec ADbase ADsim ADrec ACORDEbase ACORDEsim ACORDErec HMPIDbase HMPIDrec HMPIDsim ZDCbase ZDCsim ZDCrec VZERObase VZEROsim VZEROrec MFTbase MFTsim MFTrec EMCALUtils EMCALbase EMCALsim EMCALrec EMCALraw BCM STRUCT T0base T0sim T0rec FASTSIM HLTbase HLTshuttle TRIGGERbase STEER STAT CDB AOD  STEERBase ESD ANALYSIS RAWDatasim RAWDatarec RAWDatabase Rint MLP Gui Physics Geom EG Hist MathCore VMC Matrix Minuit Gpad Graf GenVector Proof Spectrum ANALYSISalice Graf3d HistPainter XMLParser AliHLTHOMER GeomPainter Thread ANALYSISaliceBase XMLIO)
if(TARGET microcern)
  target_link_libraries(aliroot microcern)
endif()

# Installation
install(TARGETS aliroot RUNTIME DESTINATION bin)

# The VDT wrapper
include_directories(${ROOT_INCLUDE_DIR})
add_library(alivdtwrapper SHARED alivdtwrapper.cxx)
target_link_libraries(alivdtwrapper dl m)
install(TARGETS alivdtwrapper DESTINATION lib)
