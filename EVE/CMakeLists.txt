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

# EVE libraries
add_subdirectory(EveBase)
add_subdirectory(EveDet)
add_subdirectory(EveHLT)
add_subdirectory(Converter)
add_subdirectory(ColladaConverter)
add_subdirectory(aliconvertertest)

# EVE: alieve & alihlteve executables
add_subdirectory(alieve_main)

# Installing the macros
install(DIRECTORY macros DESTINATION EVE)
install(DIRECTORY alice-macros DESTINATION EVE)
install(DIRECTORY alice-data DESTINATION EVE)
install(DIRECTORY hlt-macros DESTINATION EVE)
install(DIRECTORY icons DESTINATION EVE)
install(DIRECTORY ColladaConverter/examples DESTINATION EVE/ColladaConverter)

message(STATUS "EVE enabled")

