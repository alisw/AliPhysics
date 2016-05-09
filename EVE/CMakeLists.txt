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
add_subdirectory(converters/ConverterCollada)
if(ROOT_HTTP_FOUND)
add_subdirectory(converters/ConverterJSON)
endif(ROOT_HTTP_FOUND)

# EVE executables
add_subdirectory(executables)

# Installing macros and examples
install(DIRECTORY macros DESTINATION EVE)
install(DIRECTORY resources DESTINATION EVE)
install(DIRECTORY converters/ConverterCollada/examples DESTINATION EVE/ConverterCollada)

message(STATUS "EVE enabled")

