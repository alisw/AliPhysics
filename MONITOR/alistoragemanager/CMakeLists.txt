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
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/alistoragemanager
                   )

# Additional include folders in alphabetical order
include_directories(${AliRoot_SOURCE_DIR}/MONITOR/MONITOR
                    ${AliRoot_SOURCE_DIR}/MONITOR/MONITORzmq
                    ${AliRoot_SOURCE_DIR}/STEER/ESD
                    ${AliRoot_SOURCE_DIR}/STEER/STEERBase
                   )

add_definitions(-DZMQ)
include_directories(${ZEROMQ_INCLUDE_DIR})

# alistorage executable
add_executable(alistorage alistorage.cxx)
target_link_libraries(alistorage MONITOR MONITORzmq)
if(PYTHIA6)
    target_link_libraries(alistorage EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alistorage RUNTIME DESTINATION bin)

# Header files: set
set(HDRS
    AliStorageTypes.h
   )

# Header files: install
install(FILES ${HDRS} DESTINATION include)

# alistorageadmin executable
add_executable(alistorageadmin alistorageadmin.cxx)
target_link_libraries(alistorageadmin MONITOR MONITORzmq)
if(PYTHIA6)
    target_link_libraries(alistorageadmin EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alistorageadmin RUNTIME DESTINATION bin)


# alifakedisplay executable
add_executable(alifakedisplay alifakedisplay.cxx)
target_link_libraries(alifakedisplay MONITOR MONITORzmq)
#target_link_libraries(alifakedisplay -lzmq)
if(PYTHIA6)
    target_link_libraries(alifakedisplay EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alifakedisplay RUNTIME DESTINATION bin)


# alistorageminimalexample executable
add_executable(alistorageminimalexample alistorageminimalexample.cxx)
target_link_libraries(alistorageminimalexample MONITOR MONITORzmq)
if(PYTHIA6)
    target_link_libraries(alistorageminimalexample EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alistorageminimalexample RUNTIME DESTINATION bin)


# alistoragegetpoints executable
add_executable(alistoragegetpoints alistoragegetpoints.cxx)
target_link_libraries(alistoragegetpoints MONITOR MONITORzmq XMLIO)
if(PYTHIA6)
    target_link_libraries(alistoragegetpoints EVGEN HepMC HepMCParser)
endif(PYTHIA6)

# Installation
install(TARGETS alistoragegetpoints RUNTIME DESTINATION bin)

# alistorage2 executable
add_executable(alistorage2 alistorage2.cxx)
target_link_libraries(alistorage2 MONITOR MONITORzmq)
#target_link_libraries(alistorage2 -lzmq)

# Installation
install(TARGETS alistorage2 RUNTIME DESTINATION bin)


