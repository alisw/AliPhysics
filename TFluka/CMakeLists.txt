# -*- mode: cmake -*-
# Create a library called "lib<name>" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/TFluka
${CMAKE_SOURCE_DIR}/include 
${ROOT_INCLUDE_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

SetModule()

set(TFluka_SRCS
TFluka.cxx
TFlukaMCGeometry.cxx
TFlukaCerenkov.cxx
TFlukaConfigOption.cxx
TFlukaScoringOption.cxx
magfld.cxx
source.cxx
mgdraw.cxx
bxdraw.cxx
eedraw.cxx
endraw.cxx
sodraw.cxx
usdraw.cxx
stupre.cxx
stuprf.cxx
abscff.cxx
dffcff.cxx
queffc.cxx
rflctv.cxx
rfrndx.cxx
)
set(TFluka_FSRCS
FLUKA_input.f
crnkvp.f
)

If(RULE_CHECKER_FOUND)
  CHECK_RULES("${TFluka_SRCS}" "${INCLUDE_DIRECTORIES}" TFluka_RULES)
endIf(RULE_CHECKER_FOUND)

# fill list of header files from list of source files
# by exchanging the file extension
CHANGE_FILE_EXTENSION(*.cxx *.h TFluka_HEADERS "${TFluka_SRCS}")

set(TFluka_LINKDEF TFlukaLinkDef.h)
set(TFluka_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/TFlukaDict.cxx) 

ROOT_GENERATE_DICTIONARY("${TFluka_HEADERS}" "${TFluka_LINKDEF}" "${TFluka_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(TFluka_SRCS ${TFluka_SRCS} ${TFluka_DICTIONARY})

add_library(TFluka SHARED ${TFluka_SRCS}  ${TFluka_FSRCS})
target_link_libraries(TFluka ${ROOT_LIBRARIES})
#set_target_properties(TFluka PROPERTIES  ${CBMROOT_LIBRARY_PROPERTIES})
SET_TARGET_PROPERTIES(TFluka PROPERTIES LINKER_LANGUAGE CXX)

################ install ###################
install(TARGETS TFluka DESTINATION ${CMAKE_BINARY_DIR}/lib)

