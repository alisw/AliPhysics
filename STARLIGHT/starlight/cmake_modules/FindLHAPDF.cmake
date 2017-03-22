#Will look for LHAPDF libraries in /usr/lib/ or /usr/local/lib/. If LHAPDF is not found, then look for dummy lhapdf in pythia8.

FIND_PATH(LHAPDF_INCLUDE_DIR Pythia.h /usr/include/ /usr/local/include/ $ENV{PYTHIADIR}/include)

FIND_LIBRARY(LHAPDF_LIBRARIES NAMES lhapdfdummy PATHS /usr/lib /usr/local/lib $ENV{PYTHIADIR}/lib/archive/)

IF (LHAPDF_INCLUDE_DIR AND LHAPDF_LIBRARIES)
   SET(LHAPDF_FOUND TRUE)
ENDIF (LHAPDF_INCLUDE_DIR AND LHAPDF_LIBRARIES)

IF (LHAPDF_FOUND)
   IF (NOT Pythia8_FIND_QUIETLY)
      MESSAGE(STATUS "Found LHAPDF: ${LHAPDF_LIBRARIES}")
   ENDIF (NOT Pythia8_FIND_QUIETLY)
ELSE (LHAPDF_FOUND)
   IF (Pythia8_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find LHAPDF")
   ENDIF (Pythia8_FIND_REQUIRED)
ENDIF (LHAPDF_FOUND)

