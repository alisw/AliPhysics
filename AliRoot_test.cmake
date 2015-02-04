 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
SET (CTEST_SOURCE_DIRECTORY $ENV{SOURCEDIR})
SET (CTEST_BINARY_DIRECTORY $ENV{BUILDDIR})
SET (CTEST_ROOTSYS_DIRECTORY $ENV{ROOTSYS})
SET (CTEST_SITE $ENV{SITE})
SET (CTEST_BUILD_NAME $ENV{LABEL})
SET (CTEST_CMAKE_GENERATOR "Unix Makefiles")
SET (CTEST_PROJECT_NAME "ALIROOT")

Find_program(CTEST_GIT_COMMAND NAMES git)
Set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

#If($ENV{ctest_model} MATCHES Continuous)
#  Set(CTEST_SVN_UPDATE_OPTIONS "$ENV{REVISION}")
#EndIf($ENV{ctest_model} MATCHES Continuous)


Set (CTEST_CONFIGURE_COMMAND " \"${CMAKE_EXECUTABLE_NAME}\"   \"-DROOTSYS=${CTEST_ROOTSYS_DIRECTORY}\"  \"${CTEST_SOURCE_DIRECTORY}\" ")

SET (BUILD_COMMAND "make")
SET (CTEST_BUILD_COMMAND "${BUILD_COMMAND} -j$ENV{number_of_processors}")

If($ENV{ctest_model} MATCHES Nightly)

  Set (CTEST_CONFIGURE_COMMAND " \"${CMAKE_EXECUTABLE_NAME}\" \"-DCMAKE_BUILD_TYPE=NIGHTLY\" \"-G${CTEST_CMAKE_GENERATOR}\" \"${CTEST_SOURCE_DIRECTORY}\" \"-DROOTSYS=${CTEST_ROOTSYS_DIRECTORY}\"  ")
  CTEST_EMPTY_BINARY_DIRECTORY(${CTEST_BINARY_DIRECTORY})

endif($ENV{ctest_model} MATCHES Nightly)


configure_file(${CTEST_SOURCE_DIRECTORY}/CTestCustom.cmake
               ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake
              )
ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")

CTEST_START ($ENV{ctest_model})
CTEST_UPDATE (SOURCE "${CTEST_SOURCE_DIRECTORY}")
CTEST_CONFIGURE (BUILD "${CTEST_BINARY_DIRECTORY}")
CTEST_BUILD (BUILD "${CTEST_BINARY_DIRECTORY}")
CTEST_TEST (BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL $ENV{number_of_processors})
CTEST_SUBMIT ()
 
