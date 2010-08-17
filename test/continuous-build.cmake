#****Macros****
find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

#****Configuration Variables****
#Set Project Name
set(CTEST_PROJECT_NAME "AliRoot")
set(ALICE_ROOT $ENV{ALICE_ROOT})
set(ALICE_INSTALL $ENV{ALICE_INSTALL})

if(NOT ALICE_INSTALL)
message(FATAL_ERROR "Please set environment variable ALICE_INSTALL to the AliRoot installation directory")
endif(NOT ALICE_INSTALL)
if(NOT ALICE_ROOT)
message(FATAL_ERROR "Please set environment variable ALICE_ROOT to the AliRoot source directory")
endif(NOT ALICE_ROOT)

#Set hostname of current machine
find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
set(CTEST_SITE "${HOSTNAME}")

#Set build name
getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)
set(CTEST_BUILD_NAME        "${osname}-${cpu}-prod")

#Detect SVN and Make
find_program(CTEST_SVN_COMMAND NAMES svn)
find_program(MAKE NAMES make)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-limit=yes --leak-resolution=low")

#AliRoot Repository and build settings
set(ALICE_REPO "https://alisoft.cern.ch/AliRoot/trunk")
set(MODEL "Continuous")


#Set Nightly build start time
set(CTEST_NIGHTLY_START_TIME "00:00:00 EST")
set(CTEST_TEST_TIMEOUT 0)
#Default drop site
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "alirootbuildcmake.cern.ch")
set(CTEST_DROP_LOCATION "/~johny/submit.php?project=AliRoot")
set(CTEST_DROP_SITE_CDASH TRUE)

set(OLDREV 0)
set(NEWREV 0)
set(CYCLE 0)
set(CLEARCYCLE 5)
set(TESTCYCLE 2)
#Get Revision
execute_process(COMMAND svn info ${ALICE_ROOT} OUTPUT_VARIABLE _svn_out )
string(REGEX REPLACE "^.*Revision: ([^\n]*).*$" "\\1" ALIROOT_SVN_REVISION ${_svn_out})
set(OLDREV ${ALIROOT_SVN_REVISION})       
set(NEWREV ${ALIROOT_SVN_REVISION})       
message("Starting build of revision ${OLDREV}")
#Setup Directories
set(CTEST_BINARY_DIRECTORY  "${ALICE_INSTALL}")
set(CTEST_SOURCE_DIRECTORY  "${ALICE_ROOT}")
set(CTEST_DASHBOARD_ROOT    "$ENV{HOME}/Dashboards")
#  set(CTEST_CHECKOUT_COMMAND "${CTEST_SVN_COMMAND} co ${ALICE_REPO} $ENV{ALICE_ROOT}")
set(CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")
set(CTEST_CMAKE_COMMAND "\"${CMAKE_EXECUTABLE_NAME}\" -D Continuous")

#CMake Generator
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

message("Continuous Build Starting")

#****Testing Process****
while (${CTEST_ELAPSED_TIME} LESS 36000)
  
  set (START_TIME ${CTEST_ELAPSED_TIME})    
  ctest_start (Continuous)
      
  #Build only if new revision of source is available
  while(${OLDREV} EQUAL ${NEWREV} AND ${CYCLE} GREATER 0)
    message("Updating source")
    message("Current revision : ${OLDREV}")
    ctest_update(BUILD "${CTEST_SOURCE_DIRECTORY}")
    execute_process(COMMAND svn info ${ALICE_ROOT} OUTPUT_VARIABLE _svn_out)
    string(REGEX REPLACE "^.*Revision: ([^\n]*).*$" "\\1" ALIROOT_SVN_REVISION ${_svn_out})
    set(NEWREV ${ALIROOT_SVN_REVISION})             
    message("New revision : ${NEWREV}")
    if(${OLDREV} EQUAL ${NEWREV})
      message("Source not updated, no need to build source")
      ctest_sleep(300)
    endif(${OLDREV} EQUAL ${NEWREV})
  endwhile()
  #    ctest_submit (PARTS Update)
  set(OLDREV ${NEWREV})
  math(EXPR CYCLE "${CYCLE} + 1")
  message("Build #${CYCLE} ")
  message("Revision : ${OLDREV}")
  set (START_TIME ${CTEST_ELAPSED_TIME})    
  math(EXPR CLEAR "${CYCLE}%${CLEARCYCLE}")
  
  #Clear build directory every CLEARCYCLE build
  if(${CLEAR} EQUAL 0)
    message("Clearing build directory")
    ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
  endif(${CLEAR} EQUAL 0) 
  message(STATUS "Building source")
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}")
  
  #Test every other build
  math(EXPR TEST "${CYCLE}%${TESTCYCLE}")
  if(${TEST} EQUAL 0)
    message(STATUS "Running tests")
    ctest_test (BUILD "${CTEST_BINARY_DIRECTORY}")
  endif(${TEST} EQUAL 0)
  ctest_submit ()
  message("Build ${CYCLE} completed")
  ctest_sleep(300)
endwhile()




