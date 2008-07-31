# -*- mode: cmake -*-

  ###########################################
  #
  #       Usefull macros
  #
  ###########################################

  ###############################################################
  #
  # Exchange file extention of LIST from
  # FILE_EXT1 to FILE_EXT2 and assign the
  # newly created list to OUTVAR. The input
  # list LIST is not changed at all
  # Ex: CHANGE_FILE_EXTENSION(*.cxx *.h TRD_HEADERS "${TRD_SRCS}")
  #
  ################################################################

MACRO (CHANGE_FILE_EXTENSION FILE_EXT1 FILE_EXT2 OUTVAR LIST)

   SET(BLA)

   IF (${FILE_EXT1} MATCHES "^[*][.]+.*$")
     STRING(REGEX REPLACE "^[*]+([.].*)$" "\\1" FILE_EXT1_NEW ${FILE_EXT1}) 
   ENDIF  (${FILE_EXT1} MATCHES "^[*][.]+.*$")

   IF (${FILE_EXT2} MATCHES "^[*][.]+.*$")
     STRING(REGEX REPLACE "^[*]+([.].*)" "\\1" FILE_EXT2_NEW ${FILE_EXT2}) 
   ENDIF  (${FILE_EXT2} MATCHES "^[*][.]+.*$")

   foreach (_current_FILE ${LIST})

#     GET_FILENAME_COMPONENT(name_wo_path ${_current_FILE}
#                            NAME
#                            )

     STRING(REGEX REPLACE "^(.*)${FILE_EXT1_NEW}$" "\\1${FILE_EXT2_NEW}" test ${_current_FILE})
     SET (BLA ${BLA} ${test})

   endforeach (_current_FILE ${LIST})
   
   SET (${OUTVAR} ${BLA})



ENDMACRO (CHANGE_FILE_EXTENSION)

  ######################################################
  # 
  # Macro get string with a colon seperated string of
  # pathes or any other colon sperated list.
  # First the string is seperated  and the entries are
  # filled into a list.
  # Loop over the list and searches for the occurence
  # of keywords which are provided as a list.
  # If the keyword occurs this path (entry) is
  # deleted from the list. Returns the list of pathes
  # (entries) wich survives the loop.
  #
  # PATH: colon separated string of pathes or other 
  #       input entries
  # LIST_OF_KEYWORDS: list of the keywords which 
  #                   should be excluded in the output
  # OUTVAR: name of the variable which should be filled
  #         with the resulting output list
  #
  ######################################################

MACRO (CLEAN_PATH_LIST PATH LIST_OF_KEYWORDS OUTVAR)

  SET(BLA "")

  STRING(REGEX MATCHALL "[^:]+" PATH1 ${PATH})

  FOREACH(_current_PATH ${PATH1})
 
    SET(KEYWORD_FOUND FALSE)

    FOREACH(_current_KEYWORD ${LIST_OF_KEYWORDS})

      IF (${_current_PATH} MATCHES "${_current_KEYWORD}")
        SET(KEYWORD_FOUND TRUE)
      ENDIF (${_current_PATH} MATCHES "${_current_KEYWORD}")

    ENDFOREACH(_current_KEYWORD ${LIST_OF_KEYWORDS})
    
    IF (NOT KEYWORD_FOUND)
      SET(BLA ${BLA} ${_current_PATH})
    ENDIF (NOT KEYWORD_FOUND)  

  ENDFOREACH(_current_PATH ${PATH1})

  UNIQUE(${OUTVAR} "${BLA}")

ENDMACRO (CLEAN_PATH_LIST)

  ##########################################################
  #
  # The macro checks if the user wants to build the project
  # in the source directory and if so stop the execution
  # of cmake with an error message.
  #
  ##########################################################

MACRO (CHECK_OUT_OF_SOURCE_BUILD)

   STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
   IF(insource)
      FILE(REMOVE_RECURSE ${CMAKE_SOURCE_DIR}/Testing)
      FILE(REMOVE ${CMAKE_SOURCE_DIR}/DartConfiguration.tcl)
      MESSAGE(FATAL_ERROR "FAIRROOT should be installed as an out of source build, to keep the source directory clean. Please create a extra build directory and run the command 'cmake path_to_source_dir' in this newly created directory. You have also to delete the directory CMakeFiles and the file CMakeCache.txt in the source directory. Otherwise cmake will complain even if you run it from an out-of-source directory.") 
   ENDIF(insource)

ENDMACRO (CHECK_OUT_OF_SOURCE_BUILD)

MACRO(UNIQUE var_name list)

  #######################################################################
  # Make the given list have only one instance of each unique element and
  # store it in var_name.
  #######################################################################

  SET(unique_tmp "")
  FOREACH(l ${list})
    STRING(REGEX REPLACE "[+]" "\\\\+" l1 ${l})
    IF(NOT "${unique_tmp}" MATCHES "(^|;)${l1}(;|$)")
      SET(unique_tmp ${unique_tmp} ${l})
    ENDIF(NOT "${unique_tmp}" MATCHES "(^|;)${l1}(;|$)")
  ENDFOREACH(l)
  SET(${var_name} ${unique_tmp})
ENDMACRO(UNIQUE)

MACRO(CHECK_CMAKE_VERSION)

  #################################################################
  #Check if CMake has at least version 2.4.3
  # This has to be done before anything else, because some commands
  # are not known in older versions.
  #################################################################

  EXEC_PROGRAM( cmake ARGS "--version" OUTPUT_VARIABLE CMAKE_VERSION )
  IF (${CMAKE_MAJOR_VERSION} LESS 2 )
      MESSAGE("You are using CMake version ${CMAKE_VERSION} .")
      MESSAGE( FATAL_ERROR "This Cmake version is to old. At leasts version 2.4.3 is needed.")
  ENDIF (${CMAKE_MAJOR_VERSION} LESS 2 )
  IF (${CMAKE_MINOR_VERSION} LESS 4 )
      MESSAGE("You are using CMake version ${CMAKE_VERSION} .")
      MESSAGE( FATAL_ERROR "This Cmake version is to old. At leasts version 2.4.3 is needed.")
  ENDIF (${CMAKE_MINOR_VERSION} LESS 4 )
  IF (${CMAKE_PATCH_VERSION} LESS 3 )
      MESSAGE("You are using CMake version ${CMAKE_VERSION} .")
      MESSAGE( FATAL_ERROR "This Cmake version is to old. At leasts version 2.4.3 is needed.")
  ENDIF (${CMAKE_PATCH_VERSION} LESS 3 )

ENDMACRO(CHECK_CMAKE_VERSION)

MACRO(CHECK_GSI)

  #################################################################
  # Check if /misc/cbmsoft/config exist, because than we are at GSI
  # If we are not at GSI some of the macros may not work.
  #################################################################

  SET(GSI GSI-NOTFOUND)
  FIND_PATH(GSI NAMES config PATHS
    /misc/cbmsoft
    NO_DEFAULT_PATH
  ) 

  EXEC_PROGRAM(uname ARGS "-m" OUTPUT_VARIABLE Machine)
#  SET(Machine ${CMAKE_SYSTEM_PROCESSOR})
  #${CMAKE_SYSTEM_PROCESSOR} is not set correctly on 64bit
  #Debian Etch, so use workaround

  IF(GSI)
    MESSAGE("-- You're using the GSI installation of the external packages.")

    SET(KEYWORDS
    cbmsoft
    debian
    globus
    gcc
    binutils
    cbm 
    panda
    )
    IF(DEFINED ENV{PATH})
      CLEAN_PATH_LIST("$ENV{PATH}" "${KEYWORDS}" PATH)
      #MESSAGE ("PATH: ${PATH}")
      #set(ENV{PATH} ${PATH})
    ENDIF(DEFINED ENV{PATH})


    # Create a clean LD_LIBRARY_PATH without any pathes to libraries set
    # later on
    IF(DEFINED ENV{LD_LIBRARY_PATH})
      SET(KEYWORDS ${KEYWORDS} ${CMAKE_BINARY_DIR})
      CLEAN_PATH_LIST("$ENV{LD_LIBRARY_PATH}" "${KEYWORDS}" LD_LIBRARY_PATH)
      #MESSAGE("LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}")
    ENDIF(DEFINED ENV{LD_LIBRARY_PATH})

    # Check if you're on an 32 or 64 bit environment
    IF(${Machine} MATCHES "i686" OR 
       ${Machine} MATCHES "i386" )
      MESSAGE("-- This is an 32 bit machine")
      SET(GSI_SIMPATH /misc/cbmsoft/Debian3.1/new/)
    ELSE(${Machine} MATCHES "i686" OR 
         ${Machine} MATCHES "i386" )
      IF(${Machine} MATCHES "x86_64")
        MESSAGE("-- This is an 64 bit machine")
        SET(GSI_SIMPATH /misc/cbmsoft/Debian64/new/)
      ELSE(${Machine} MATCHES "x86_64")
        MESSAGE(FATAL_ERROR "This is a UNIX machine wich is neither 32 nor 64 bit. Don't know what to do, so break. ")
      ENDIF(${Machine} MATCHES "x86_64")
    ENDIF(${Machine} MATCHES "i686" OR 
          ${Machine} MATCHES "i386" )

  ELSE(GSI)
    MESSAGE("-- You're not using the GSI installation of the external packages.")
    MESSAGE("-- If any problems occur this can be due to the fact that the macros")
    MESSAGE("-- were mostly tested at GSI. If you find any problems please contact")
    MESSAGE("-- f.uhlig@gsi.de")
    STRING(REGEX MATCHALL "[^:]+" PATH $ENV{PATH})
    #MESSAGE("PATH: ${PATH}")
  ENDIF(GSI)

ENDMACRO(CHECK_GSI)

MACRO(CHECK_EXTERNAL_PACKAGES_INSTALLATION)
 
  #############################################################
  # Check if cmake can find the root of the simulation packages
  # check if the installation is as at GSI
  # This is needed later on to check for the required packages
  # Should be done in a better way that you find the packages
  # independent of the way where they are installed
  #############################################################

  SET(SIMPATH SIMPATH-NOTFOUND)
  SET(SIMPATH1 SIMPATH1-NOTFOUND)
  SET(SIMPATH2 SIMPATH2-NOTFOUND)
  SET(SIMPATH3 SIMPATH3-NOTFOUND)
  SET(SIMPATH4 SIMPATH4-NOTFOUND)

  FIND_PATH(SIMPATH1 NAMES transport PATHS
    $ENV{SIMPATH}
    ${GSI_SIMPATH}
    NO_DEFAULT_PATH
  )
  FIND_PATH(SIMPATH2 NAMES tools PATHS
    $ENV{SIMPATH}
    ${GSI_SIMPATH}
    NO_DEFAULT_PATH
  )
  FIND_PATH(SIMPATH3 NAMES cern PATHS
    $ENV{SIMPATH}
    ${GSI_SIMPATH}
    NO_DEFAULT_PATH
  )
  FIND_PATH(SIMPATH4 NAMES generators PATHS
    $ENV{SIMPATH}
    ${GSI_SIMPATH}
    NO_DEFAULT_PATH
  )
  
  IF(SIMPATH1 AND SIMPATH2 AND SIMPATH3 AND SIMPATH4)
    SET(SIMPATH ${SIMPATH1})
  ENDIF(SIMPATH1 AND SIMPATH2 AND SIMPATH3 AND SIMPATH4)
  
      
  IF (${SIMPATH} MATCHES "SIMPATH-NOTFOUND")
    MESSAGE( FATAL_ERROR "Could not find the root of the simulation software. Please define SIMPATH as the path to the root of your instalation of simulation software.")
  ELSE (${SIMPATH} MATCHES "SIMPATH-NOTFOUND")
    MESSAGE( STATUS "Found root of the simulation software at ${SIMPATH}")
    SET( ENV{SIMPATH} ${SIMPATH})
    SET(ROOTSYS ${SIMPATH}/tools/root)
    SET( ENV{ROOTSYS} ${SIMPATH}/tools/root)
  ENDIF (${SIMPATH} MATCHES "SIMPATH-NOTFOUND")

ENDMACRO(CHECK_EXTERNAL_PACKAGES_INSTALLATION)

###################################################
# Creates a variable which stores the intersection 
# between two lists
####################################################

MACRO(INTERSECTION var_name list1 list2)
  # Store the intersection between the two given lists in var_name.
  SET(intersect_tmp "")
  FOREACH(l ${list1})
    IF("${list2}" MATCHES "(^|;)${l}(;|$)")
      SET(intersect_tmp ${intersect_tmp} ${l})
    ENDIF("${list2}" MATCHES "(^|;)${l}(;|$)")
  ENDFOREACH(l)
  SET(${var_name} ${intersect_tmp})
ENDMACRO(INTERSECTION)

MACRO(REMOVE_FROM_LIST var_name list1 list2)
  # Remove elements in list2 from list1 and store the result in var_name.
  SET(filter_tmp "")
  FOREACH(l ${list1})
    IF(NOT "${list2}" MATCHES "(^|;)${l}(;|$)")
      SET(filter_tmp ${filter_tmp} ${l})
    ENDIF(NOT "${list2}" MATCHES "(^|;)${l}(;|$)")
  ENDFOREACH(l)
  SET(${var_name} ${filter_tmp})
ENDMACRO(REMOVE_FROM_LIST)





