
# 1. Extracts versioning information from the Git repository
# 2. Enables rerun of cmake configuration on each pull: GetGitRevisionDescription
#  - ALIPHYSICS_VERSION - branch/tag name or short hash if detached at randon hash
#  - ALIPHYSICS_REVISION - short sha1
#  - ALIPHYSICS_SERIAL - number of commits
#  - ALIPHYSISC_VERSION_RPM - name of the branch/tag in rpm format, - replaced with .
#  - ALIPHYSICS_GITREFTYPE - BRANCH/TAG/DETACHED

# Setting default values
set(ALIPHYSICS_VERSION "")
set(ALIPHYSICS_REVISION "")
set(ALIPHYSICS_SERIAL 0)
set(ALIPHYSICS_GITREFTYPE "")

# Checks if the sources where cloned as a full git repository
if((EXISTS ${AliPhysics_SOURCE_DIR}/.git/) OR (EXISTS ${AliPhysics_SOURCE_DIR}/.git))
  # Git installation mandatory
  find_package(Git REQUIRED)

  # The simple include will not trigger the reconfiguration
  # get_git_head_revision has to be called at least once
  include(GetGitRevisionDescription)
  # GIT_SHA1 - current long hash
  # GIT_REFSPEC
  #     1. branches: refs/heads/master
  #     2. detached mode(tags or hashes) empty string
  get_git_head_revision(GIT_REFSPEC GIT_SHA1)

  if(CMAKEDEBUG)
    message(STATUS "DEBUG: GIT_REFSPEC = \"${GIT_REFSPEC}\", GIT_SHA1 = \"${GIT_SHA1}\"")
  endif(CMAKEDEBUG)

  # Setting ALIPHYSICS_REVISION as the long hash
  set(ALIPHYSICS_REVISION ${GIT_SHA1})

  # Generate the short version of the revision hash
  execute_process(COMMAND git rev-parse --short ${GIT_SHA1} 
                  WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                  OUTPUT_VARIABLE GIT_SHORT_SHA1
                  RESULT_VARIABLE process_result
                  ERROR_VARIABLE process_error
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE
                )
  # Set ALIPHYSICS_REVISION to short hash if no error
  if(process_result EQUAL 0)
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: Short SHA1 = \"${GIT_SHORT_SHA1}\"")
    endif(CMAKEDEBUG)

    set(ALIPHYSICS_REVISION ${GIT_SHORT_SHA1})
  else()
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: result = \"${process_result}\",  parse-rev error : ${ERROR_VARIABLE}")
    endif()

    message(WARNING "Could not retrieve the short hash, using the long version : \"${ALIPHYSICS_REVISION}\"")
  endif()

  # Generate ALIPHYSICS_VERSION
  # 1. Branch -> Branch name
  # 2. Detached mode
  #    2.1 Tags -> Tag name
  #    2.2 Detached hash -> Short hash
  
  # Check if dettached mode
  # rev-parse will return:
  # 1. Branch -> Branch name, ex: master
  # 2. Detached mode: "HEAD" for both tags and random hashes
  execute_process(COMMAND git rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                  OUTPUT_VARIABLE ref_output
                  RESULT_VARIABLE ref_result
                  ERROR_VARIABLE ref_error
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE
                )

  if(ref_result EQUAL 0)
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: rev-parse HEAD result = \"${ref_output}\"")
    endif()

    # detached mode
    if(ref_output STREQUAL "HEAD")
      # Checking if this is a tag in detached mode
      #  1. If tag the OUTPUT_VARIABLE will contain the tag name
      #  2. If random hash the RESULT_VARIABLE is 128 and ERROR_VARIABLE contains the error message
      execute_process(COMMAND git describe --exact-match
                      WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                      OUTPUT_VARIABLE tag_output
                      RESULT_VARIABLE tag_result
                      ERROR_VARIABLE tag_error
                      OUTPUT_STRIP_TRAILING_WHITESPACE
                      ERROR_STRIP_TRAILING_WHITESPACE
                    )

      if(tag_result EQUAL 0)
      
        if(CMAKEDEBUG)
          message(STATUS "DEBUG: git describe tag_result = ${tag_output}")
        endif()

        set(ALIPHYSICS_VERSION ${tag_output})
        set(ALIPHYSICS_GITREFTYPE "TAG")
      else()
        # Detached at a random hash, the version is the short SHA1
        if(CMAKEDEBUG)
          message(STATUS "DEBUG: git describe tar_error = ${tag_error}")
        endif()  

        set(ALIPHYSICS_VERSION ${ALIPHYSICS_REVISION})
        set(ALIPHYSICS_GITREFTYPE "DETACHED")
      endif()
    else()
      # Branch
      set(ALIPHYSICS_VERSION ${ref_output})
      set(ALIPHYSICS_GITREFTYPE "BRANCH")
    endif()
  else(ref_result EQUAL 0)
    message(FATAL_ERROR "Could not retreive information about the current git hash: ${ref_error}")
  endif(ref_result EQUAL 0)
  
  # Generating the ALIPHYSICS_SERIAL using git rev-list
  # Older Git version < 1.7.3 do not have --count option for rev-list
  # We use simple rev-list and count the lines of the output 

  if(CMAKEDEBUG)
    message(STATUS "DEBUG: GIT_VERSION_STRING = ${GIT_VERSION_STRING}")
  endif()
  
  # extract major minor and patch from Git version
  string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" GIT_VERSION_MAJOR "${GIT_VERSION_STRING}")
  string(REGEX REPLACE "^[0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" GIT_VERSION_MINOR "${GIT_VERSION_STRING}")
  string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" GIT_VERSION_PATCH "${GIT_VERSION_STRING}")

  if(CMAKEDEBUG)
    message(STATUS "DEBUG: GIT_VERSION_MAJOR = ${GIT_VERSION_MAJOR}, GIT_VERSION_MINOR = ${GIT_VERSION_MINOR}, GIT_VERSION_PATCH = ${GIT_VERSION_PATCH}")
  endif()

  if((${GIT_VERSION_MAJOR} EQUAL 1) AND (${GIT_VERSION_MINOR} LESS 8) AND (${GIT_VERSION_PATCH} LESS 3))
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: cmake version less that 1.7.3!")
    endif()
    
    execute_process(COMMAND git rev-list ${GIT_SHA1}
                    COMMAND wc -l
                    WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                    RESULT_VARIABLE serial_result
                    ERROR_VARIABLE serial_error
                    OUTPUT_VARIABLE serial_output
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
                  )
  else()
    execute_process(COMMAND git rev-list --count ${GIT_SHA1}
                    WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                    RESULT_VARIABLE serial_result
                    ERROR_VARIABLE serial_error
                    OUTPUT_VARIABLE serial_output
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
      )

  endif()

  if(serial_result EQUAL 0)
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: AliPhysics serial: ${serial_output}")
    endif()
    
    set(ALIPHYSICS_SERIAL ${serial_output})
  else()
    message(FATAL_ERROR "Could not retrieve serial number: ${serial_error}")
  endif()

  if(${ALIPHYSICS_GITREFTYPE} STREQUAL "DETACHED")
    message(STATUS "Found AliPhysics in detached mode, hash \"${ALIPHYSICS_REVISION}\", serial \"${ALIPHYSICS_SERIAL}\"")
  elseif(${ALIPHYSICS_GITREFTYPE} STREQUAL "BRANCH")
    message(STATUS "Found AliPhysics branch \"${ALIPHYSICS_VERSION}\", hash \"${ALIPHYSICS_REVISION}\", serial \"${ALIPHYSICS_SERIAL}\"")
  elseif(${ALIPHYSICS_GITREFTYPE} STREQUAL "TAG")
    message(STATUS "Found AliPhysics tag \"${ALIPHYSICS_VERSION}\", hash \"${ALIPHYSICS_REVISION}\", serial \"${ALIPHYSICS_SERIAL}\"")
  else()
    # it does not get here
    message(FATAL_ERROR "Git type error")
  endif()
else()
    message(WARNING "AliPhysics sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif()

# ALIPHYSICS_VERSION_RPM
# Replacing -/ with . , normally it should not contain / 
# - and / are forbidden characters in rpm creation
string(REPLACE "-" "." ALIPHYSICS_VERSION_RPM ${ALIPHYSICS_VERSION})
string(REPLACE "/" "." ALIPHYSICS_VERSION_RPM ${ALIPHYSICS_VERSION_RPM})
if(CMAKEDEBUG)
  message(STATUS "DEBUG: ALIPHYSICS_VERSION_RPM = ${ALIPHYSICS_VERSION_RPM}")
endif()

# Generating APVersion.h from APVersion.h.tmp
set(ALIPHYSICS_AR_VERSION ${ALIPHYSICS_VERSION})
if(${ALIPHYSICS_GITREFTYPE} STREQUAL "BRANCH")
  set(ALIPHYSICS_AR_REVISION "")
  set(ALIPHYSICS_AR_SERIAL 0)
else()
  set(ALIPHYSICS_AR_REVISION ${ALIPHYSICS_REVISION})
  set(ALIPHYSICS_AR_SERIAL ${ALIPHYSICS_SERIAL})
endif()
configure_file(${PROJECT_SOURCE_DIR}/cmake/APVersion.h.in ${CMAKE_BINARY_DIR}/version/APVersion.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/version/APVersion.h DESTINATION include)
