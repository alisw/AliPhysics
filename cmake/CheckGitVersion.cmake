# 1. Extracts versioning information from the Git repository
# 2. Enables rerun of cmake configuration on each pull: GetGitRevisionDescription
#  - ALIROOT_VERSION - branch/tag name or short hash if detached at randon hash
#  - ALIROOT_REVISION - short sha1
#  - ALIROOT_SERIAL - number of commits
#  - ALIROOT_GIT_DESCRIBE - output of "git describe --tags"
#  - ALIROOT_VERSION_RPM - output of "git describe --tags" suitable for using it
#                          as a version in a RPM spec, i.e. with all "weird"
#                          characters replaced with dashes
#  - ALIROOT_GITREFTYPE - BRANCH/TAG/DETACHED

# Setting default values
set(ALIROOT_VERSION "")
set(ALIROOT_REVISION "")
set(ALIROOT_SERIAL 0)
set(ALIROOT_GITREFTYPE "")

# Checks if the sources where cloned as a full git repository
if((EXISTS ${AliRoot_SOURCE_DIR}/.git/) OR (EXISTS ${AliRoot_SOURCE_DIR}/.git))
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

  # Get output of git describe
  execute_process(COMMAND sh -c "git describe --tags || (echo -n none-0-g ; git rev-parse --short HEAD)"
                  WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                  OUTPUT_VARIABLE ALIROOT_GIT_DESCRIBE
                  ERROR_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(CMAKEDEBUG)
    message(STATUS "DEBUG: GIT_REFSPEC = \"${GIT_REFSPEC}\", GIT_SHA1 = \"${GIT_SHA1}\"")
  endif(CMAKEDEBUG)

  # Setting ALIROOT_REVISION as the long hash
  set(ALIROOT_REVISION ${GIT_SHA1})

  # Generate the short version of the revision hash
  execute_process(COMMAND git rev-parse --short ${GIT_SHA1} 
                  WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                  OUTPUT_VARIABLE GIT_SHORT_SHA1
                  RESULT_VARIABLE process_result
                  ERROR_VARIABLE process_error
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE
                )
  # Set ALIROOT_REVISION to short hash if no error
  if(process_result EQUAL 0)
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: Short SHA1 = \"${GIT_SHORT_SHA1}\"")
    endif(CMAKEDEBUG)

    set(ALIROOT_REVISION ${GIT_SHORT_SHA1})
  else()
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: result = \"${process_result}\",  parse-rev error : ${ERROR_VARIABLE}")
    endif()

    message(WARNING "Could not retrieve the short hash, using the long version : \"${ALIROOT_REVISION}\"")
  endif()

  # Generate ALIROOT_VERSION
  # 1. Branch -> Branch name
  # 2. Detached mode
  #    2.1 Tags -> Tag name
  #    2.2 Detached hash -> Short hash
  
  # Check if dettached mode
  # rev-parse will return:
  # 1. Branch -> Branch name, ex: master
  # 2. Detached mode: "HEAD" for both tags and random hashes
  execute_process(COMMAND git rev-parse --abbrev-ref HEAD
                  WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
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
                      WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
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

        set(ALIROOT_VERSION ${tag_output})
        set(ALIROOT_GITREFTYPE "TAG")
      else()
        # Detached at a random hash, the version is the short SHA1
        if(CMAKEDEBUG)
          message(STATUS "DEBUG: git describe tar_error = ${tag_error}")
        endif()  

        set(ALIROOT_VERSION ${ALIROOT_REVISION})
        set(ALIROOT_GITREFTYPE "DETACHED")
      endif()
    else()
      # Branch
      set(ALIROOT_VERSION ${ref_output})
      set(ALIROOT_GITREFTYPE "BRANCH")
    endif()
  else(ref_result EQUAL 0)
    message(FATAL_ERROR "Could not retreive information about the current git hash: ${ref_error}")
  endif(ref_result EQUAL 0)
  
  # Generating the ALIROOT_SERIAL using git rev-list
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
                    WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                    RESULT_VARIABLE serial_result
                    ERROR_VARIABLE serial_error
                    OUTPUT_VARIABLE serial_output
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
                  )
  else()
    execute_process(COMMAND git rev-list --count ${GIT_SHA1}
                    WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                    RESULT_VARIABLE serial_result
                    ERROR_VARIABLE serial_error
                    OUTPUT_VARIABLE serial_output
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
      )

  endif()

  if(serial_result EQUAL 0)
    if(CMAKEDEBUG)
      message(STATUS "DEBUG: AliRoot serial: ${serial_output}")
    endif()
    
    set(ALIROOT_SERIAL ${serial_output})
  else()
    message(FATAL_ERROR "Could not retrieve serial number: ${serial_error}")
  endif()

  if(${ALIROOT_GITREFTYPE} STREQUAL "DETACHED")
    message(STATUS "Found AliRoot in detached mode, hash \"${ALIROOT_REVISION}\", serial \"${ALIROOT_SERIAL}\", describe \"${ALIROOT_GIT_DESCRIBE}\"")
  elseif(${ALIROOT_GITREFTYPE} STREQUAL "BRANCH")
    message(STATUS "Found AliRoot branch \"${ALIROOT_VERSION}\", hash \"${ALIROOT_REVISION}\", serial \"${ALIROOT_SERIAL}\", describe \"${ALIROOT_GIT_DESCRIBE}\"")
  elseif(${ALIROOT_GITREFTYPE} STREQUAL "TAG")
    message(STATUS "Found AliRoot tag \"${ALIROOT_VERSION}\", hash \"${ALIROOT_REVISION}\", serial \"${ALIROOT_SERIAL}\", describe \"${ALIROOT_GIT_DESCRIBE}\"")
  else()
    # it does not get here
    message(FATAL_ERROR "Git type error")
  endif()
else()
    message(WARNING "AliRoot sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif()

# Use the output of git describe for RPM versions, but strip invalid characters.
string(REGEX REPLACE "[^A-Za-z0-9_.]" "." ALIROOT_VERSION_RPM ${ALIROOT_GIT_DESCRIBE})
string(REGEX REPLACE "\\.([0-9]+)\\.(g[a-f0-9]+)$" "_\\1_\\2" ALIROOT_VERSION_RPM ${ALIROOT_VERSION_RPM})
if(CMAKEDEBUG)
  message(STATUS "DEBUG: ALIROOT_VERSION_RPM = ${ALIROOT_VERSION_RPM}")
endif()

# Generating ARVersion.h from ARVersion.h.tmp
set(ALIROOT_AR_VERSION ${ALIROOT_VERSION})
if(${ALIROOT_GITREFTYPE} STREQUAL "BRANCH")
  set(ALIROOT_AR_REVISION "")
  set(ALIROOT_AR_SERIAL 0)
else()
  set(ALIROOT_AR_REVISION ${ALIROOT_REVISION})
  set(ALIROOT_AR_SERIAL ${ALIROOT_SERIAL})
endif()
configure_file(${PROJECT_SOURCE_DIR}/cmake/ARVersion.h.tmp ${CMAKE_BINARY_DIR}/version/ARVersion.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/version/ARVersion.h DESTINATION include)
