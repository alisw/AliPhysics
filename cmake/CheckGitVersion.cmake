# Configure ARVerion.h using Git informatiion
# Sets 4 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIPHYSICS_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
#
#  - ALIPHYSICS_VERSION - name of the branch/tag
#  - ALIPHYSICS_VERSION_RPM - name of the branch/tag in rpm format, - replaced with .
#  - ALIPHYSICS_REVISION - short sha1
if(EXISTS ${PROJECT_SOURCE_DIR}/.git/)
    include(GetGitRevisionDescription)
    
    find_package(Git)
    
    if(GIT_FOUND)
        message(STATUS "Git version = ${GIT_VERSION_STRING}")
        
        get_git_head_revision(GIT_REFSPEC GIT_SHA1)

        # generate the short version of the revision hash
        execute_process(COMMAND git rev-parse --short ${GIT_SHA1}
                          WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                          OUTPUT_STRIP_TRAILING_WHITESPACE
                          RESULT_VARIABLE res
                          OUTPUT_VARIABLE GIT_SHORT_SHA1)

        # if the rev-parse fails we set the short sha to the long initial one
        if(NOT res EQUAL 0)
            set(GIT_SHORT_SHA1 ${GIT_SHA1})
        endif()

        # Older Git version < 1.7.3 do not have --count option for rev-list
        # We use simple rev-list and we count the lines of the output
        string(COMPARE GREATER "${GIT_VERSION_STRING}" "1.7.3" NEWGIT)
        
        if(NEWGIT)
            # generate the short version of the revision hash using --count
            execute_process(COMMAND git rev-list --count ${GIT_SHA1}
                            WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            RESULT_VARIABLE revcount
                            OUTPUT_VARIABLE ALIPHYSICS_SERIAL_ORIGINAL)
        else()
            # generate the short version of the revision hash using -wc -l
            execute_process(COMMAND git rev-list ${GIT_SHA1}
                            COMMAND wc -l
                            WORKING_DIRECTORY ${AliPhysics_SOURCE_DIR}
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            RESULT_VARIABLE revcount
                            OUTPUT_VARIABLE ALIPHYSICS_SERIAL_ORIGINAL)
        endif()

        # GIT_REFSPEC is empty for detached mode = tags in detached mode or checkout to specific hash

        # returns the closest reference to the current hash
        # name of the current tag or heads/branch in the case of branches
        # Older git version of Git report only the name of the tag
        # Newer version report tags/vAN-20141215
        # Just in case we replace tags/ with nothing
        git_describe(ALIPHYSICS_GIT_TAG "--all" "--abbrev=0")
        
        if(ALIPHYSICS_GIT_TAG)
            string(STRIP ${ALIPHYSICS_GIT_TAG} ALIPHYSICS_GIT_TAG)
            string(REPLACE "tags/" ""  ALIPHYSICS_GIT_TAG ${ALIPHYSICS_GIT_TAG})
        endif(ALIPHYSICS_GIT_TAG)
        
        # using the closest tag for branches
        git_describe(ALIPHYSICS_CLOSEST_GIT_TAG "--abbrev=0")

        if(ALIPHYSICS_GIT_TAG)
            string(STRIP ${ALIPHYSICS_GIT_TAG} ALIPHYSICS_GIT_TAG)
            string(REPLACE "tags/" ""  ALIPHYSICS_GIT_TAG ${ALIPHYSICS_GIT_TAG})
        endif(ALIPHYSICS_GIT_TAG)
        
        STRING(REGEX REPLACE "^(.+/)(.+)/(.*)$" "\\2" BRANCH_TYPE "${GIT_REFSPEC}" )
        
        # the revision is not set in the case of a branch, it means we are doing development
        # and the revision will trigger a reconfiguration
        if(BRANCH_TYPE STREQUAL "heads")
            set(ALIPHYSICS_REVISION "ThisIsaBranchNoRevisionProvided")
            set(ALIPHYSICS_SERIAL 0)
            set(ALIPHYSICS_GIT_TAG ${ALIPHYSICS_CLOSEST_GIT_TAG})
            STRING(REGEX REPLACE "^(.+/)(.+/)(.*)$" "\\3" SHORT_BRANCH "${GIT_REFSPEC}" )
            message(STATUS "This is a working branch, ARVersion will not contain the revision and the serial number")
        else()
            set(BRANCH_TYPE "tags")
            set(SHORT_BRANCH ${ALIPHYSICS_GIT_TAG})
            set(ALIPHYSICS_REVISION ${GIT_SHORT_SHA1})
            set(ALIPHYSICS_SERIAL ${ALIPHYSICS_SERIAL_ORIGINAL})
        endif()

        set(ALIPHYSICS_BRANCH ${SHORT_BRANCH})
        set(ALIPHYSICS_VERSION ${SHORT_BRANCH})
        
        message(STATUS "ALIPHYSICS branch/tag: \"${ALIPHYSICS_VERSION}\" - Revision:  \"${GIT_SHORT_SHA1}\" - Serial: \"${ALIPHYSICS_SERIAL_ORIGINAL}\"")

    else()
        message(STATUS "Git not installed. I can't tell you which revision you are using!")
    endif(GIT_FOUND)
else()
    message(STATUS "ALIPHYSICS sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git/)

configure_file(${PROJECT_SOURCE_DIR}/cmake/APVersion.h.in ${CMAKE_BINARY_DIR}/version/APVersion.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/version/APVersion.h DESTINATION include)