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

# Configure ARVerion.h using Git informatiion
# Sets 4 git variables
#  - GIT_REFSPEC - complete name of the current reference
#  - ALIROOT_BRANCH - name of the branch or tag extracted from the current reference
#  - GIT_SHA1 - current hash in the long format
#  - GIT_SHORT_SHA1 - current hash in the short format
#
#  - ALIROOT_VERSION - name of the branch/tag
#  - ALIROOT_VERSION_RPM - name of the branch/tag in rpm format, - replaced with .
#  - ALIROOT_REVISION - short sha1
if(EXISTS ${PROJECT_SOURCE_DIR}/.git/)
    include(GetGitRevisionDescription)
    
    find_package(Git)
    
    if(GIT_FOUND)
        message(STATUS "Git version = ${GIT_VERSION_STRING}")
        
        get_git_head_revision(GIT_REFSPEC GIT_SHA1)

        # generate the short version of the revision hash
        execute_process(COMMAND git rev-parse --short ${GIT_SHA1}
                          WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
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
                            WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            RESULT_VARIABLE revcount
                            OUTPUT_VARIABLE ALIROOT_SERIAL_ORIGINAL)
        else()
            # generate the short version of the revision hash using -wc -l
            execute_process(COMMAND git rev-list ${GIT_SHA1}
                            COMMAND wc -l
                            WORKING_DIRECTORY ${AliRoot_SOURCE_DIR}
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            RESULT_VARIABLE revcount
                            OUTPUT_VARIABLE ALIROOT_SERIAL_ORIGINAL)
        endif()

        # GIT_REFSPEC is empty for detached mode = tags in detached mode or checkout to specific hash

        # returns the closest reference to the current hash
        # name of the current tag or heads/branch in the case of branches
        # Older git version of Git report only the name of the tag
        # Newer version report tags/vAN-20141215
        # Just in case we replace tags/ with nothing
        git_describe(ALIROOT_GIT_TAG "--all" "--abbrev=0")
        
        if(ALIROOT_GIT_TAG)
            string(STRIP ${ALIROOT_GIT_TAG} ALIROOT_GIT_TAG)
            string(REPLACE "tags/" ""  ALIROOT_GIT_TAG ${ALIROOT_GIT_TAG})
        endif(ALIROOT_GIT_TAG)
        
        # using the closest tag for branches
        git_describe(ALIROOT_CLOSEST_GIT_TAG "--abbrev=0")

        if(ALIROOT_GIT_TAG)
            string(STRIP ${ALIROOT_GIT_TAG} ALIROOT_GIT_TAG)
            string(REPLACE "tags/" ""  ALIROOT_GIT_TAG ${ALIROOT_GIT_TAG})
        endif(ALIROOT_GIT_TAG)
        
        STRING(REGEX REPLACE "^(.+/)(.+)/(.*)$" "\\2" BRANCH_TYPE "${GIT_REFSPEC}" )
        
        # the revision is not set in the case of a branch, it means we are doing development
        # and the revision will trigger a reconfiguration
        if(BRANCH_TYPE STREQUAL "heads")
            set(ALIROOT_REVISION "ThisIsaBranchNoRevisionProvided")
            set(ALIROOT_SERIAL 0)
            set(ALIROOT_GIT_TAG ${ALIROOT_CLOSEST_GIT_TAG})
            STRING(REGEX REPLACE "^(.+/)(.+/)(.*)$" "\\3" SHORT_BRANCH "${GIT_REFSPEC}" )
            message(STATUS "This is a working branch, ARVersion will not contain the revision and the serial number")
        else()
            set(SHORT_BRANCH ${ALIROOT_GIT_TAG})
            set(ALIROOT_REVISION ${GIT_SHORT_SHA1})
            set(ALIROOT_SERIAL ${ALIROOT_SERIAL_ORIGINAL})
        endif()

        set(ALIROOT_BRANCH ${SHORT_BRANCH})
        set(ALIROOT_VERSION ${SHORT_BRANCH})
        
        # extract major minor and patch from AliRoot tag
        string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+" "\\1" ALIROOT_VERSION_MAJOR "${ALIROOT_GIT_TAG}")
        string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+" "\\1" ALIROOT_VERSION_MINOR "${ALIROOT_GIT_TAG}")
        string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)" "\\1" ALIROOT_VERSION_PATCH "${ALIROOT_GIT_TAG}")
        message(STATUS "Found ALIROOT version ${ALIROOT_VERSION_MAJOR}.${ALIROOT_VERSION_MINOR}.${ALIROOT_VERSION_PATCH}")
        
        # Replace - with . for rpm creation
        string(REPLACE "-" "." ALIROOT_VERSION_RPM ${ALIROOT_VERSION})

        message(STATUS "Aliroot branch/tag: \"${ALIROOT_VERSION}\" - Revision:  \"${GIT_SHORT_SHA1}\" - Serial: \"${ALIROOT_SERIAL_ORIGINAL}\"")

    else()
        message(STATUS "Git not installed. I can't tell you which revision you are using!")
    endif(GIT_FOUND)
else()
    message("AliRoot sources not downloaded from a Version Control System. I can't tell which revision you are using!")
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git/)

configure_file(${PROJECT_SOURCE_DIR}/cmake/ARVersion.h.tmp ${CMAKE_BINARY_DIR}/version/ARVersion.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/version/ARVersion.h DESTINATION include)