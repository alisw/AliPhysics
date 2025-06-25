####################################################################################
# Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 #
# All rights reserved.                                                             #
#                                                                                  #
# Redistribution and use in source and binary forms, with or without               #
# modification, are permitted provided that the following conditions are met:      #
#     * Redistributions of source code must retain the above copyright             #
#       notice, this list of conditions and the following disclaimer.              #
#     * Redistributions in binary form must reproduce the above copyright          #
#       notice, this list of conditions and the following disclaimer in the        #
#       documentation and/or other materials provided with the distribution.       #
#     * Neither the name of the <organization> nor the                             #
#       names of its contributors may be used to endorse or promote products       #
#       derived from this software without specific prior written permission.      #
#                                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           #
# DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      #
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       #
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     #
####################################################################################
# Find module for KFParticle 
# Copied and adapted from FindRooUnfold.cmake by Dario Berzano <dario.berzano@cern.ch>
#
# From your `CMakeLists.txt` use it as simply as:
#
#   find_package(KFParticle)
#
# It expects `-DKFPARTICLE=<prefix>` to be passed as argument, and it falls back
# to the environment variable `KFPARTICLE_ROOT`.
#
# A target is exported if KFParticle is found. The target takes care of defining
# all the library paths, compile flags, and include paths. Usage is very simple:
#
#   if(KFParticle_FOUND)
#     target_link_library(MyAnalysisClass MyDep1 MyDep2 KFParticle::KFParticle)
#   endif()
#
# No additional command is required (no `include_directories()`, etc.) and no
# variable other than `KFParticle_FOUND` is exported.

include(FindPackageHandleStandardArgs)

if(NOT DEFINED KFPARTICLE)
  set(KFPARTICLE "$ENV{KFPARTICLE_ROOT}")
endif()

find_library(KFPARTICLE_LIBRARY KFParticle
             PATHS "${KFPARTICLE}/lib" NO_DEFAULT_PATH)
find_path(KFPARTICLE_INCLUDE_DIR KFParticle.h
          PATH_SUFFIXES "include" "include/KFParticle"
          HINTS "${KFPARTICLE}" NO_DEFAULT_PATH)

find_package_handle_standard_args(KFParticle DEFAULT_MSG
                                  KFPARTICLE_LIBRARY KFPARTICLE_INCLUDE_DIR)

# Set KFParticle::KFParticle target
add_library(KFParticle::KFParticle SHARED IMPORTED)
set_target_properties(KFParticle::KFParticle PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${KFPARTICLE_INCLUDE_DIR}"
  IMPORTED_LOCATION             "${KFPARTICLE_LIBRARY}"
  INTERFACE_COMPILE_DEFINITIONS "WITH_KFPARTICLE")

# Unset RooUnfold variables
unset(KFPARTICLE_LIBPATH)
unset(KFPARTICLE_INCLUDE_DIR)