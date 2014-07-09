# AliRoot Build System CMakeLists for EVE
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

if(ZEROMQ_FOUND)
set(PACKAGES "${PACKAGES}"
    CMakebinalifakedisplay.pkg
    CMakebinalistorageadmin.pkg
    CMakebinalistorage.pkg
    CMakelibStorage.pkg
)
ENDIF(ZEROMQ_FOUND)

ALICE_BuildModule()
