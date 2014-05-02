# AliRoot Build System CMakeLists for MONITOR
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This will only build MONITOR only if ZeroMQ is found
set (PACKAGES
  CMakebinalitestproc.pkg
  CMakebinderoot.pkg
  CMakebinmonitorCheck.pkg
  CMakebinmonitorGDC.pkg
  CMakelibMONITOR.pkg
)

if(ZEROMQ_FOUND)
set(PACKAGES "${PACKAGES}"
  CMakelibMONITORzmq.pkg
  CMakebinalieventserver.pkg
)
ENDIF(ZEROMQ_FOUND)

ALICE_BuildModule()
