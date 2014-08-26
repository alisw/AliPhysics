# AliRoot Build System CMakeLists for MONITOR
#
# Author: Johny Jose m(johny.jose@cern.ch)
#         Port of previous Makefile build to cmake

cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This will only build MONITOR only if ZeroMQ is found

set (PACKAGES
  CMakelibMONITOR.pkg
  alitestproc/CMakebinalitestproc.pkg
  deroot/CMakebinderoot.pkg
  monitorCheck/CMakebinmonitorCheck.pkg
  monitorGDC/CMakebinmonitorGDC.pkg 
)

if(ZEROMQ_FOUND)
set(PACKAGES "${PACKAGES}"
  CMakelibMONITORzmq.pkg
  alistoragemanager/CMakebinalistorage.pkg
  alistoragemanager/CMakebinalistorageadmin.pkg
  alistoragemanager/CMakebinalifakedisplay.pkg
  alistoragemanager/CMakebinalistorageminimalexample.pkg
  alistoragemanager/CMakebinalistoragegetpoints.pkg
  alieventserver/CMakebinalieventserver.pkg
)
ENDIF(ZEROMQ_FOUND)

ALICE_BuildModule()
