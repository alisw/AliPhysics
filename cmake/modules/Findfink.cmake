# -*- mode: cmake -*-
# - Finds fink instalation

If(APPLE)
  Find_program(FINK fink)

  If (${FINK} MATCHES "ROOT_CONFIG-NOTFOUND")
    Set(FINK_FOUND FALSE)
    Message(SEND_ERROR "Could not find fink, will continue but probably fail...")
    Set(FINK_ROOT "/usr/local")
  Else (${FINK} MATCHES "ROOT_CONFIG-NOTFOUND")
    String(REPLACE "/bin/fink" "" FINK_ROOT ${FINK})
    Message(STATUS "Found fink in ${FINK_ROOT}")
  EndIf (${FINK} MATCHES "ROOT_CONFIG-NOTFOUND")

Else(APPLE)
  Message(FATAL "fink only needed on Mac, something is very wrong...")
Endif(APPLE)

