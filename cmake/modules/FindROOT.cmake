# -*- mode: cmake -*-
# - Finds ROOT instalation
# This module sets up ROOT information 
# We suppose root-config to be in the PATH. Otherwise we stop.

Find_program(ROOT_CONFIG root-config)

If (${ROOT_CONFIG} MATCHES "ROOT_CONFIG-NOTFOUND")
  Set(ROOT_FOUND FALSE)
  Message(STATUS "Install Root and make sure it is in the PATH")

Else (${ROOT_CONFIG} MATCHES "ROOT_CONFIG-NOTFOUND")  
  
  Set(ROOT_FOUND TRUE)

  Execute_process(
    COMMAND root-config --prefix 
    OUTPUT_VARIABLE ROOTSYS 
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --arch
    OUTPUT_VARIABLE ALICE_TARGET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --f77 
    OUTPUT_VARIABLE _f77 
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  Set(ENV{F77} ${_f77})

  Execute_process(
    COMMAND root-config --cc
    OUTPUT_VARIABLE _cc 
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  Set(ENV{CC} ${_cc})

  Execute_process(
    COMMAND root-config --cxx
    OUTPUT_VARIABLE _cxx
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  Set(ENV{CXX} ${_cxx})

  Execute_process(
    COMMAND root-config --version 
    OUTPUT_VARIABLE ROOT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --glibs
    OUTPUT_VARIABLE ROOT_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Find_program(ROOTCINT rootcint)
  If(NOT ROOTCINT)
    Message(FATAL_ERROR "Found ROOT but not rootcint, your ROOT installation is corrupted")
  EndIf(NOT ROOTCINT)

  Set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lThread -lMinuit -lHtml -lVMC -lEG -lGeom -lTreePlayer -lXMLIO -lProof)
  Set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lProofPlayer -lMLP -lSpectrum -lEve -lRGL -lGed -lXMLParser -lPhysics)
  Set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)

  # Make variables changeble to the advanced user
  Mark_as_advanced(ROOT_LIBRARY_DIR ROOT_INCLUDE_DIR ROOT_DEFINITIONS)

  Set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${ROOT_LIBRARY_DIR})

  Message(STATUS "Found Root ${ROOT_VERSION} in ${ROOTSYS}/bin/root")   

# we need at least version 5.00/00
  If (NOT ROOT_MIN_VERSION)
    Set(ROOT_MIN_VERSION "5.00/00")
  Endif (NOT ROOT_MIN_VERSION)
   
  # now parse the parts of the user given version string into variables
  String(REGEX REPLACE "^([0-9]+)\\.[0-9][0-9]+\\/[0-9][0-9]+" "\\1" req_root_major_vers "${ROOT_MIN_VERSION}")
  String(REGEX REPLACE "^[0-9]+\\.([0-9][0-9])+\\/[0-9][0-9]+.*" "\\1" req_root_minor_vers "${ROOT_MIN_VERSION}")
  String(REGEX REPLACE "^[0-9]+\\.[0-9][0-9]+\\/([0-9][0-9]+)" "\\1" req_root_patch_vers "${ROOT_MIN_VERSION}")
   
  # and now the version string given by qmake
  String(REGEX REPLACE "^([0-9]+)\\.[0-9][0-9]+\\/[0-9][0-9]+.*" "\\1" found_root_major_vers "${ROOT_VERSION}")
  String(REGEX REPLACE "^[0-9]+\\.([0-9][0-9])+\\/[0-9][0-9]+.*" "\\1" found_root_minor_vers "${ROOT_VERSION}")
  String(REGEX REPLACE "^[0-9]+\\.[0-9][0-9]+\\/([0-9][0-9]+).*" "\\1" found_root_patch_vers "${ROOT_VERSION}")

  If (found_root_major_vers LESS 5)
    Message(FATAL_ERROR "Invalid ROOT version \"${ROOT_VERSION}\", at least major version 4 is required, e.g. \"5.00/00\"")
  Endif(found_root_major_vers LESS 5)

  # compute an overall version number which can be compared at once
  Math(EXPR req_vers "${req_root_major_vers}*10000 + ${req_root_minor_vers}*100 + ${req_root_patch_vers}")
  Math(EXPR found_vers "${found_root_major_vers}*10000 + ${found_root_minor_vers}*100 + ${found_root_patch_vers}")
   
  If (found_vers LESS req_vers)
    Set(ROOT_FOUND FALSE)
    Set(ROOT_INSTALLED_VERSION_TOO_OLD TRUE)
  Else (found_vers LESS req_vers)
    Set(ROOT_FOUND TRUE)
  Endif (found_vers LESS req_vers)

Endif (${ROOT_CONFIG} MATCHES "ROOT_CONFIG-NOTFOUND")  


#####################################################################################


Macro(ROOT_GENERATE_DICTIONARY INFILES LINKDEF_FILE OUTFILE INCLUDE_DIRS_IN)
 
  Set(_special_settings "${ARGV4}")
  Set(INCLUDE_DIRS)
  Set(infiles_nopath)

  Foreach (_current_FILE ${INCLUDE_DIRS_IN})
    Set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})   
  Endforeach (_current_FILE ${INCLUDE_DIRS_IN})
 
  String(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" bla "${OUTFILE}")
  Set(OUTFILES ${OUTFILE} ${bla})

  Foreach (_current_FILE ${INFILES})
    Get_filename_component(name_wo_path ${_current_FILE} NAME)
    Set(infiles_nopath ${infiles_nopath} ${name_wo_path})   
  Endforeach (_current_FILE ${INFILES})

  Get_property(_defs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY COMPILE_DEFINITIONS)
  Set(_ddefs)
  Foreach (_def ${_defs})
    Set(_ddefs "${_ddefs} -D${_def}")
  Endforeach (_def ${_defs})
  Separate_arguments(_ddefs)

  Add_custom_command(OUTPUT ${OUTFILES}
     COMMAND DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:${ROOT_LIBRARY_DIR} ${ROOTCINT}
     ARGS -f ${OUTFILE} -c -DHAVE_CONFIG_H ${_ddefs} ${_special_settings} ${INCLUDE_DIRS} ${infiles_nopath} ${LINKDEF_FILE} 
     DEPENDS ${INFILES} ${LINKDEF_FILE})

Endmacro(ROOT_GENERATE_DICTIONARY)
