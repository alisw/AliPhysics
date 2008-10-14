# -*- mode: cmake -*-

# - Find IRST Code Analysis Tool (RuleChecker)
# Finds if the RuleChecker is installed and sets the following variables:
#
# RULE_CHECKER_FOUND = Boolean defining if Rule Checker is available
# RULE_CHECKER_PATH  = Path to the java class containing the ALICE Rule Checker
# RULE_CHECKER       = Rules (defaults to ALICE rules!)
# RULE_CHECKER_PATCH = Perl script to be executed to patch preprocessor files
#

# We suppose that perl and java are present and working

Find_File(RULE_CHECKER_PATCH patch4alice.prl
  /usr/local/IRST/patch ${ALICE}/local/IRST/patch /opt/IRST/patch)

If(RULE_CHECKER_PATCH)
  Set(RULE_CHECKER_FOUND YES)
  Set(RULE_CHECKER rules.ALICE.ALICERuleChecker)
  String(REPLACE "/patch/patch4alice.prl" "" 
    RULE_CHECKER_PATH "${RULE_CHECKER_PATCH}")
  
Else(RULE_CHECKER_PATCH)
  Set(RULE_CHECKER_FOUND NO)
  Message(STATUS "RuleChecker not installed on this system")
  
EndIf(RULE_CHECKER_PATCH)


#_______________________________________________________________________________
Function(CheckViols LIB SRCS)

  If(RULE_CHECKER_FOUND)
    String(REGEX MATCHALL "[^ ]*.cxx" CXXSRCS "${SRCS}")
    
    Set(_checkDir "${CMAKE_CURRENT_BINARY_DIR}/check")
    
    If(NOT EXISTS ${_checkDir})
      File(MAKE_DIRECTORY ${_checkDir})
    EndIf(NOT EXISTS ${_checkDir})
    
    Set(_inc_dirs)
    Foreach(_dir ${INCLUDE_DIRECTORIES})
      Set(_inc_dirs ${_inc_dirs} -I${_dir})
    EndForeach(_dir ${INCLUDE_DIRECTORIES})
    
    Set(VIOLS)
    Foreach(_checkFile ${CXXSRCS})
      Get_Filename_Component(_checkFilePath ${_checkFile} PATH)
      If(NOT EXISTS ${_checkDir}/${_checkFilePath})
	File(MAKE_DIRECTORY ${_checkDir}/${_checkFilePath})
      EndIf(NOT EXISTS ${_checkDir}/${_checkFilePath})
      String(REGEX REPLACE "([^;]*).cxx" "${_checkDir}/\\1.viol" _violFile "${_checkFile}")
      String(REGEX REPLACE "([^;]*).cxx" "\\1.i" _tempFile_in ${_checkFile})
      Add_Custom_Command(
	OUTPUT ${_violFile}
	COMMAND ${CMAKE_CXX_COMPILER} 
        ARGS -E ${CMAKE_CXX_FLAGS} ${_inc_dirs} -I. ${CMAKE_CURRENT_SOURCE_DIR}/${_checkFile} > ${_tempFile_in}
	COMMAND ${RULE_CHECKER_PATCH} ARGS ${_tempFile_in}
	COMMAND CLASSPATH=${RULE_CHECKER_PATH} java ${RULE_CHECKER} 
        ARGS ${_tempFile_in} ${CMAKE_CURRENT_SOURCE_DIR} > ${_violFile}
	DEPENDS ${_checkFile}
	WORKING_DIRECTORY ${_checkDir})
      Set(VIOLS ${VIOLS} ${_violFile})
    EndForeach(_checkFile ${CXXSRCS})
    
    Add_Custom_Target(check-${LIB} DEPENDS ${VIOLS})
    Add_Dependencies(check-all check-${LIB})

  EndIf(RULE_CHECKER_FOUND)
    
EndFunction (CheckViols)

#_______________________________________________________________________________
Function(CheckViolsNew LIB SRCS)

  If(RULE_CHECKER_FOUND)
    Set(_FactExt $ENV{ALICE}/local/ALICENewRuleChecker/FactExtractor/FactExtractor.jar)
    Set(_RuleCheck $ENV{ALICE}/local/ALICENewRuleChecker/NewRuleChecker/NewRuleChecker.jar)

    String(REGEX MATCHALL "[^ ]*.cxx" CXXSRCS "${SRCS}")
    
    Set(_checkDir "${CMAKE_CURRENT_BINARY_DIR}/check_new")
    
    If(NOT EXISTS ${_checkDir})
      File(MAKE_DIRECTORY ${_checkDir})
    EndIf(NOT EXISTS ${_checkDir})
    
    Set(_inc_dirs)
    Foreach(_dir ${INCLUDE_DIRECTORIES})
      Set(_inc_dirs ${_inc_dirs} -I${_dir})
    EndForeach(_dir ${INCLUDE_DIRECTORIES})
    
    Set(VIOLS)
    Foreach(_checkFile ${CXXSRCS})
      Get_Filename_Component(_checkFilePath ${_checkFile} PATH)
      If(NOT EXISTS ${_checkDir}/${_checkFilePath})
	File(MAKE_DIRECTORY ${_checkDir}/${_checkFilePath})
      EndIf(NOT EXISTS ${_checkDir}/${_checkFilePath})
      String(REGEX REPLACE "([^;]*).cxx" "${_checkDir}/\\1.cxx.xml" _srcxml "${_checkFile}")
      String(REGEX REPLACE "([^;]*).cxx" "${_checkDir}/\\1.h.xml"   _hdrxml "${_checkFile}")
      String(REGEX REPLACE "([^;]*).cxx" "${_checkDir}/\\1.h.viol"   _violFile "${_checkFile}")
      String(REGEX REPLACE "([^;]*).cxx" "\\1.h" _checkHead "${_checkFile}")
      Add_Custom_Command(
	OUTPUT ${_violFile}
	COMMAND src2srcml 
        ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${_checkFile} ${_srcxml}
	COMMAND src2srcml 
        ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${_checkHead} ${_hdrxml}
	COMMAND java
	ARGS -jar ${_FactExt} 
	COMMAND ${RULE_CHECKER_PATCH} ARGS ${_tempFile_in}
	COMMAND CLASSPATH=${RULE_CHECKER_PATH} java ${RULE_CHECKER} 
        ARGS ${_tempFile_in} ${CMAKE_CURRENT_SOURCE_DIR} > ${_violFile}
	DEPENDS ${_checkFile}
	WORKING_DIRECTORY ${_checkDir})
      Set(VIOLS ${VIOLS} ${_violFile})
    EndForeach(_checkFile ${CXXSRCS})
    
    Add_Custom_Target(check-${LIB} DEPENDS ${VIOLS})
    Add_Dependencies(check-all check-${LIB})

  EndIf(RULE_CHECKER_FOUND)
    
EndFunction (CheckViolsNew)

