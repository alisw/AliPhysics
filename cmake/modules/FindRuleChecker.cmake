# -*- mode: cmake -*-

# - Find IRST Code Analysis Tool (RuleChecker)
# Finds if the RuleChecker is installed and sets the following variables:
#
# RULE_CHECKER_FOUND = Boolean defining if Rule Checker is available
# RULE_CHECKER_SRCML = PATH TO src2srcml executable
# RULE_CHECKER_JAVA  = PATH TO java executable 
# RULE_CHECKER_JAR   = PATH TO NewRuleChecker.jar 
# RULE_CHECKER_RULES = Rules (defaults to ALICE rules!)

# We suppose that perl and java are present and working
Find_Program(RULE_CHECKER_SRCML src2srcml $ENV{PATH})
Find_Program(RULE_CHECKER_JAVA java $ENV{PATH})
#Message(STATUS "${RULE_CHECKER_SRCML} ${RULE_CHECKER_JAVA}")
Set(RULE_CHECKER_FOUND NO)
If(RULE_CHECKER_SRCML AND RULE_CHECKER_JAVA)
  Find_File(RULE_CHECKER_JAR NAMES NewRuleChecker.jar PATHS $ENV{ALICE}/local/ALICENewRuleChecker/NewRuleChecker)
  Find_File(RULE_CHECKER_RULES AliceCodingConventions.xml $ENV{ALICE}/local/ALICENewRuleChecker/NewRuleChecker/config)
  #Message(STATUS "${RULE_CHECKER_JAR} ${RULE_CHECKER_RULES}")
  If(RULE_CHECKER_JAR AND RULE_CHECKER_RULES)
    Set(RULE_CHECKER_FOUND YES)
    Message(STATUS "RuleChecker installed on this system")
  EndIf(RULE_CHECKER_JAR AND RULE_CHECKER_RULES)
EndIf(RULE_CHECKER_SRCML AND RULE_CHECKER_JAVA)
If(NOT RULE_CHECKER_FOUND)
  Message(STATUS "RuleChecker not installed on this system")
EndIf(NOT RULE_CHECKER_FOUND)

#______________________________________________________________________________________________________________
Function(CheckViols LIB SRCS)

  If(RULE_CHECKER_FOUND)
    Set(_FactExt $ENV{ALICE}/local/ALICENewRuleChecker/FactExtractor/FactExtractor.jar)
    Set(_RuleCheck ${RULE_CHECKER_JAR})

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
      String(REGEX REPLACE "([^;]*).cxx" "${_checkDir}/\\1.viol"    _violFile "${_checkFile}")
      String(REGEX REPLACE "([^;]*).cxx" "\\1.h" _checkHead "${_checkFile}")
      Add_Custom_Command(
	OUTPUT ${_violFile}
	COMMAND src2srcml 
        ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${_checkFile} ${_srcxml}
	COMMAND src2srcml 
        ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${_checkHead} ${_hdrxml}
	COMMAND java -jar
	ARGS ${_FactExt} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND java -Xmx500m -jar ${RULE_CHECKER_JAR} 
        ARGS ${_srcxml} ${_hdrxml} ${CMAKE_CURRENT_SOURCE_DIR}/factFile.xml ${RULE_CHECKER_RULES} > ${_violFile}
	DEPENDS ${_checkFile}
	WORKING_DIRECTORY ${_checkDir})
      Set(VIOLS ${VIOLS} ${_violFile})
    EndForeach(_checkFile ${CXXSRCS})
    
    Add_Custom_Target(check-${LIB} DEPENDS ${VIOLS})
    Add_Dependencies(check-all check-${LIB})

  EndIf(RULE_CHECKER_FOUND)
    
EndFunction (CheckViols)

