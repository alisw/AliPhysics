Macro(CreateHLTSubCMakeFiles)

  set(_dir HLT)
  
  # Check if ROOT was compile with xml and alien support
  # This is needed later to set properly the definitions during
  # compilation
  Execute_process(
    COMMAND root-config --has-xml
    OUTPUT_VARIABLE ROOT_HAS_XML
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --has-alien
    OUTPUT_VARIABLE ROOT_HAS_ALIEN
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Loop over the subdirectories, extract the package files and
  # call for each of the package files a macro which creates the
  # corresponding CMake input files in the subdirectory

  file(GLOB Package RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${_dir}/*.pkg)
  get_filename_component(__path ${_dir} NAME)
  set(EINCLUDE_SUBDIR ${__path} STEER)
  set(PACKCXXFLAGS_SUBDIR)
  set(PACKAGES)

  ParseHLTPackageFile(HLT/hlt.conf)
  set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
#  MESSAGE("${Package}")

#  set(Package HLT/libHLTrec.pkg)   
  Foreach(_pack ${Package})
    GenerateHLTInputFiles(${_pack})
    set(EINCLUDE_SUBDIR ${EINCLUDE_SUBDIR} ${EINCLUDE_PACKAGE})
    set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
  EndForeach(_pack ${Package})
  list(REMOVE_DUPLICATES EINCLUDE_SUBDIR)
  CreateHLTMainCMakeFile(${__path})

EndMacro(CreateHLTSubCMakeFiles)

Macro(GenerateHLTInputFiles Package)

  get_filename_component(_path ${Package} PATH)
  get_filename_component(_name ${Package} NAME)

  STRING(REGEX REPLACE "^(lib.*).pkg$" "\\1" _lib "${_name}")
  STRING(REGEX REPLACE "^(bin.*).pkg$" "\\1" _bin "${_name}")

  ParseHLTPackageFile(${Package})

  If(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_lib})
#    SpecialHLTSettingsMadeByHand(${_lib}) # Correct the Parser
    CreateHLTLibPackageFile(${_path}  ${_lib})
#    Message("Would create ${_lib}")
  Else(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_bin})
    CreateHLTBinPackageFile(${_path} ${_bin})
#     Message("Would create ${_bin}")
  EndIf(_name MATCHES "^lib.*$") 

EndMacro(GenerateHLTInputFiles Package)

Macro(ParseHLTPackageFile Package)

  Set(CXXSRCS_PACKAGE)
  Set(CSRCS_PACKAGE)
  Set(FSRCS_PACKAGE)
  Set(EINCLUDE_PACKAGE)
  Set(DHDR_PACKAGE)
  Set(ELIBS_PACKAGE)
  Set(HDRS_PACKAGE)
  Set(PACKCXXFLAGS_PACKAGE)
  Set(PACKFFLAGS_PACKAGE)
  Set(ADD_HEADERS)
  Set(EXCLUDE_PART FALSE)
  Set(MATCHED FALSE)
  Set(_file_glob FALSE)
  Set(_file_glob_dir)
  Set(SRCS_PACKAGE) 
  
  # Read the complete text file into variable contents

  FILE(READ "${Package}" contents)

  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  # Replace as first step ";" to "\". One "\" in a regex must be "\\\\"
  # After that replace line breaks by a semicolon, which is in cmake the
  # separator between list elements. The nice feature is that if there is a
  # follow up line this is indicated with an "\" at the end of the line
  # directly before the line break. In the generated string this two letters
  # together becomes "\;" which is not the separator between two list elements
  # but a single ";" in a liste element. With this trick one gets all
  # follow up lines in one list element which can be parsed 

  STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
  STRING(REGEX REPLACE "\n" ";" contents "${contents}")

  # Iterate over the lines of the input file
  Foreach(line ${contents})


    # Simple technique to handle if statements in the package files
    # If you hit an endif or an else at the begining of a line read
    # again lines below this point. In case of else check if it is a
    # nested "if else if else endif endif" structure. If the correct
    # case is already found don't read the following lines

    STRING(REGEX REPLACE "\t" "" line "${line}")

    If(line MATCHES "^[#].*$")
      Set(EXCLUDE_COMMENT TRUE)
#      MESSAGE("This is a comment")
    Else(line MATCHES "^[#].*$")
      Set(EXCLUDE_COMMENT FALSE)
#      MESSAGE("This is not a comment") 
    EndIf(line MATCHES "^[#].*$")

    If(line MATCHES "^endif.*$")
      Set(EXCLUDE_PART FALSE)
    EndIf(line MATCHES "^endif.*$")

    If(line MATCHES "^else.*$")
      If(NOT MATCHED)
        Set(EXCLUDE_PART FALSE)
      Else(NOT MATCHED)
        Set(EXCLUDE_PART TRUE)
      EndIf(NOT MATCHED)
    EndIf(line MATCHES "^else.*$")

    # Special handling in case ther is a if statement

    If(line MATCHES "^if.*$")

      # Special handling in case ther is a ifeq statement
 
      If(line MATCHES "^ifeq.*$")

        Set(EXCLUDE_PART TRUE)

        # Remove ifeq from text. The rest is the interesting part of
        # the statement
  
        STRING(REGEX REPLACE "ifeq" "" statement "${line}")

        # check if "ifeq" compares for ALICE_TARGET. If required target is
        # equal to ALICE_TARGET, then read the following lines until endif
        # or else is found
        # Since only one target is possible mark that a target was already
        # found. A short summary of regular expressions can be found at the
        # end of this file.

        If(line MATCHES "^.*ALICE_TARGET.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]ALICE_TARGET[)][)].*$" "\\1" _result "${line}")
          If(_result STREQUAL ${ALICE_TARGET})
            Set(MATCHED TRUE)
            Set(EXCLUDE_PART FALSE)
          Else(_result STREQUAL ${ALICE_TARGET})
            Set(EXCLUDE_PART TRUE)
          EndIf(_result STREQUAL ${ALICE_TARGET})
        EndIf(line MATCHES "^.*ALICE_TARGET.*$")


        # check if "if" statement checks for Alien or XML installation.
        # If ROOT is installed with Alien or XML switch on the compile
        # flags

        If(line MATCHES "^.*CHECK.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]CHECKXML[)][)].*$" "\\1" _result_xml "${line}") 
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]CHECKALIEN[)][)].*$" "\\1" _result_alien "${line}") 
          If(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
            Set(EXCLUDE_PART FALSE)
          Else(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
#          MESSAGE("HIER")
        EndIf(line MATCHES "^.*CHECK.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")
        
        If(line MATCHES "^.*CCM.*$")
          STRING(REGEX REPLACE "^.*[(][$][(]CCMAJORV[)],(.*)[)].*$" "\\1" _result_maj "${line}")
          STRING(REGEX REPLACE "^.*[(][$][(]CCMINORV[)],(.*)[)].*$" "\\1" _result_min "${line}")
          If(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
        EndIf(line MATCHES "^.*CCM.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

        If(line MATCHES "^.*F77.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]F77[)][)].*$" "\\1" _result_f77 "${line}")
          If(${_result_f77} STREQUAL g95)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_f77} STREQUAL g95)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_f77} STREQUAL g95)
        EndIf(line MATCHES "^.*F77.*$")


#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

        If(line MATCHES "^.*MACOSX_MINOR.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]MACOSX_MINOR[)][)].*$" "\\1" _result_mac "${line}")
          If(${_result_mac} EQUAL 5)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_mac} EQUAL 5)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_mac} EQUAL 5)
        EndIf(line MATCHES "^.*MACOSX_MINOR.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

      Else(line MATCHES "^ifeq.*$")

        If(line MATCHES "^ifdef.*$")

          Set(EXCLUDE_PART TRUE)

          # line starts with if but not ifeq        
          STRING(REGEX REPLACE "ifdef" "" statement "${line}")
 
          # Parse DATE dependent part of if statements
          If(line MATCHES "^.*DATE_ROOT.*$")
            STRING(REGEX REPLACE "^.*(DATE_ROOT).*$" "\\1" _result_date "${statement}")
            If(${_result_date} STREQUAL DATE_ROOT)
              If(DATE_FOUND)
                Set(EXCLUDE_PART FALSE)
              Else(DATE_FOUND)
                Set(EXCLUDE_PART TRUE)
              EndIf(DATE_FOUND)
            EndIf(${_result_date} STREQUAL DATE_ROOT)
          EndIf(line MATCHES "^.*DATE_ROOT.*$")

#          MESSAGE("EXCLUDE1: ${EXCLUDE_PART}")


          If(line MATCHES "^.*ALIDEBUG.*$")
            If($ENV{ALIDEBUG})
              Set(EXCLUDE_PART FALSE)
            Else($ENV{ALIDEBUG})
              Set(EXCLUDE_PART TRUE)
            EndIf($ENV{ALIDEBUG})
          EndIf(line MATCHES "^.*ALIDEBUG.*$")

#        MESSAGE("EXCLUDE2: ${EXCLUDE_PART}")
        
          If(line MATCHES "^.*ALIHLT_MLUCDIR.*$")
            If($ENV{ALIHLT_MLUCDIR})
              Set(EXCLUDE_PART FALSE)
            Else(STREQUAL $ENV{ALIHLT_MLUCDIR})
              Set(EXCLUDE_PART TRUE)
            EndIf($ENV{ALIHLT_MLUCDIR})
          EndIf(line MATCHES "^.*ALIHLT_MLUCDIR.*$")

#        MESSAGE("EXCLUDE3: ${EXCLUDE_PART}")
           

        Else(line MATCHES "^ifdef.*$")
          If(line MATCHES "^ifneq.*$")
  
            If(line MATCHES "^.*FASTJET_ROOT.*$")
              STRING(REGEX REPLACE "^.*[(][$][(]FASTJET_ROOT[)],(.*)[)].*$" "\\1" _result_fastjet "${line}")
              If(NOT ${_length_fastjet})
                Set(EXCLUDE_PART FALSE)
              Else(NOT ${_length_fastjet})
                Set(EXCLUDE_PART TRUE)
              EndIf(NOT ${_length_fastjet})

            Else(line MATCHES "^.*FASTJET_ROOT.*$")
              If(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              
                STRING(REGEX REPLACE "^.*findstring[ ](,*),[$][(]ALICE_TARGET[)].*$" "\\1" _result_macox "${line}")
                If(${_result_macox} MATCHES ".*macosx.*")
                  Set(EXCLUDE_PART FALSE)
                Else(${_result_macox} MATCHES ".*macosx.*")
                  Set(EXCLUDE_PART FALSE)
                EndIf(${_result_macox} MATCHES ".*macosx.*")
                 
              Else(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              

                MESSAGE(FATAL_ERROR "There is no rule how to handle ifneq statement in ${line}")

              EndIf(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              
            EndIf(line MATCHES "^.*FASTJET_ROOT.*$")
            
            

          Else(line MATCHES "^ifneq.*$")
            Set(EXCLUDE_PART TRUE)
            MESSAGE(FATAL_ERROR "There is no rule how to handle if statement in ${line}")
          EndIf(line MATCHES "^ifneq.*$")

        EndIf(line MATCHES "^ifdef.*$")
      EndIf(line MATCHES "^ifeq.*$")
    EndIf(line MATCHES "^if.*$")

    # If the lines are inside an if statement which is not true exclude this
    # part. In other words read only interesting part of of the file

#    MESSAGE("EXCLUDE: ${EXCLUDE_PART}, ${EXCLUDE_COMMENT}")

    if(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)
#      MESSAGE("Hier")
      If(line MATCHES "^ORGSRCS.*$") 
        STRING(REGEX REPLACE "^.*[$][(]wildcard[ ](.*)[*].*$" "\\1" orgsrc "${line}")
        Set(_file_glob_dir ${_file_glob_dir} ${orgsrc})
      EndIf(line MATCHES "^ORGSRCS.*$") 

      If(line MATCHES "^MODULE_SRCS.*$") 
        STRING(REGEX REPLACE "MODULE_SRCS" "" CXXSRCS "${line}")
         # Check if list of source files should be build from
         # list of header files. Add additional source files to
         # the list if they are defined. The parser does not 
         If("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")
            set(ADD_CXXSRCS TRUE)
         Else("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")
           if(ADD_CXXSRCS)
             set(ADD_CXXSRCS TRUE)
            else(ADD_CXXSRSC)
              set(ADD_CXXSRCS FALSE)
            endif(ADD_CXXSRCS)
         EndIf("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")


        STRING(REGEX REPLACE "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE ":=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "^;" "" CXXSRCS "${CXXSRCS}")
        SET(CXXSRCS_PACKAGE ${CXXSRCS_PACKAGE} ${CXXSRCS}) 
      EndIf(line MATCHES "^MODULE_SRCS.*$") 
  
      If(line MATCHES "^SRCS.*$")       
        If(line MATCHES patsubst)
          STRING(REGEX REPLACE "^.*[%][,](.*)[/][%][,].*$" "\\1" PACKAGE_DIR "${line}")
        Else(line MATCHES patsubst)
          STRING(REGEX REPLACE "SRCS" "" SRCS "${line}")
          STRING(REGEX REPLACE ":=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "[+]=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "[ ]+" ";" SRCS "${SRCS}")
          STRING(REGEX REPLACE "^;" "" SRCS "${SRCS}")
          SET(SRCS_PACKAGE ${SRCS_PACKAGE} ${SRCS}) 
        EndIf(line MATCHES patsubst)
      EndIf(line MATCHES "^SRCS.*$") 

      If(line MATCHES "^CSRCS.*$")       
        STRING(REGEX REPLACE "CSRCS" "" CSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "^;" "" CSRCS "${CSRCS}")
        SET(CSRCS_PACKAGE ${CSRCS_PACKAGE} ${CSRCS}) 
      EndIf(line MATCHES "^CSRCS.*$") 
 
      If(line MATCHES "^FSRCS.*$") 
        STRING(REGEX REPLACE "FSRCS" "" FSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "[+]=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "^;" "" FSRCS "${FSRCS}")
        SET(FSRCS_PACKAGE ${FSRCS_PACKAGE} ${FSRCS}) 
      EndIf(line MATCHES "^FSRCS.*$") 
  
      If(line MATCHES "^EINCLUDE.*$") 
#        MESSAGE("BLA: ${line}")
        STRING(REGEX REPLACE "EINCLUDE" "" EINCLUDE "${line}")
        STRING(REGEX REPLACE ":=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "[+]=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "[ ]+" ";" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "^;" "" EINCLUDE "${EINCLUDE}")
        SET(EINCLUDE_PACKAGE ${EINCLUDE_PACKAGE} ${EINCLUDE}) 
#        MESSAGE("EINCLUDE_PACKAGE: ${EINCLUDE_PACKAGE}")
      EndIf(line MATCHES "^EINCLUDE.*$") 
  
      If(line MATCHES "^MODULE_DHDR.*$") 
        STRING(REGEX REPLACE "MODULE_DHDR" "" DHDR "${line}")
        STRING(REGEX REPLACE "\t" "" DHDR "${DHDR}")
        STRING(STRIP ${DHDR} DHDR)
        STRING(REGEX REPLACE ":=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[+]=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[ ]+" ";" DHDR "${DHDR}")
        STRING(REGEX REPLACE "^;" "" DHDR "${DHDR}")
        SET(DHDR_PACKAGE ${DHDR_PACKAGE} ${DHDR}) 
#        MESSAGE("DHDR_PACKAGE: ${DHDR_PACKAGE}")
      EndIf(line MATCHES "^MODULE_DHDR.*$") 
  
      If(line MATCHES "^LIBHLT.*$") 
        STRING(REGEX REPLACE "^.*LIBHLT(.*)[_]VERSION.*$" "\\1" _result_library "${line}")
        STRING(REGEX REPLACE "^.*LIBHLT.*VERSION" "" LIBHLT "${line}")
        STRING(REGEX REPLACE ":=" "" LIBHLT "${LIBHLT}")
        STRING(STRIP ${LIBHLT} LIBHLT)
        set(LIBHLT_LIB_VERSION ${LIBHLT})
      EndIf(line MATCHES "^LIBHLT.*$") 

      If(line MATCHES "^PACKCXXFLAGS.*$") 
        STRING(REGEX REPLACE "PACKCXXFLAGS" "" PACKCXXFLAGS "${line}")
        STRING(REGEX REPLACE ":=" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "[+]=" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "^[ ]+[=]" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "[ ]+" ";" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "^;" ""  PACKCXXFLAGS "${PACKCXXFLAGS}")
        foreach(elem ${PACKCXXFLAGS}) 
          STRING(REGEX MATCH "^[-]D.*$" _match_result ${elem})
          if(${elem} STREQUAL "${_match_result}")
#            STRING(REGEX REPLACE "\"" "\\\\\"" PACKCXXFLAGS "${PACKCXXFLAGS}") 
            If(${elem} MATCHES LIBHLT)
#              STRING(REGEX REPLACE "[$][(].*[)]" "\\\\\"${LIBHLT_LIB_VERSION}\\\\\"" elem "${elem}")
              STRING(REGEX REPLACE "[$][(].*[)]" "${LIBHLT_LIB_VERSION}" elem "${elem}")
            EndIf(${elem} MATCHES LIBHLT)
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKCXXFLAGS}) 
  #      MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKCXXFLAGS.*$") 

     If(line MATCHES "^HLTDEFS.*$")       
        STRING(REGEX REPLACE "HLTDEFS" "" HLTDEFS "${line}")
        STRING(REGEX REPLACE ":=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "[+]=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "[ ]+" ";" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "^;" "" HLTDEFS "${HLTDEFS}")
        foreach(elem ${HLTDEFS}) 
          STRING(REGEX MATCH "^[-]D.*$" _match_result ${elem})
          If(${elem} STREQUAL "${_match_result}")
            STRING(REGEX REPLACE "VERSION" "VERSION=" elem "${elem}")
            If(${elem} MATCHES LIBHLT)
              STRING(REGEX REPLACE "[$][(].*[)]" "${LIBHLT_LIB_VERSION}" elem "${elem}")
            EndIf(${elem} MATCHES LIBHLT)
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${HLTDEFS}) 
#        MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^HLTDEFS.*$") 
   
     If(line MATCHES "^PACKFFLAGS.*$") 
        STRING(REGEX REPLACE "PACKFFLAGS" "" PACKFFLAGS "${line}")
        STRING(REGEX REPLACE ":=" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "[+]=" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "^[ ]+[=]" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "[ ]+" ";" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "^;" ""  PACKFFLAGS "${PACKFFLAGS}")
        foreach(elem ${PACKFFLAGS})
          STRING(REGEX MATCH "[-]O[%]" _match_result ${elem})
          if("${_match_result}" STREQUAL "-O%")
            set(elem "bla bla")
          endif("${_match_result}" STREQUAL "-O%")
          STRING(REGEX MATCH "^[-].*$" _match_result ${elem})
          if(${elem} STREQUAL "${_match_result}")
            SET(PACKFFLAGS_PACKAGE ${PACKFFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKFFLAGS}) 
  #      MESSAGE("FDEFINITIONS: ${PACKFFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKFFLAGS.*$") 
  
      If(line MATCHES "^ELIBS.*$") 
        If(NOT line MATCHES "^ELIBSCPP.*$")
          If(NOT line MATCHES "^ELIBSDIR.*$")
            STRING(REGEX REPLACE "ELIBS" "" ELIBS "${line}")
            STRING(REGEX REPLACE ":=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "[+]=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "[ ]+" ";" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "^;" "" ELIBS "${ELIBS}")
            SET(ELIBS_PACKAGE ${ELIBS_PACKAGE} ${ELIBS}) 
          EndIf(NOT line MATCHES "^ELIBSDIR.*$")
        EndIf(NOT line MATCHES "^ELIBSCPP.*$")
      EndIf(line MATCHES "^ELIBS.*$") 

      If(line MATCHES "^CLASS_HDRS.*$")
#          MESSAGE("HIER")
          If(NOT line MATCHES "^CLASS_HDRS_FJ.*$")
#            MESSAGE("Auch HIER")
            STRING(REGEX REPLACE "CLASS_HDRS" "" HDRS "${line}")
            STRING(REGEX REPLACE "\t" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE ":=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "[+]=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "[ ]+" ";" HDRS "${HDRS}")
            STRING(REGEX REPLACE "^;" "" HDRS "${HDRS}")
            SET(HDRS_PACKAGE ${HDRS_PACKAGE} ${HDRS}) 
          EndIf(NOT line MATCHES "^CLASS_HDRS_FJ.*$")
      EndIf(line MATCHES "^CLASS_HDRS.*$") 

#      MESSAGE("Parsed:  ${line}")
    Else(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)
#      MESSAGE("Not parsed:  ${line}")
    EndIf(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)

  EndForeach(line ${contents})
EndMacro(ParseHLTPackageFile Package)

Macro(CreateHLTMainCMakeFile path)
  set(FileName ${path}/CMakeLists.txt)
  file(WRITE ${FileName} "# -*- mode: cmake -*-\n")
  file(APPEND ${FileName} "# Create a library called \"lib<name>\" which includes the source files given in\n")
  file(APPEND ${FileName} "# the array .\n")
  file(APPEND ${FileName} "# The extension is already found.  Any number of sources could be listed here.\n\n")
  file(APPEND ${FileName} "set(INCLUDE_DIRECTORIES\n")
  foreach(file ${EINCLUDE_SUBDIR})
    If(file MATCHES "^.*[$].*$")
      STRING(REGEX REPLACE "[(]" "ENV{"  file "${file}")
      STRING(REGEX REPLACE "[)]" "}"  file "${file}")
      file(APPEND ${FileName} "${file}\n")
    Else(file MATCHES "^.*[$].*$")
      file(APPEND ${FileName} "\${CMAKE_SOURCE_DIR}/${file}\n")
    EndIf(file MATCHES "^.*[$].*$")
  endforeach(file ${EINCLUDE_SUBDIR})
  file(APPEND ${FileName} "\${ROOT_INCLUDE_DIR}\n")

  if(${_dir} STREQUAL STEER OR ${_dir} STREQUAL TPC OR ${_dir} STREQUAL ALIROOT) 
    file(APPEND ${FileName} "\${CMAKE_BINARY_DIR}/STEER\n")
  endif(${_dir} STREQUAL STEER OR ${_dir} STREQUAL TPC OR ${_dir} STREQUAL ALIROOT) 
  if(${_dir} STREQUAL RAW)
    file(APPEND ${FileName} "\${CMAKE_SOURCE_DIR}\n")
  endif(${_dir} STREQUAL RAW)


  file(APPEND ${FileName} ")\n\n")
  file(APPEND ${FileName} "include_directories(\${INCLUDE_DIRECTORIES})\n\n")
  file(APPEND ${FileName} "set(LINK_DIRECTORIES\n")
  file(APPEND ${FileName} "\${ROOT_LIBRARY_DIR}\n")
  file(APPEND ${FileName} ")\n\n")
  file(APPEND ${FileName} "link_directories(\${LINK_DIRECTORIES})\n\n")
  
  list(LENGTH PACKCXXFLAGS_SUBDIR definition_length)
  if(${definition_length} GREATER 0)
    list(REMOVE_DUPLICATES PACKCXXFLAGS_SUBDIR)
    file(APPEND  ${FileName} "Add_Definitions(\n")
    foreach(file ${PACKCXXFLAGS_SUBDIR})
#      If(${file} MATCHES LIBHLT)
#        STRING(REGEX REPLACE "[(]" "{"  file "${file}")
#        STRING(REGEX REPLACE "[)]" "}"  file "${file}")    
#        MESSAGE("BLA: ${file}")
#        file(APPEND  ${FileName} "${${file}}\n")
#
#      Else(${file} MATCHES LIBHLT)


        file(APPEND  ${FileName} "${file}\n")
#      EndIf(${file} MATCHES LIBHLT)
    endforeach(file ${CXXSRCS_SUBDIR})
    file(APPEND  ${FileName} ")\n\n")
  endif(${definition_length} GREATER 0)

  file(APPEND ${FileName} "SetModule()\n\n")
  foreach(file ${PACKAGES})
    file(APPEND ${FileName} "include (CMake_${file}.txt)\n\n")
  endforeach(file ${PACKAGES})
EndMacro(CreateHLTMainCMakeFile path)

Macro(CreateHLTLibPackageFile path lib)

  set(FileName ${path}/CMake_${lib}.txt)
  set(AddCSource FALSE)
  set(AddFortranSource FALSE)

  STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH HDRS_PACKAGE hdrs_length)
    if(${hdrs_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS\n")
      foreach(file ${HDRS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${hdrs_length} GREATER 0)


#  list(LENGTH CSRCS_PACKAGE c_length)
#  if(${c_length} GREATER 0)
#    set(AddCSource TRUE)
#    file(APPEND  ${FileName} "set(CSRCS\n")
#    foreach(file ${CSRCS_PACKAGE})
#      file(APPEND  ${FileName} "${file}\n")
#    endforeach(file ${CSRCS_PACKAGE})
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${c_length} GREATER 0)
#
#  list(LENGTH FSRCS_PACKAGE f_length)
#  if(${f_length} GREATER 0)
#    set(AddFortranSource TRUE)
#    file(APPEND  ${FileName} "set(FSRCS\n")
#    foreach(file ${FSRCS_PACKAGE})
#      file(APPEND  ${FileName} "${file}\n")
#    endforeach(file ${FSRCS_PACKAGE})
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${f_length} GREATER 0)


  if(ADD_CXXSRCS)
    file(APPEND  ${FileName} "# fill list of header files from list of source files\n")
    file(APPEND  ${FileName} "# by exchanging the file extension\n")

    file(APPEND  ${FileName} "String(REPLACE \".h\" \".cxx\" SRCS \"\${HDRS}\")\n\n")
    list(LENGTH CXXSRCS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS \${SRCS}\n")
      foreach(file ${CXXSRCS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)    
  else(ADD_CXXSRCS)
    list(LENGTH CXXSRCS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS\n")
      foreach(file ${CXXSRCS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)
  endif(ADD_CXXSRCS)

#  If(AddCSource)
#    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${CSRCS})\n\n")
#  EndIf(AddCSource)
#  If(AddFortranSource)
#    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${FSRCS})\n\n")
#  EndIf(AddFortranSource)
#
#  list(LENGTH PACKFFLAGS_PACKAGE packf_length)
#  if(${packf_length} GREATER 0)
#    file(APPEND  ${FileName} "SET_SOURCE_FILES_PROPERTIES(\n")
#    file(APPEND  ${FileName} "${FSRCS} PROPERTIES COMPILE_FLAGS\n") 
#    file(APPEND  ${FileName} "${PACKFFLAGS_PACKAGE}\n")
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${packf_length} GREATER 0)
  

#  file(APPEND  ${FileName} "AddHLTLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\" \"\${DHDR_PACKAGE}\")\n") 
  file(APPEND  ${FileName} "AddHLTLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\")\n") 

#  Message("DHDR: ${_lib}, ${DHDR_PACKAGE}")
  list(LENGTH DHDR_PACKAGE _length)
  If(${_length} EQUAL 0)
#    STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")    
    set(LinkDefFileName ${CMAKE_CURRENT_BINARY_DIR}/${path}/${_lib}LinkDef.h)
#    MESSAGE("LINKDEF: ${LinkDefFileName}")
    GenerateLinkDefFile(${LinkDefFileName} "${HDRS}")
  EndIf(${_length} EQUAL 0)

EndMacro(CreateHLTLibPackageFile path lib)

Macro(CreateHLTBinPackageFile path bin)

  set(FileName ${path}/CMake_${bin}.txt)
  STRING(REGEX REPLACE "^bin(.*)" "\\1" _bin "${bin}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH SRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(SRCS\n")
    foreach(file ${SRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${SRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  list(LENGTH ELIBS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(LIBS\n")
    foreach(file ${ELIBS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${ELIBS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  file(APPEND  ${FileName} "AddExecutable(${_bin}exe \"\${SRCS}\" \"\${LIBS}\")\n") 
EndMacro(CreateHLTBinPackageFile path bin)
 
Macro(SpecialHLTSettingsMadeByHand lib)
  If(${lib}  STREQUAL libAliengui)
    list(REMOVE_ITEM CXXSRCS_PACKAGE Aliengui/AliAnalysisGUIdummy.cxx)
  EndIf(${lib} STREQUAL libAliengui)
#  If(${lib}  STREQUAL libTPCmon)
#    list(REMOVE_ITEM CXXSRCS_PACKAGE AliTPCMonitorDateMonitor.cxx)
#    list(REMOVE_ITEM HDRS_PACKAGE AliTPCMonitorDateMonitor.h)
#  EndIf(${lib} STREQUAL libTPCmon)
  If(${lib}  STREQUAL libRAWDatabase)
    list(APPEND HDRS_PACKAGE \${ROOT_INCLUDE_DIR}/TH1F.h)
    list(APPEND H_PACKAGE \${ROOT_INCLUDE_DIR}/TH1F.h)
  EndIf(${lib} STREQUAL libRAWDatabase)
  If(${lib}  STREQUAL liblhapdf)
#    STRING(REGEX REPLACE "[=]" "\\\\=" PACKCXXFLAGS "${PACKCXXFLAGS}")
#      STRING(REGEX REPLACE "[$][(]ALICE_ROOT[)]" "\$ENV{ALICE_ROOT}" PACKCXXFLAGS "${PACKCXXFLAGS}")    
  EndIf(${lib} STREQUAL liblhapdf)
EndMacro(SpecialHLTSettingsMadeByHand lib)

Macro(GenerateLinkDefFile FileName HDRS)
  file(WRITE ${FileName} "//automatically generated ROOT DICT definition\n")
  file(APPEND ${FileName} "//!!! DO NOT EDIT THIS FILE !!!\n")
  file(APPEND ${FileName} "#ifdef __CINT__\n")
  file(APPEND ${FileName} "#pragma link off all globals;\n")
  file(APPEND ${FileName} "#pragma link off all classes;\n")
  file(APPEND ${FileName} "#pragma link off all functions;\n")
  ForEach(file ${HDRS})
    String(STRIP ${file} file)
#    MESSAGE("${file}")
    get_filename_component(_file ${file} NAME_WE)
    file(APPEND ${FileName} "#pragma link C++ class ${_file}+;\n")
  EndForEach(file ${HDRS})
  file(APPEND ${FileName} "#endif\n")
EndMacro(GenerateLinkDefFile FileName HDRS)
Macro(CreateHLTSubCMakeFiles)

  set(_dir HLT)
  
  # Check if ROOT was compile with xml and alien support
  # This is needed later to set properly the definitions during
  # compilation
  Execute_process(
    COMMAND root-config --has-xml
    OUTPUT_VARIABLE ROOT_HAS_XML
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  Execute_process(
    COMMAND root-config --has-alien
    OUTPUT_VARIABLE ROOT_HAS_ALIEN
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Loop over the subdirectories, extract the package files and
  # call for each of the package files a macro which creates the
  # corresponding CMake input files in the subdirectory

  file(GLOB Package RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${_dir}/*.pkg)
  get_filename_component(__path ${_dir} NAME)
  set(EINCLUDE_SUBDIR ${__path} STEER)
  set(PACKCXXFLAGS_SUBDIR)
  set(PACKAGES)

  ParseHLTPackageFile(HLT/hlt.conf)
  set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
#  MESSAGE("${Package}")

#  set(Package HLT/libHLTrec.pkg)   
  Foreach(_pack ${Package})
    GenerateHLTInputFiles(${_pack})
    set(EINCLUDE_SUBDIR ${EINCLUDE_SUBDIR} ${EINCLUDE_PACKAGE})
    set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
  EndForeach(_pack ${Package})
  list(REMOVE_DUPLICATES EINCLUDE_SUBDIR)
  CreateHLTMainCMakeFile(${__path})

EndMacro(CreateHLTSubCMakeFiles)

Macro(GenerateHLTInputFiles Package)

  get_filename_component(_path ${Package} PATH)
  get_filename_component(_name ${Package} NAME)

  STRING(REGEX REPLACE "^(lib.*).pkg$" "\\1" _lib "${_name}")
  STRING(REGEX REPLACE "^(bin.*).pkg$" "\\1" _bin "${_name}")

  ParseHLTPackageFile(${Package})

  If(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_lib})
#    SpecialHLTSettingsMadeByHand(${_lib}) # Correct the Parser
    CreateHLTLibPackageFile(${_path}  ${_lib})
#    Message("Would create ${_lib}")
  Else(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_bin})
    CreateHLTBinPackageFile(${_path} ${_bin})
#     Message("Would create ${_bin}")
  EndIf(_name MATCHES "^lib.*$") 

EndMacro(GenerateHLTInputFiles Package)

Macro(ParseHLTPackageFile Package)

  Set(CXXSRCS_PACKAGE)
  Set(CSRCS_PACKAGE)
  Set(FSRCS_PACKAGE)
  Set(EINCLUDE_PACKAGE)
  Set(DHDR_PACKAGE)
  Set(ELIBS_PACKAGE)
  Set(HDRS_PACKAGE)
  Set(PACKCXXFLAGS_PACKAGE)
  Set(PACKFFLAGS_PACKAGE)
  Set(ADD_HEADERS)
  Set(EXCLUDE_PART FALSE)
  Set(MATCHED FALSE)
  Set(_file_glob FALSE)
  Set(_file_glob_dir)
  Set(SRCS_PACKAGE) 
  
  # Read the complete text file into variable contents

  FILE(READ "${Package}" contents)

  # Convert file contents into a CMake list (where each element in the list
  # is one line of the file)
  # Replace as first step ";" to "\". One "\" in a regex must be "\\\\"
  # After that replace line breaks by a semicolon, which is in cmake the
  # separator between list elements. The nice feature is that if there is a
  # follow up line this is indicated with an "\" at the end of the line
  # directly before the line break. In the generated string this two letters
  # together becomes "\;" which is not the separator between two list elements
  # but a single ";" in a liste element. With this trick one gets all
  # follow up lines in one list element which can be parsed 

  STRING(REGEX REPLACE ";" "\\\\;" contents "${contents}")
  STRING(REGEX REPLACE "\n" ";" contents "${contents}")

  # Iterate over the lines of the input file
  Foreach(line ${contents})


    # Simple technique to handle if statements in the package files
    # If you hit an endif or an else at the begining of a line read
    # again lines below this point. In case of else check if it is a
    # nested "if else if else endif endif" structure. If the correct
    # case is already found don't read the following lines

    STRING(REGEX REPLACE "\t" "" line "${line}")

    If(line MATCHES "^[#].*$")
      Set(EXCLUDE_COMMENT TRUE)
#      MESSAGE("This is a comment")
    Else(line MATCHES "^[#].*$")
      Set(EXCLUDE_COMMENT FALSE)
#      MESSAGE("This is not a comment") 
    EndIf(line MATCHES "^[#].*$")

    If(line MATCHES "^endif.*$")
      Set(EXCLUDE_PART FALSE)
    EndIf(line MATCHES "^endif.*$")

    If(line MATCHES "^else.*$")
      If(NOT MATCHED)
        Set(EXCLUDE_PART FALSE)
      Else(NOT MATCHED)
        Set(EXCLUDE_PART TRUE)
      EndIf(NOT MATCHED)
    EndIf(line MATCHES "^else.*$")

    # Special handling in case ther is a if statement

    If(line MATCHES "^if.*$")

      # Special handling in case ther is a ifeq statement
 
      If(line MATCHES "^ifeq.*$")

        Set(EXCLUDE_PART TRUE)

        # Remove ifeq from text. The rest is the interesting part of
        # the statement
  
        STRING(REGEX REPLACE "ifeq" "" statement "${line}")

        # check if "ifeq" compares for ALICE_TARGET. If required target is
        # equal to ALICE_TARGET, then read the following lines until endif
        # or else is found
        # Since only one target is possible mark that a target was already
        # found. A short summary of regular expressions can be found at the
        # end of this file.

        If(line MATCHES "^.*ALICE_TARGET.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]ALICE_TARGET[)][)].*$" "\\1" _result "${line}")
          If(_result STREQUAL ${ALICE_TARGET})
            Set(MATCHED TRUE)
            Set(EXCLUDE_PART FALSE)
          Else(_result STREQUAL ${ALICE_TARGET})
            Set(EXCLUDE_PART TRUE)
          EndIf(_result STREQUAL ${ALICE_TARGET})
        EndIf(line MATCHES "^.*ALICE_TARGET.*$")


        # check if "if" statement checks for Alien or XML installation.
        # If ROOT is installed with Alien or XML switch on the compile
        # flags

        If(line MATCHES "^.*CHECK.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]CHECKXML[)][)].*$" "\\1" _result_xml "${line}") 
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]CHECKALIEN[)][)].*$" "\\1" _result_alien "${line}") 
          If(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
            Set(EXCLUDE_PART FALSE)
          Else(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_xml} STREQUAL ${ROOT_HAS_XML} OR ${_result_alien} STREQUAL ${ROOT_HAS_ALIEN})
#          MESSAGE("HIER")
        EndIf(line MATCHES "^.*CHECK.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")
        
        If(line MATCHES "^.*CCM.*$")
          STRING(REGEX REPLACE "^.*[(][$][(]CCMAJORV[)],(.*)[)].*$" "\\1" _result_maj "${line}")
          STRING(REGEX REPLACE "^.*[(][$][(]CCMINORV[)],(.*)[)].*$" "\\1" _result_min "${line}")
          If(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_maj} EQUAL 4 OR ${_result_min} EQUAL 3)
        EndIf(line MATCHES "^.*CCM.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

        If(line MATCHES "^.*F77.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]F77[)][)].*$" "\\1" _result_f77 "${line}")
          If(${_result_f77} STREQUAL g95)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_f77} STREQUAL g95)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_f77} STREQUAL g95)
        EndIf(line MATCHES "^.*F77.*$")


#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

        If(line MATCHES "^.*MACOSX_MINOR.*$")
          STRING(REGEX REPLACE "^.*[(](.*),[$][(]MACOSX_MINOR[)][)].*$" "\\1" _result_mac "${line}")
          If(${_result_mac} EQUAL 5)
            Set(EXCLUDE_PART FALSE)
          Else(${_result_mac} EQUAL 5)
            Set(EXCLUDE_PART TRUE)
          EndIf(${_result_mac} EQUAL 5)
        EndIf(line MATCHES "^.*MACOSX_MINOR.*$")

#        MESSAGE("EXCLUDE: ${EXCLUDE_PART}")

      Else(line MATCHES "^ifeq.*$")

        If(line MATCHES "^ifdef.*$")

          Set(EXCLUDE_PART TRUE)

          # line starts with if but not ifeq        
          STRING(REGEX REPLACE "ifdef" "" statement "${line}")
 
          # Parse DATE dependent part of if statements
          If(line MATCHES "^.*DATE_ROOT.*$")
            STRING(REGEX REPLACE "^.*(DATE_ROOT).*$" "\\1" _result_date "${statement}")
            If(${_result_date} STREQUAL DATE_ROOT)
              If(DATE_FOUND)
                Set(EXCLUDE_PART FALSE)
              Else(DATE_FOUND)
                Set(EXCLUDE_PART TRUE)
              EndIf(DATE_FOUND)
            EndIf(${_result_date} STREQUAL DATE_ROOT)
          EndIf(line MATCHES "^.*DATE_ROOT.*$")

#          MESSAGE("EXCLUDE1: ${EXCLUDE_PART}")


          If(line MATCHES "^.*ALIDEBUG.*$")
            If($ENV{ALIDEBUG})
              Set(EXCLUDE_PART FALSE)
            Else($ENV{ALIDEBUG})
              Set(EXCLUDE_PART TRUE)
            EndIf($ENV{ALIDEBUG})
          EndIf(line MATCHES "^.*ALIDEBUG.*$")

#        MESSAGE("EXCLUDE2: ${EXCLUDE_PART}")
        
          If(line MATCHES "^.*ALIHLT_MLUCDIR.*$")
            If($ENV{ALIHLT_MLUCDIR})
              Set(EXCLUDE_PART FALSE)
            Else(STREQUAL $ENV{ALIHLT_MLUCDIR})
              Set(EXCLUDE_PART TRUE)
            EndIf($ENV{ALIHLT_MLUCDIR})
          EndIf(line MATCHES "^.*ALIHLT_MLUCDIR.*$")

#        MESSAGE("EXCLUDE3: ${EXCLUDE_PART}")
           

        Else(line MATCHES "^ifdef.*$")
          If(line MATCHES "^ifneq.*$")
  
            If(line MATCHES "^.*FASTJET_ROOT.*$")
              STRING(REGEX REPLACE "^.*[(][$][(]FASTJET_ROOT[)],(.*)[)].*$" "\\1" _result_fastjet "${line}")
              If(NOT ${_length_fastjet})
                Set(EXCLUDE_PART FALSE)
              Else(NOT ${_length_fastjet})
                Set(EXCLUDE_PART TRUE)
              EndIf(NOT ${_length_fastjet})

            Else(line MATCHES "^.*FASTJET_ROOT.*$")
              If(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              
                STRING(REGEX REPLACE "^.*findstring[ ](,*),[$][(]ALICE_TARGET[)].*$" "\\1" _result_macox "${line}")
                If(${_result_macox} MATCHES ".*macosx.*")
                  Set(EXCLUDE_PART FALSE)
                Else(${_result_macox} MATCHES ".*macosx.*")
                  Set(EXCLUDE_PART FALSE)
                EndIf(${_result_macox} MATCHES ".*macosx.*")
                 
              Else(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              

                MESSAGE(FATAL_ERROR "There is no rule how to handle ifneq statement in ${line}")

              EndIf(line MATCHES "^.*findstring macosx,[$][(]ALICE_TARGET[)].*$")              
            EndIf(line MATCHES "^.*FASTJET_ROOT.*$")
            
            

          Else(line MATCHES "^ifneq.*$")
            Set(EXCLUDE_PART TRUE)
            MESSAGE(FATAL_ERROR "There is no rule how to handle if statement in ${line}")
          EndIf(line MATCHES "^ifneq.*$")

        EndIf(line MATCHES "^ifdef.*$")
      EndIf(line MATCHES "^ifeq.*$")
    EndIf(line MATCHES "^if.*$")

    # If the lines are inside an if statement which is not true exclude this
    # part. In other words read only interesting part of of the file

#    MESSAGE("EXCLUDE: ${EXCLUDE_PART}, ${EXCLUDE_COMMENT}")

    if(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)
#      MESSAGE("Hier")
      If(line MATCHES "^ORGSRCS.*$") 
        STRING(REGEX REPLACE "^.*[$][(]wildcard[ ](.*)[*].*$" "\\1" orgsrc "${line}")
        Set(_file_glob_dir ${_file_glob_dir} ${orgsrc})
      EndIf(line MATCHES "^ORGSRCS.*$") 

      If(line MATCHES "^MODULE_SRCS.*$") 
        STRING(REGEX REPLACE "MODULE_SRCS" "" CXXSRCS "${line}")
         # Check if list of source files should be build from
         # list of header files. Add additional source files to
         # the list if they are defined. The parser does not 
         If("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")
            set(ADD_CXXSRCS TRUE)
         Else("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")
           if(ADD_CXXSRCS)
             set(ADD_CXXSRCS TRUE)
            else(ADD_CXXSRSC)
              set(ADD_CXXSRCS FALSE)
            endif(ADD_CXXSRCS)
         EndIf("${CXXSRCS}" MATCHES "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]")


        STRING(REGEX REPLACE "[$][(]CLASS_HDRS:[.]h=[.]cxx[)]" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE ":=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "^;" "" CXXSRCS "${CXXSRCS}")
        SET(CXXSRCS_PACKAGE ${CXXSRCS_PACKAGE} ${CXXSRCS}) 
      EndIf(line MATCHES "^MODULE_SRCS.*$") 
  
      If(line MATCHES "^SRCS.*$")       
        If(line MATCHES patsubst)
          STRING(REGEX REPLACE "^.*[%][,](.*)[/][%][,].*$" "\\1" PACKAGE_DIR "${line}")
        Else(line MATCHES patsubst)
          STRING(REGEX REPLACE "SRCS" "" SRCS "${line}")
          STRING(REGEX REPLACE ":=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "[+]=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "=" "" SRCS "${SRCS}")
          STRING(REGEX REPLACE "[ ]+" ";" SRCS "${SRCS}")
          STRING(REGEX REPLACE "^;" "" SRCS "${SRCS}")
          SET(SRCS_PACKAGE ${SRCS_PACKAGE} ${SRCS}) 
        EndIf(line MATCHES patsubst)
      EndIf(line MATCHES "^SRCS.*$") 

      If(line MATCHES "^CSRCS.*$")       
        STRING(REGEX REPLACE "CSRCS" "" CSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "=" "" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CSRCS "${CSRCS}")
        STRING(REGEX REPLACE "^;" "" CSRCS "${CSRCS}")
        SET(CSRCS_PACKAGE ${CSRCS_PACKAGE} ${CSRCS}) 
      EndIf(line MATCHES "^CSRCS.*$") 
 
      If(line MATCHES "^FSRCS.*$") 
        STRING(REGEX REPLACE "FSRCS" "" FSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "[+]=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "=" "" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" FSRCS "${FSRCS}")
        STRING(REGEX REPLACE "^;" "" FSRCS "${FSRCS}")
        SET(FSRCS_PACKAGE ${FSRCS_PACKAGE} ${FSRCS}) 
      EndIf(line MATCHES "^FSRCS.*$") 
  
      If(line MATCHES "^EINCLUDE.*$") 
#        MESSAGE("BLA: ${line}")
        STRING(REGEX REPLACE "EINCLUDE" "" EINCLUDE "${line}")
        STRING(REGEX REPLACE ":=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "[+]=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "=" "" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "[ ]+" ";" EINCLUDE "${EINCLUDE}")
        STRING(REGEX REPLACE "^;" "" EINCLUDE "${EINCLUDE}")
        SET(EINCLUDE_PACKAGE ${EINCLUDE_PACKAGE} ${EINCLUDE}) 
#        MESSAGE("EINCLUDE_PACKAGE: ${EINCLUDE_PACKAGE}")
      EndIf(line MATCHES "^EINCLUDE.*$") 
  
      If(line MATCHES "^MODULE_DHDR.*$") 
        STRING(REGEX REPLACE "MODULE_DHDR" "" DHDR "${line}")
        STRING(REGEX REPLACE "\t" "" DHDR "${DHDR}")
        STRING(STRIP ${DHDR} DHDR)
        STRING(REGEX REPLACE ":=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[+]=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[ ]+" ";" DHDR "${DHDR}")
        STRING(REGEX REPLACE "^;" "" DHDR "${DHDR}")
        SET(DHDR_PACKAGE ${DHDR_PACKAGE} ${DHDR}) 
#        MESSAGE("DHDR_PACKAGE: ${DHDR_PACKAGE}")
      EndIf(line MATCHES "^MODULE_DHDR.*$") 
  
      If(line MATCHES "^LIBHLT.*$") 
        STRING(REGEX REPLACE "^.*LIBHLT(.*)[_]VERSION.*$" "\\1" _result_library "${line}")
        STRING(REGEX REPLACE "^.*LIBHLT.*VERSION" "" LIBHLT "${line}")
        STRING(REGEX REPLACE ":=" "" LIBHLT "${LIBHLT}")
        STRING(STRIP ${LIBHLT} LIBHLT)
        set(LIBHLT_LIB_VERSION ${LIBHLT})
      EndIf(line MATCHES "^LIBHLT.*$") 

      If(line MATCHES "^PACKCXXFLAGS.*$") 
        STRING(REGEX REPLACE "PACKCXXFLAGS" "" PACKCXXFLAGS "${line}")
        STRING(REGEX REPLACE ":=" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "[+]=" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "^[ ]+[=]" "" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "[ ]+" ";" PACKCXXFLAGS "${PACKCXXFLAGS}")
        STRING(REGEX REPLACE "^;" ""  PACKCXXFLAGS "${PACKCXXFLAGS}")
        foreach(elem ${PACKCXXFLAGS}) 
          STRING(REGEX MATCH "^[-]D.*$" _match_result ${elem})
          if(${elem} STREQUAL "${_match_result}")
#            STRING(REGEX REPLACE "\"" "\\\\\"" PACKCXXFLAGS "${PACKCXXFLAGS}") 
            If(${elem} MATCHES LIBHLT)
#              STRING(REGEX REPLACE "[$][(].*[)]" "\\\\\"${LIBHLT_LIB_VERSION}\\\\\"" elem "${elem}")
              STRING(REGEX REPLACE "[$][(].*[)]" "${LIBHLT_LIB_VERSION}" elem "${elem}")
            EndIf(${elem} MATCHES LIBHLT)
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKCXXFLAGS}) 
  #      MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKCXXFLAGS.*$") 

     If(line MATCHES "^HLTDEFS.*$")       
        STRING(REGEX REPLACE "HLTDEFS" "" HLTDEFS "${line}")
        STRING(REGEX REPLACE ":=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "[+]=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "=" "" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "[ ]+" ";" HLTDEFS "${HLTDEFS}")
        STRING(REGEX REPLACE "^;" "" HLTDEFS "${HLTDEFS}")
        foreach(elem ${HLTDEFS}) 
          STRING(REGEX MATCH "^[-]D.*$" _match_result ${elem})
          If(${elem} STREQUAL "${_match_result}")
            STRING(REGEX REPLACE "VERSION" "VERSION=" elem "${elem}")
            If(${elem} MATCHES LIBHLT)
              STRING(REGEX REPLACE "[$][(].*[)]" "${LIBHLT_LIB_VERSION}" elem "${elem}")
            EndIf(${elem} MATCHES LIBHLT)
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${HLTDEFS}) 
#        MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^HLTDEFS.*$") 
   
     If(line MATCHES "^PACKFFLAGS.*$") 
        STRING(REGEX REPLACE "PACKFFLAGS" "" PACKFFLAGS "${line}")
        STRING(REGEX REPLACE ":=" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "[+]=" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "^[ ]+[=]" "" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "[ ]+" ";" PACKFFLAGS "${PACKFFLAGS}")
        STRING(REGEX REPLACE "^;" ""  PACKFFLAGS "${PACKFFLAGS}")
        foreach(elem ${PACKFFLAGS})
          STRING(REGEX MATCH "[-]O[%]" _match_result ${elem})
          if("${_match_result}" STREQUAL "-O%")
            set(elem "bla bla")
          endif("${_match_result}" STREQUAL "-O%")
          STRING(REGEX MATCH "^[-].*$" _match_result ${elem})
          if(${elem} STREQUAL "${_match_result}")
            SET(PACKFFLAGS_PACKAGE ${PACKFFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKFFLAGS}) 
  #      MESSAGE("FDEFINITIONS: ${PACKFFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKFFLAGS.*$") 
  
      If(line MATCHES "^ELIBS.*$") 
        If(NOT line MATCHES "^ELIBSCPP.*$")
          If(NOT line MATCHES "^ELIBSDIR.*$")
            STRING(REGEX REPLACE "ELIBS" "" ELIBS "${line}")
            STRING(REGEX REPLACE ":=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "[+]=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "=" "" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "[ ]+" ";" ELIBS "${ELIBS}")
            STRING(REGEX REPLACE "^;" "" ELIBS "${ELIBS}")
            SET(ELIBS_PACKAGE ${ELIBS_PACKAGE} ${ELIBS}) 
          EndIf(NOT line MATCHES "^ELIBSDIR.*$")
        EndIf(NOT line MATCHES "^ELIBSCPP.*$")
      EndIf(line MATCHES "^ELIBS.*$") 

      If(line MATCHES "^CLASS_HDRS.*$")
#          MESSAGE("HIER")
          If(NOT line MATCHES "^CLASS_HDRS_FJ.*$")
#            MESSAGE("Auch HIER")
            STRING(REGEX REPLACE "CLASS_HDRS" "" HDRS "${line}")
            STRING(REGEX REPLACE "\t" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE ":=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "[+]=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "=" "" HDRS "${HDRS}")
            STRING(REGEX REPLACE "[ ]+" ";" HDRS "${HDRS}")
            STRING(REGEX REPLACE "^;" "" HDRS "${HDRS}")
            SET(HDRS_PACKAGE ${HDRS_PACKAGE} ${HDRS}) 
          EndIf(NOT line MATCHES "^CLASS_HDRS_FJ.*$")
      EndIf(line MATCHES "^CLASS_HDRS.*$") 

#      MESSAGE("Parsed:  ${line}")
    Else(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)
#      MESSAGE("Not parsed:  ${line}")
    EndIf(NOT EXCLUDE_PART AND NOT EXCLUDE_COMMENT)

  EndForeach(line ${contents})
EndMacro(ParseHLTPackageFile Package)

Macro(CreateHLTMainCMakeFile path)
  set(FileName ${path}/CMakeLists.txt)
  file(WRITE ${FileName} "# -*- mode: cmake -*-\n")
  file(APPEND ${FileName} "# Create a library called \"lib<name>\" which includes the source files given in\n")
  file(APPEND ${FileName} "# the array .\n")
  file(APPEND ${FileName} "# The extension is already found.  Any number of sources could be listed here.\n\n")
  file(APPEND ${FileName} "set(INCLUDE_DIRECTORIES\n")
  foreach(file ${EINCLUDE_SUBDIR})
    If(file MATCHES "^.*[$].*$")
      STRING(REGEX REPLACE "[(]" "ENV{"  file "${file}")
      STRING(REGEX REPLACE "[)]" "}"  file "${file}")
      file(APPEND ${FileName} "${file}\n")
    Else(file MATCHES "^.*[$].*$")
      file(APPEND ${FileName} "\${CMAKE_SOURCE_DIR}/${file}\n")
    EndIf(file MATCHES "^.*[$].*$")
  endforeach(file ${EINCLUDE_SUBDIR})
  file(APPEND ${FileName} "\${ROOT_INCLUDE_DIR}\n")

  if(${_dir} STREQUAL STEER OR ${_dir} STREQUAL TPC OR ${_dir} STREQUAL ALIROOT) 
    file(APPEND ${FileName} "\${CMAKE_BINARY_DIR}/STEER\n")
  endif(${_dir} STREQUAL STEER OR ${_dir} STREQUAL TPC OR ${_dir} STREQUAL ALIROOT) 
  if(${_dir} STREQUAL RAW)
    file(APPEND ${FileName} "\${CMAKE_SOURCE_DIR}\n")
  endif(${_dir} STREQUAL RAW)


  file(APPEND ${FileName} ")\n\n")
  file(APPEND ${FileName} "include_directories(\${INCLUDE_DIRECTORIES})\n\n")
  file(APPEND ${FileName} "set(LINK_DIRECTORIES\n")
  file(APPEND ${FileName} "\${ROOT_LIBRARY_DIR}\n")
  file(APPEND ${FileName} ")\n\n")
  file(APPEND ${FileName} "link_directories(\${LINK_DIRECTORIES})\n\n")
  
  list(LENGTH PACKCXXFLAGS_SUBDIR definition_length)
  if(${definition_length} GREATER 0)
    list(REMOVE_DUPLICATES PACKCXXFLAGS_SUBDIR)
    file(APPEND  ${FileName} "Add_Definitions(\n")
    foreach(file ${PACKCXXFLAGS_SUBDIR})
#      If(${file} MATCHES LIBHLT)
#        STRING(REGEX REPLACE "[(]" "{"  file "${file}")
#        STRING(REGEX REPLACE "[)]" "}"  file "${file}")    
#        MESSAGE("BLA: ${file}")
#        file(APPEND  ${FileName} "${${file}}\n")
#
#      Else(${file} MATCHES LIBHLT)


        file(APPEND  ${FileName} "${file}\n")
#      EndIf(${file} MATCHES LIBHLT)
    endforeach(file ${CXXSRCS_SUBDIR})
    file(APPEND  ${FileName} ")\n\n")
  endif(${definition_length} GREATER 0)

  file(APPEND ${FileName} "SetModule()\n\n")
  foreach(file ${PACKAGES})
    file(APPEND ${FileName} "include (CMake_${file}.txt)\n\n")
  endforeach(file ${PACKAGES})
EndMacro(CreateHLTMainCMakeFile path)

Macro(CreateHLTLibPackageFile path lib)

  set(FileName ${path}/CMake_${lib}.txt)
  set(AddCSource FALSE)
  set(AddFortranSource FALSE)

  STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH HDRS_PACKAGE hdrs_length)
    if(${hdrs_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS\n")
      foreach(file ${HDRS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${hdrs_length} GREATER 0)


#  list(LENGTH CSRCS_PACKAGE c_length)
#  if(${c_length} GREATER 0)
#    set(AddCSource TRUE)
#    file(APPEND  ${FileName} "set(CSRCS\n")
#    foreach(file ${CSRCS_PACKAGE})
#      file(APPEND  ${FileName} "${file}\n")
#    endforeach(file ${CSRCS_PACKAGE})
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${c_length} GREATER 0)
#
#  list(LENGTH FSRCS_PACKAGE f_length)
#  if(${f_length} GREATER 0)
#    set(AddFortranSource TRUE)
#    file(APPEND  ${FileName} "set(FSRCS\n")
#    foreach(file ${FSRCS_PACKAGE})
#      file(APPEND  ${FileName} "${file}\n")
#    endforeach(file ${FSRCS_PACKAGE})
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${f_length} GREATER 0)


  if(ADD_CXXSRCS)
    file(APPEND  ${FileName} "# fill list of header files from list of source files\n")
    file(APPEND  ${FileName} "# by exchanging the file extension\n")

    file(APPEND  ${FileName} "String(REPLACE \".h\" \".cxx\" SRCS \"\${HDRS}\")\n\n")
    list(LENGTH CXXSRCS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS \${SRCS}\n")
      foreach(file ${CXXSRCS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)    
  else(ADD_CXXSRCS)
    list(LENGTH CXXSRCS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS\n")
      foreach(file ${CXXSRCS_PACKAGE})
        String(STRIP ${file} file)
        file(APPEND  ${FileName} "${PACKAGE_DIR}/${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)
  endif(ADD_CXXSRCS)

#  If(AddCSource)
#    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${CSRCS})\n\n")
#  EndIf(AddCSource)
#  If(AddFortranSource)
#    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${FSRCS})\n\n")
#  EndIf(AddFortranSource)
#
#  list(LENGTH PACKFFLAGS_PACKAGE packf_length)
#  if(${packf_length} GREATER 0)
#    file(APPEND  ${FileName} "SET_SOURCE_FILES_PROPERTIES(\n")
#    file(APPEND  ${FileName} "${FSRCS} PROPERTIES COMPILE_FLAGS\n") 
#    file(APPEND  ${FileName} "${PACKFFLAGS_PACKAGE}\n")
#    file(APPEND  ${FileName} ")\n\n")
#  endif(${packf_length} GREATER 0)
  

#  file(APPEND  ${FileName} "AddHLTLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\" \"\${DHDR_PACKAGE}\")\n") 
  file(APPEND  ${FileName} "AddHLTLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\")\n") 

#  Message("DHDR: ${_lib}, ${DHDR_PACKAGE}")
  list(LENGTH DHDR_PACKAGE _length)
  If(${_length} EQUAL 0)
#    STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")    
    set(LinkDefFileName ${CMAKE_CURRENT_BINARY_DIR}/${path}/${_lib}LinkDef.h)
#    MESSAGE("LINKDEF: ${LinkDefFileName}")
    GenerateLinkDefFile(${LinkDefFileName} "${HDRS}")
  EndIf(${_length} EQUAL 0)

EndMacro(CreateHLTLibPackageFile path lib)

Macro(CreateHLTBinPackageFile path bin)

  set(FileName ${path}/CMake_${bin}.txt)
  STRING(REGEX REPLACE "^bin(.*)" "\\1" _bin "${bin}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH SRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(SRCS\n")
    foreach(file ${SRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${SRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  list(LENGTH ELIBS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(LIBS\n")
    foreach(file ${ELIBS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${ELIBS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  file(APPEND  ${FileName} "AddExecutable(${_bin}exe \"\${SRCS}\" \"\${LIBS}\")\n") 
EndMacro(CreateHLTBinPackageFile path bin)
 
Macro(SpecialHLTSettingsMadeByHand lib)
  If(${lib}  STREQUAL libAliengui)
    list(REMOVE_ITEM CXXSRCS_PACKAGE Aliengui/AliAnalysisGUIdummy.cxx)
  EndIf(${lib} STREQUAL libAliengui)
#  If(${lib}  STREQUAL libTPCmon)
#    list(REMOVE_ITEM CXXSRCS_PACKAGE AliTPCMonitorDateMonitor.cxx)
#    list(REMOVE_ITEM HDRS_PACKAGE AliTPCMonitorDateMonitor.h)
#  EndIf(${lib} STREQUAL libTPCmon)
  If(${lib}  STREQUAL libRAWDatabase)
    list(APPEND HDRS_PACKAGE \${ROOT_INCLUDE_DIR}/TH1F.h)
    list(APPEND H_PACKAGE \${ROOT_INCLUDE_DIR}/TH1F.h)
  EndIf(${lib} STREQUAL libRAWDatabase)
  If(${lib}  STREQUAL liblhapdf)
#    STRING(REGEX REPLACE "[=]" "\\\\=" PACKCXXFLAGS "${PACKCXXFLAGS}")
#      STRING(REGEX REPLACE "[$][(]ALICE_ROOT[)]" "\$ENV{ALICE_ROOT}" PACKCXXFLAGS "${PACKCXXFLAGS}")    
  EndIf(${lib} STREQUAL liblhapdf)
EndMacro(SpecialHLTSettingsMadeByHand lib)

Macro(GenerateLinkDefFile FileName HDRS)
  file(WRITE ${FileName} "//automatically generated ROOT DICT definition\n")
  file(APPEND ${FileName} "//!!! DO NOT EDIT THIS FILE !!!\n")
  file(APPEND ${FileName} "#ifdef __CINT__\n")
  file(APPEND ${FileName} "#pragma link off all globals;\n")
  file(APPEND ${FileName} "#pragma link off all classes;\n")
  file(APPEND ${FileName} "#pragma link off all functions;\n")
  ForEach(file ${HDRS})
    String(STRIP ${file} file)
#    MESSAGE("${file}")
    get_filename_component(_file ${file} NAME_WE)
    file(APPEND ${FileName} "#pragma link C++ class ${_file}+;\n")
  EndForEach(file ${HDRS})
  file(APPEND ${FileName} "#endif\n")
EndMacro(GenerateLinkDefFile FileName HDRS)
