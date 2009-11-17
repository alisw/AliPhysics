Macro(CreateSubCMakeFiles)

  # Extract list of subdirectories for which the script has to
  # create CMakeLists.txt out of the information of the main
  # CMakeLists.txt

  MESSAGE(STATUS "Generating CMakeLists.txt in subdirectories from the package files.")

  file(STRINGS ${CMAKE_SOURCE_DIR}/CMakeLists.txt bla 
              REGEX "^Add_subdirectory(.*)"
      )

  Foreach(_dir ${bla})
    STRING(REGEX REPLACE "^Add_subdirectory\((.*)\)$" "\\1" _OutDir "${_dir}")
    STRING(STRIP ${_OutDir} _OutDir1)
    string(LENGTH ${_OutDir1} _Length)
    math(EXPR _Length1 ${_Length}-2)
    string(SUBSTRING ${_OutDir1} 1 ${_Length1} _OutDir)
    List(APPEND Directories ${_OutDir})
  EndForeach(_dir ${bla})

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

  Foreach(_dir ${Directories})
    if(NOT ${_dir} STREQUAL HLT) 
      file(GLOB Package RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${_dir}/*.pkg)
      get_filename_component(__path ${_dir} NAME)
      set(EINCLUDE_SUBDIR ${__path} STEER)
      set(PACKCXXFLAGS_SUBDIR)
      set(PACKAGES)
      Foreach(_pack ${Package})
        GenerateInputFiles(${_pack})
        set(EINCLUDE_SUBDIR ${EINCLUDE_SUBDIR} ${EINCLUDE_PACKAGE})
        set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
      EndForeach(_pack ${Package})
      list(REMOVE_DUPLICATES EINCLUDE_SUBDIR)
      CreateMainCMakeFile(${__path})
    else(NOT ${_dir} STREQUAL HLT) 
#      MESSAGE(STATUS "Don't generate files for HLT do to complete different syntax of package files")
       Include(CreateHLTSubCMakeFiles)
       CreateHLTSubCMakeFiles()
    endif(NOT ${_dir} STREQUAL HLT) 
  EndForeach(_dir ${Directories})


EndMacro(CreateSubCMakeFiles)

Macro(GenerateInputFiles Package)

  get_filename_component(_path ${Package} PATH)
  get_filename_component(_name ${Package} NAME)

  STRING(REGEX REPLACE "^(lib.*).pkg$" "\\1" _lib "${_name}")
  STRING(REGEX REPLACE "^(bin.*).pkg$" "\\1" _bin "${_name}")

  ParsePackageFile(${Package})

  If(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_lib})
    SpecialSettingsMadeByHand(${_lib}) # Correct the Parser
    CreateLibPackageFile(${_path}  ${_lib})
  Else(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_bin})
    CreateBinPackageFile(${_path} ${_bin})
  EndIf(_name MATCHES "^lib.*$") 

EndMacro(GenerateInputFiles Package)

Macro(ParsePackageFile Package)

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
  Set(DIMDIR $ENV{DIMDIR})
  
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

 #   If(line MATCHES "^[#].*$")
 #     Set(EXCLUDE_PART TRUE)
 #   EndIf(line MATCHES "^[#].*$")

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
          #MESSAGE("LINE: ${line}")

          # line starts with if but not ifeq        

          # Parse DATE dependent part of if statements
          STRING(REGEX REPLACE "^.*(DATE_ROOT).*$" "\\1" _result_date "${line}")
          If(${_result_date} STREQUAL DATE_ROOT)
            If(DATE_FOUND)
              Set(EXCLUDE_PART FALSE)
            Else(DATE_FOUND)
              Set(EXCLUDE_PART TRUE)
            EndIf(DATE_FOUND)
          EndIf(${_result_date} STREQUAL DATE_ROOT)

          STRING(REGEX REPLACE "^.*(DIMDIR).*$" "\\1" _result_date "${line}")
          If(${_result_date} STREQUAL DIMDIR)
            If(DIMDIR)
              Set(EXCLUDE_PART FALSE)
            Else(DIMDIR)
              Set(EXCLUDE_PART TRUE)
            EndIf(DIMDIR)
          EndIf(${_result_date} STREQUAL DIMDIR)

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

#    MESSAGE("EXCLUDE*: ${EXCLUDE_PART}")

    if(NOT EXCLUDE_PART)

      If(line MATCHES "^ORGSRCS.*$") 
        STRING(REGEX REPLACE "^.*[$][(]wildcard[ ](.*)[*].*$" "\\1" orgsrc "${line}")
        Set(_file_glob_dir ${_file_glob_dir} ${orgsrc})
      EndIf(line MATCHES "^ORGSRCS.*$") 

      If(line MATCHES "^SRCS.*$") 
        STRING(REGEX REPLACE "SRCS" "" CXXSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "^;" "" CXXSRCS "${CXXSRCS}")
        If(CXXSRCS MATCHES "^.*patsubst.*$") 
          set(_file_glob TRUE)
         Else(CXXSRCS MATCHES "^.*patsubst.*$") 
          SET(CXXSRCS_PACKAGE ${CXXSRCS_PACKAGE} ${CXXSRCS}) 
        EndIf(CXXSRCS MATCHES "^.*patsubst.*$") 
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
  
      If(line MATCHES "^DHDR.*$") 
        STRING(REGEX REPLACE "DHDR" "" DHDR "${line}")
        STRING(REGEX REPLACE ":=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[+]=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[ ]+" ";" DHDR "${DHDR}")
        STRING(REGEX REPLACE "^;" "" DHDR "${DHDR}")
        SET(DHDR_PACKAGE ${DHDR_PACKAGE} ${DHDR}) 
      EndIf(line MATCHES "^DHDR.*$") 
  
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
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKCXXFLAGS}) 
  #      MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKCXXFLAGS.*$") 
  
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

      If(line MATCHES "^HDRS.*$")
          STRING(REGEX REPLACE "HDRS" "" HDRS "${line}")
    
          # Check if list of header files should be build from
          # list of source files. Add additional header files to
          # the list if they are defined. The parser does not 
          IF("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
            set(ADD_HEADERS TRUE)
          ELSE("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
             if(ADD_HEADERS)
               set(ADD_HEADERS TRUE)
             else(ADD_HEADERS)
               set(ADD_HEADERS FALSE)
             endif(ADD_HEADERS)
          ENDIF("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
          STRING(REGEX REPLACE "[$][(]SRCS:[.]cxx=[.]h[)]" "" HDRS "${HDRS}")

#          STRING(REGEX REPLACE "[$][(]filter-out(.*),.*[)]" "//1" _exclude_h "${HDRS}")
#          STRING(LENGTH ${_exclude_h} _exclude_h_length)
#          If(${_exclude_h_length} GREATER 0)
#            String(STRIP ${_exclude_h} _exclude_h) 
#            list(REMOVE_ITEM HDRS ${_exclude_h})
#          EndIf(${_exclude_h_length} GREATER 0)
 
          STRING(REGEX REPLACE "[$][(]filter-out.*[)]" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE ":=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "[+]=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "[ ]+" ";" HDRS "${HDRS}")
          STRING(REGEX REPLACE "^;" "" HDRS "${HDRS}")
          SET(HDRS_PACKAGE ${HDRS_PACKAGE} ${HDRS}) 
      EndIf(line MATCHES "^HDRS.*$") 

#      MESSAGE("Parsed:  ${line}")
    Else(NOT EXCLUDE_PART)
#      MESSAGE("Not parsed:  ${line}")
    EndIf(NOT EXCLUDE_PART)

  EndForeach(line ${contents})
EndMacro(ParsePackageFile Package)

Macro(CreateMainCMakeFile path)
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
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_SUBDIR})
    file(APPEND  ${FileName} ")\n\n")
  endif(${definition_length} GREATER 0)

  file(APPEND ${FileName} "SetModule()\n\n")
  foreach(file ${PACKAGES})
    file(APPEND ${FileName} "include (CMake_${file}.txt)\n\n")
  endforeach(file ${PACKAGES})
EndMacro(CreateMainCMakeFile path)

Macro(CreateLibPackageFile path lib)

  set(FileName ${path}/CMake_${lib}.txt)
  set(AddCSource FALSE)
  set(AddFortranSource FALSE)

  STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  If(_file_glob)
    set(_counter 0)
    ForEach(_file ${_file_glob_dir})
      file(APPEND  ${FileName} "File(GLOB SRCS${_counter} RELATIVE \${CMAKE_CURRENT_SOURCE_DIR} \${CMAKE_SOURCE_DIR}/${_file}*.cxx)\n\n")
      Math(EXPR _counter ${_counter}+1)
    EndForEach(_file ${_file_glob_dir})
    Math(EXPR _counter ${_counter}-1)
    file(APPEND  ${FileName} "set(SRCS)\n")
    ForEach(_loop RANGE ${_counter})
      file(APPEND  ${FileName} "set(SRCS \${SRCS} \${SRCS${_loop}})\n")
    EndForEach(_loop RANGE ${_counter})
  Else(_file_glob)
    list(LENGTH CXXSRCS_PACKAGE cxx_length)
    if(${cxx_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS\n")
      foreach(file ${CXXSRCS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${cxx_length} GREATER 0)
  EndIf(_file_glob)

  list(LENGTH CSRCS_PACKAGE c_length)
  if(${c_length} GREATER 0)
    set(AddCSource TRUE)
    file(APPEND  ${FileName} "set(CSRCS\n")
    foreach(file ${CSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${c_length} GREATER 0)

  list(LENGTH FSRCS_PACKAGE f_length)
  if(${f_length} GREATER 0)
    set(AddFortranSource TRUE)
    file(APPEND  ${FileName} "set(FSRCS\n")
    foreach(file ${FSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${FSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${f_length} GREATER 0)


  if(ADD_HEADERS)
    file(APPEND  ${FileName} "# fill list of header files from list of source files\n")
    file(APPEND  ${FileName} "# by exchanging the file extension\n")

    file(APPEND  ${FileName} "String(REPLACE \".cxx\" \".h\" HDRS \"\${SRCS}\")\n\n")
    list(LENGTH HDRS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS \${HDRS}\n")
      foreach(file ${HDRS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)    
  else(ADD_HEADERS)
    list(LENGTH HDRS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS\n")
      foreach(file ${HDRS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)
  endif(ADD_HEADERS)

  If(AddCSource)
    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${CSRCS})\n\n")
  EndIf(AddCSource)
  If(AddFortranSource)
    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${FSRCS})\n\n")
  EndIf(AddFortranSource)

  list(LENGTH PACKFFLAGS_PACKAGE packf_length)
  if(${packf_length} GREATER 0)
    file(APPEND  ${FileName} "SET_SOURCE_FILES_PROPERTIES(\n")
    file(APPEND  ${FileName} "${FSRCS} PROPERTIES COMPILE_FLAGS\n") 
    file(APPEND  ${FileName} "${PACKFFLAGS_PACKAGE}\n")
    file(APPEND  ${FileName} ")\n\n")
  endif(${packf_length} GREATER 0)
  

  file(APPEND  ${FileName} "AddLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\")\n") 

  
  # If package file is empty remove the CMake input file for that packge and
  # remove the package from the list.
#  If( ${cxx_length} EQUAL 0 AND ${c_length} EQUAL 0 AND ${f_length} EQUAL 0)
#    file(REMOVE ${FileName})
#    list(REMOVE_ITEM PACKAGES lib${_lib})
#  EndIf( ${cxx_length} EQUAL 0 AND ${c_length} EQUAL 0 AND ${f_length} EQUAL 0)

EndMacro(CreateLibPackageFile path lib)

Macro(CreateBinPackageFile path bin)

  set(FileName ${path}/CMake_${bin}.txt)
  STRING(REGEX REPLACE "^bin(.*)" "\\1" _bin "${bin}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH CXXSRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(SRCS\n")
    foreach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  list(LENGTH CSRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(CSRCS\n")
    foreach(file ${CSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_PACKAGE})
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
EndMacro(CreateBinPackageFile path bin)
 
Macro(SpecialSettingsMadeByHand lib)
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
EndMacro(SpecialSettingsMadeByHand lib)
Macro(CreateSubCMakeFiles)

  # Extract list of subdirectories for which the script has to
  # create CMakeLists.txt out of the information of the main
  # CMakeLists.txt

  MESSAGE(STATUS "Generating CMakeLists.txt in subdirectories from the package files.")

  file(STRINGS ${CMAKE_SOURCE_DIR}/CMakeLists.txt bla 
              REGEX "^Add_subdirectory(.*)"
      )

  Foreach(_dir ${bla})
    STRING(REGEX REPLACE "^Add_subdirectory\((.*)\)$" "\\1" _OutDir "${_dir}")
    STRING(STRIP ${_OutDir} _OutDir1)
    string(LENGTH ${_OutDir1} _Length)
    math(EXPR _Length1 ${_Length}-2)
    string(SUBSTRING ${_OutDir1} 1 ${_Length1} _OutDir)
    List(APPEND Directories ${_OutDir})
  EndForeach(_dir ${bla})

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

  Foreach(_dir ${Directories})
    if(NOT ${_dir} STREQUAL HLT) 
      file(GLOB Package RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${_dir}/*.pkg)
      get_filename_component(__path ${_dir} NAME)
      set(EINCLUDE_SUBDIR ${__path} STEER)
      set(PACKCXXFLAGS_SUBDIR)
      set(PACKAGES)
      Foreach(_pack ${Package})
        GenerateInputFiles(${_pack})
        set(EINCLUDE_SUBDIR ${EINCLUDE_SUBDIR} ${EINCLUDE_PACKAGE})
        set(PACKCXXFLAGS_SUBDIR ${PACKCXXFLAGS_SUBDIR} ${PACKCXXFLAGS_PACKAGE})
      EndForeach(_pack ${Package})
      list(REMOVE_DUPLICATES EINCLUDE_SUBDIR)
      CreateMainCMakeFile(${__path})
    else(NOT ${_dir} STREQUAL HLT) 
#      MESSAGE(STATUS "Don't generate files for HLT do to complete different syntax of package files")
       Include(CreateHLTSubCMakeFiles)
       CreateHLTSubCMakeFiles()
    endif(NOT ${_dir} STREQUAL HLT) 
  EndForeach(_dir ${Directories})


EndMacro(CreateSubCMakeFiles)

Macro(GenerateInputFiles Package)

  get_filename_component(_path ${Package} PATH)
  get_filename_component(_name ${Package} NAME)

  STRING(REGEX REPLACE "^(lib.*).pkg$" "\\1" _lib "${_name}")
  STRING(REGEX REPLACE "^(bin.*).pkg$" "\\1" _bin "${_name}")

  ParsePackageFile(${Package})

  If(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_lib})
    SpecialSettingsMadeByHand(${_lib}) # Correct the Parser
    CreateLibPackageFile(${_path}  ${_lib})
  Else(_name MATCHES "^lib.*$") 
    Set(PACKAGES ${PACKAGES} ${_bin})
    CreateBinPackageFile(${_path} ${_bin})
  EndIf(_name MATCHES "^lib.*$") 

EndMacro(GenerateInputFiles Package)

Macro(ParsePackageFile Package)

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
  Set(DIMDIR $ENV{DIMDIR})
  
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

 #   If(line MATCHES "^[#].*$")
 #     Set(EXCLUDE_PART TRUE)
 #   EndIf(line MATCHES "^[#].*$")

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
          #MESSAGE("LINE: ${line}")

          # line starts with if but not ifeq        

          # Parse DATE dependent part of if statements
          STRING(REGEX REPLACE "^.*(DATE_ROOT).*$" "\\1" _result_date "${line}")
          If(${_result_date} STREQUAL DATE_ROOT)
            If(DATE_FOUND)
              Set(EXCLUDE_PART FALSE)
            Else(DATE_FOUND)
              Set(EXCLUDE_PART TRUE)
            EndIf(DATE_FOUND)
          EndIf(${_result_date} STREQUAL DATE_ROOT)

          STRING(REGEX REPLACE "^.*(DIMDIR).*$" "\\1" _result_date "${line}")
          If(${_result_date} STREQUAL DIMDIR)
            If(DIMDIR)
              Set(EXCLUDE_PART FALSE)
            Else(DIMDIR)
              Set(EXCLUDE_PART TRUE)
            EndIf(DIMDIR)
          EndIf(${_result_date} STREQUAL DIMDIR)

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

#    MESSAGE("EXCLUDE*: ${EXCLUDE_PART}")

    if(NOT EXCLUDE_PART)

      If(line MATCHES "^ORGSRCS.*$") 
        STRING(REGEX REPLACE "^.*[$][(]wildcard[ ](.*)[*].*$" "\\1" orgsrc "${line}")
        Set(_file_glob_dir ${_file_glob_dir} ${orgsrc})
      EndIf(line MATCHES "^ORGSRCS.*$") 

      If(line MATCHES "^SRCS.*$") 
        STRING(REGEX REPLACE "SRCS" "" CXXSRCS "${line}")
        STRING(REGEX REPLACE ":=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[+]=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "=" "" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "[ ]+" ";" CXXSRCS "${CXXSRCS}")
        STRING(REGEX REPLACE "^;" "" CXXSRCS "${CXXSRCS}")
        If(CXXSRCS MATCHES "^.*patsubst.*$") 
          set(_file_glob TRUE)
         Else(CXXSRCS MATCHES "^.*patsubst.*$") 
          SET(CXXSRCS_PACKAGE ${CXXSRCS_PACKAGE} ${CXXSRCS}) 
        EndIf(CXXSRCS MATCHES "^.*patsubst.*$") 
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
  
      If(line MATCHES "^DHDR.*$") 
        STRING(REGEX REPLACE "DHDR" "" DHDR "${line}")
        STRING(REGEX REPLACE ":=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[+]=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "=" "" DHDR "${DHDR}")
        STRING(REGEX REPLACE "[ ]+" ";" DHDR "${DHDR}")
        STRING(REGEX REPLACE "^;" "" DHDR "${DHDR}")
        SET(DHDR_PACKAGE ${DHDR_PACKAGE} ${DHDR}) 
      EndIf(line MATCHES "^DHDR.*$") 
  
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
            SET(PACKCXXFLAGS_PACKAGE ${PACKCXXFLAGS_PACKAGE} ${elem}) 
          endif(${elem} STREQUAL "${_match_result}")
        endforeach(elem ${PACKCXXFLAGS}) 
  #      MESSAGE("DEFINITIONS: ${PACKCXXFLAGS_PACKAGE}")
      EndIf(line MATCHES "^PACKCXXFLAGS.*$") 
  
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

      If(line MATCHES "^HDRS.*$")
          STRING(REGEX REPLACE "HDRS" "" HDRS "${line}")
    
          # Check if list of header files should be build from
          # list of source files. Add additional header files to
          # the list if they are defined. The parser does not 
          IF("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
            set(ADD_HEADERS TRUE)
          ELSE("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
             if(ADD_HEADERS)
               set(ADD_HEADERS TRUE)
             else(ADD_HEADERS)
               set(ADD_HEADERS FALSE)
             endif(ADD_HEADERS)
          ENDIF("${HDRS}" MATCHES "[$][(]SRCS:[.]cxx=[.]h[)]")
          STRING(REGEX REPLACE "[$][(]SRCS:[.]cxx=[.]h[)]" "" HDRS "${HDRS}")

#          STRING(REGEX REPLACE "[$][(]filter-out(.*),.*[)]" "//1" _exclude_h "${HDRS}")
#          STRING(LENGTH ${_exclude_h} _exclude_h_length)
#          If(${_exclude_h_length} GREATER 0)
#            String(STRIP ${_exclude_h} _exclude_h) 
#            list(REMOVE_ITEM HDRS ${_exclude_h})
#          EndIf(${_exclude_h_length} GREATER 0)
 
          STRING(REGEX REPLACE "[$][(]filter-out.*[)]" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE ":=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "[+]=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "=" "" HDRS "${HDRS}")
          STRING(REGEX REPLACE "[ ]+" ";" HDRS "${HDRS}")
          STRING(REGEX REPLACE "^;" "" HDRS "${HDRS}")
          SET(HDRS_PACKAGE ${HDRS_PACKAGE} ${HDRS}) 
      EndIf(line MATCHES "^HDRS.*$") 

#      MESSAGE("Parsed:  ${line}")
    Else(NOT EXCLUDE_PART)
#      MESSAGE("Not parsed:  ${line}")
    EndIf(NOT EXCLUDE_PART)

  EndForeach(line ${contents})
EndMacro(ParsePackageFile Package)

Macro(CreateMainCMakeFile path)
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
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_SUBDIR})
    file(APPEND  ${FileName} ")\n\n")
  endif(${definition_length} GREATER 0)

  file(APPEND ${FileName} "SetModule()\n\n")
  foreach(file ${PACKAGES})
    file(APPEND ${FileName} "include (CMake_${file}.txt)\n\n")
  endforeach(file ${PACKAGES})
EndMacro(CreateMainCMakeFile path)

Macro(CreateLibPackageFile path lib)

  set(FileName ${path}/CMake_${lib}.txt)
  set(AddCSource FALSE)
  set(AddFortranSource FALSE)

  STRING(REGEX REPLACE "^lib(.*)" "\\1" _lib "${lib}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  If(_file_glob)
    set(_counter 0)
    ForEach(_file ${_file_glob_dir})
      file(APPEND  ${FileName} "File(GLOB SRCS${_counter} RELATIVE \${CMAKE_CURRENT_SOURCE_DIR} \${CMAKE_SOURCE_DIR}/${_file}*.cxx)\n\n")
      Math(EXPR _counter ${_counter}+1)
    EndForEach(_file ${_file_glob_dir})
    Math(EXPR _counter ${_counter}-1)
    file(APPEND  ${FileName} "set(SRCS)\n")
    ForEach(_loop RANGE ${_counter})
      file(APPEND  ${FileName} "set(SRCS \${SRCS} \${SRCS${_loop}})\n")
    EndForEach(_loop RANGE ${_counter})
  Else(_file_glob)
    list(LENGTH CXXSRCS_PACKAGE cxx_length)
    if(${cxx_length} GREATER 0)
      file(APPEND  ${FileName} "set(SRCS\n")
      foreach(file ${CXXSRCS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${cxx_length} GREATER 0)
  EndIf(_file_glob)

  list(LENGTH CSRCS_PACKAGE c_length)
  if(${c_length} GREATER 0)
    set(AddCSource TRUE)
    file(APPEND  ${FileName} "set(CSRCS\n")
    foreach(file ${CSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${c_length} GREATER 0)

  list(LENGTH FSRCS_PACKAGE f_length)
  if(${f_length} GREATER 0)
    set(AddFortranSource TRUE)
    file(APPEND  ${FileName} "set(FSRCS\n")
    foreach(file ${FSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${FSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${f_length} GREATER 0)


  if(ADD_HEADERS)
    file(APPEND  ${FileName} "# fill list of header files from list of source files\n")
    file(APPEND  ${FileName} "# by exchanging the file extension\n")

    file(APPEND  ${FileName} "String(REPLACE \".cxx\" \".h\" HDRS \"\${SRCS}\")\n\n")
    list(LENGTH HDRS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS \${HDRS}\n")
      foreach(file ${HDRS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)    
  else(ADD_HEADERS)
    list(LENGTH HDRS_PACKAGE _length)
    if(${_length} GREATER 0)
      file(APPEND  ${FileName} "set(HDRS\n")
      foreach(file ${HDRS_PACKAGE})
        file(APPEND  ${FileName} "${file}\n")
      endforeach(file ${HDRS_PACKAGE})
      file(APPEND  ${FileName} ")\n\n")
    endif(${_length} GREATER 0)
  endif(ADD_HEADERS)

  If(AddCSource)
    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${CSRCS})\n\n")
  EndIf(AddCSource)
  If(AddFortranSource)
    file(APPEND  ${FileName} "Set(SRCS \${SRCS} \${FSRCS})\n\n")
  EndIf(AddFortranSource)

  list(LENGTH PACKFFLAGS_PACKAGE packf_length)
  if(${packf_length} GREATER 0)
    file(APPEND  ${FileName} "SET_SOURCE_FILES_PROPERTIES(\n")
    file(APPEND  ${FileName} "${FSRCS} PROPERTIES COMPILE_FLAGS\n") 
    file(APPEND  ${FileName} "${PACKFFLAGS_PACKAGE}\n")
    file(APPEND  ${FileName} ")\n\n")
  endif(${packf_length} GREATER 0)
  

  file(APPEND  ${FileName} "AddLibrary(${_lib} \"\${SRCS}\" \"\${HDRS}\")\n") 

  
  # If package file is empty remove the CMake input file for that packge and
  # remove the package from the list.
#  If( ${cxx_length} EQUAL 0 AND ${c_length} EQUAL 0 AND ${f_length} EQUAL 0)
#    file(REMOVE ${FileName})
#    list(REMOVE_ITEM PACKAGES lib${_lib})
#  EndIf( ${cxx_length} EQUAL 0 AND ${c_length} EQUAL 0 AND ${f_length} EQUAL 0)

EndMacro(CreateLibPackageFile path lib)

Macro(CreateBinPackageFile path bin)

  set(FileName ${path}/CMake_${bin}.txt)
  STRING(REGEX REPLACE "^bin(.*)" "\\1" _bin "${bin}")

  file(WRITE ${FileName} "# -*- mode: cmake -*-\n\n")

  file(APPEND  ${FileName} "set(SRCS)\n\n")

  list(LENGTH CXXSRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(SRCS\n")
    foreach(file ${CXXSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_PACKAGE})
    file(APPEND  ${FileName} ")\n\n")
  endif(${_length} GREATER 0)

  list(LENGTH CSRCS_PACKAGE _length)
  if(${_length} GREATER 0)
    file(APPEND  ${FileName} "set(CSRCS\n")
    foreach(file ${CSRCS_PACKAGE})
      file(APPEND  ${FileName} "${file}\n")
    endforeach(file ${CXXSRCS_PACKAGE})
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
EndMacro(CreateBinPackageFile path bin)
 
Macro(SpecialSettingsMadeByHand lib)
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
EndMacro(SpecialSettingsMadeByHand lib)
