classfilename=$1

classfilename_h=$classfilename.h
classfilename_cxx=$classfilename.cxx
classname=$classfilename

if [  -a $classfilename_h ] || [ -a $classfilename_cxx ];then
    echo ERROR, $classfilename_h and $classfilename_cxx allready exist, delelte them or rename the class you want to make
else 
    echo creating new files   $classfilename_h and $classfilename_cxx
    printf "#ifndef " > $classfilename_h
    classguard=${classfilename_h/.h/_H}

    echo $classguard | tr "[:lower:]" "[:upper:]"  >> $classfilename_h
    printf "#define " >> $classfilename_h
    echo $classguard | tr "[:lower:]" "[:upper:]"  >> $classfilename_h

    printf "\n" >>$classfilename_h 

    printf "/**************************************************************************\n"  >> $classfilename_h
    printf " * This file is property of and copyright by the Experimental Nuclear     *\n"  >> $classfilename_h
    printf " * Physics Group, Dep. of Physics                                         *\n"  >> $classfilename_h
    printf " * University of Oslo, Norway, 2007                                       *\n"  >> $classfilename_h
    printf " *                                                                        *\n"  >> $classfilename_h
    printf " * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*\n"  >> $classfilename_h
    printf " * Contributors are mentioned in the code where appropriate.              *\n"  >> $classfilename_h
    printf " * Please report bugs to perthi@fys.uio.no                                *\n"  >> $classfilename_h 
    printf " *                                                                        *\n"  >> $classfilename_h
    printf " * Permission to use, copy, modify and distribute this software and its   *\n"  >> $classfilename_h
    printf " * documentation strictly for non-commercial purposes is hereby granted   *\n"  >> $classfilename_h
    printf " * without fee, provided that the above copyright notice appears in all   *\n"  >> $classfilename_h
    printf " * copies and that both the copyright notice and this permission notice   *\n"  >> $classfilename_h
    printf " * appear in the supporting documentation. The authors make no claims     *\n"  >> $classfilename_h
    printf " * about the suitability of this software for any purpose. It is          *\n"  >> $classfilename_h
    printf " * provided \"as is\" without express or implied warranty.                  *\n"  >> $classfilename_h
    printf " **************************************************************************/\n"  >> $classfilename_h

    printf "\n\n\n" >>$classfilename_h 
    printf "class " >> $classfilename_h
    printf " $classname\n{\n\t$classname();\n\tvirtual ~$classname();\n};\n\n#endif\n" >> $classfilename_h  

    printf "/**************************************************************************\n"  >  $classfilename_cxx
    printf " * This file is property of and copyright by the Experimental Nuclear     *\n"  >> $classfilename_cxx
    printf " * Physics Group, Dep. of Physics                                         *\n"  >> $classfilename_cxx
    printf " * University of Oslo, Norway, 2007                                       *\n"  >> $classfilename_cxx
    printf " *                                                                        *\n"  >> $classfilename_cxx
    printf " * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*\n"  >> $classfilename_cxx
    printf " * Contributors are mentioned in the code where appropriate.              *\n"  >> $classfilename_cxx
    printf " * Please report bugs to perthi@fys.uio.no                                *\n"  >> $classfilename_cxx 
    printf " *                                                                        *\n"  >> $classfilename_cxx
    printf " * Permission to use, copy, modify and distribute this software and its   *\n"  >> $classfilename_cxx
    printf " * documentation strictly for non-commercial purposes is hereby granted   *\n"  >> $classfilename_cxx
    printf " * without fee, provided that the above copyright notice appears in all   *\n"  >> $classfilename_cxx
    printf " * copies and that both the copyright notice and this permission notice   *\n"  >> $classfilename_cxx
    printf " * appear in the supporting documentation. The authors make no claims     *\n"  >> $classfilename_cxx
    printf " * about the suitability of this software for any purpose. It is          *\n"  >> $classfilename_cxx
    printf " * provided \"as is\" without express or implied warranty.                  *\n"  >> $classfilename_cxx
    printf " **************************************************************************/\n"  >> $classfilename_cxx
    printf "#include \"$classfilename_h\"\n\n" >>  $classfilename_cxx
    printf "$classname::$classname()\n{\n\n}\n\n" >> $classfilename_cxx
    printf "$classname::~$classname()\n{\n\n}\n\n" >> $classfilename_cxx

    emacs $classfilename_h  &
    emacs $classfilename_cxx  &

fi

