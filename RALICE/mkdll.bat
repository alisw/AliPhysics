@echo off
rem ****************************************************************************
rem *       Script to create a relocatable MSVC++ DLL from *.cxx files
rem *
rem * Usage : mkdll
rem *
rem * This will create ralice.dll from all .h and .cxx files in the current dir
rem *
rem * In view of the ROOTCINT processing, the following two standard files
rem * are always required :
rem *
rem * allhead.h ==> containing an include of all .h files
rem * linkdef.h ==> containing the #pragma's to define all classes
rem *
rem ****************************************************************************
rem *
echo .
echo === Automatic ROOT DLL production of file ralice.dll ===
echo .
rem *
rem --- The option strings for MSVC++ DLL compilation and linking ***
set mscomp=/nologo /c /TP /Za /MD /I%ROOTSYS%\include
set msdll=/nologo /TP /Za /MD /LD /GD /I%ROOTSYS%\include
set mslink=/ENTRY:_DllMainCRTStartup@12 %ROOTSYS%\lib\*.lib %MYLIBS%\*.lib
rem *
rootcint zzzralicedict.cxx -c allhead.h linkdef.h
rem *
cl %msdll% *.cxx /link %mslink% /OUT:ralice.dll
rem *
rem --- Delete all intermediate files 
del .def
del containing
del zzzralicedict.h
del zzzralicedict.cxx
del *.obj
rem *
echo *** mkdll done.
