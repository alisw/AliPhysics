@echo off
rem ****************************************************************************
rem *   Script to create an MSVC++ LIB and relocatable DLL from *.cxx files
rem *
rem * Usage :
rem * -------
rem * mklibs export : ROOT loadable DLL and export LIB are created
rem * mklibs full   : ROOT loadable DLL and full LIB version are created
rem * 
rem * Notes :
rem * -------
rem * 1) "mklibs export" is the default, enabling ROOT loadable library creation
rem *    via 'double clicking'.
rem *
rem * 2) Providing unsupported options results in displaying the help info.  
rem *
rem * This script creates ralice.lib and ralice.dll from all .h and .cxx files
rem * in the current directory.
rem *
rem * In view of the ROOTCINT processing, the following two standard files
rem * are always required :
rem *
rem * RALICEHeaders.h : containing an include of all .h files
rem * RALICELinkDef.h : containing the #pragma's to define all classes
rem *
rem * --- NvE 04-apr-2000 UU-SAP Utrecht
rem ****************************************************************************

echo .
echo === Automatic ROOT library production of files ralice.lib and ralice.dll ===
echo .

rem --- The option strings for MSVC++ DLL compilation and linking ---
set mscomp=/nologo /c /Ze /TP /MD /GR /GX /I%ROOTSYS%\include
set msdll=/nologo /Ze /TP /MD /LD /GD /GR /GX /I%ROOTSYS%\include
set mslink=/ENTRY:_DllMainCRTStartup@12 %ROOTSYS%\lib\*.lib

if "%1" == "" goto export
if "%1" == "export" goto export
if "%1" == "full" goto full

rem --- Displaying of the help info ---
echo ****************************************************************************
echo *   Script to create an MSVC++ LIB and relocatable DLL from *.cxx files
echo *
echo * Usage :
echo * -------
echo * mklibs export : ROOT loadable DLL and export LIB are created
echo * mklibs full   : ROOT loadable DLL and full LIB version are created
echo * 
echo * Notes :
echo * -------
echo * 1) "mklibs export" is the default, enabling ROOT loadable library creation
echo *    via 'double clicking'.
echo * 2) Providing unsupported options results in displaying the help info.  
echo *
echo * This script creates ralice.lib and ralice.dll from all .h and .cxx files
echo * in the current directory.
echo *
echo * In view of the ROOTCINT processing, the following two standard files
echo * are always required :
echo *
echo * RALICEHeaders.h : containing an include of all .h files
echo * RALICELinkDef.h : containing the #pragma's to define all classes
echo ****************************************************************************
goto end

:export
echo *** Creation of ROOT loadable export libraries
echo.
rem --- Set the RALICE source directory as working directory
cd ..
rem --- Creation of ROOT dictionary ---
rootcint zzzralicedict.cxx -c RALICEHeaders.h RALICELinkDef.h
rem --- Compilation step ---
cl %mscomp% *.cxx
rem --- Creation of the export LIB ---
bindexplib ralice *.obj > ralice.def
lib /nologo /machine:IX86 *.obj /def:ralice.def /out:ralice.lib
rem --- Creation of the DLL ---
link /nologo /machine:IX86 /DLL *.obj ralice.exp %mslink% /OUT:ralice.dll
rem --- Move the created libs to the SCRIPTS subdirectory
move ralice.lib .\scripts
move ralice.dll .\scripts
rem --- Delete all intermediate files --- 
del .def
del ralice.def
del ralice.exp
del zzzralicedict.h
del zzzralicedict.cxx
del *.obj
echo.
echo *** mklibs done.
goto end

:full
echo *** Creation of ROOT loadable full version libraries
echo.
rem --- Set the RALICE source directory as working directory
cd ..
rem --- Creation of ROOT dictionary ---
rootcint zzzralicedict.cxx -c RALICEHeaders.h RALICELinkDef.h
rem --- Creation of the DLL ---
cl %msdll% *.cxx /link %mslink% /OUT:ralice.dll
rem --- Creation of the full version LIB ---
lib /nologo /machine:IX86 *.obj /out:ralice.lib
rem --- Move the created libs to the SCRIPTS subdirectory
move ralice.lib .\scripts
move ralice.dll .\scripts
rem --- Delete all intermediate files --- 
del .def
del zzzralicedict.h
del zzzralicedict.cxx
del *.obj
echo.
echo *** mklibs done.
goto end

:end
rem --- Go back to original directory
cd scripts
rem --- End of script ---
