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
rem * This script creates icepack.lib and icepack.dll from all .h and .cxx files
rem * in the current directory.
rem *
rem * In view of the ROOTCINT processing, the following two standard files
rem * are always required :
rem *
rem * ICEHeaders.h : containing an include of all .h files
rem * ICELinkDef.h : containing the #pragma's to define all classes
rem *
rem * --- NvE 04-apr-2000 UU-SAP Utrecht
rem ****************************************************************************

echo .
echo === Automatic ROOT library production of files ice.lib and ice.dll ===
echo .
rem --- Set the IcePack source directory as working directory
cd %ALIROOT%\RALICE\icepack

rem --- The option strings for MSVC++ DLL compilation and linking ---
set mscomp=/nologo /c /TP /Ze /MD /GR /GX /I%ROOTSYS%\include /I%ALIROOT%\RALICE
set msdll=/nologo /TP /Ze /MD /LD /GD /GR /GX /I%ROOTSYS%\include /I%ALIROOT%\RALICE
set mslink=/ENTRY:_DllMainCRTStartup@12 %ROOTSYS%\lib\*.lib %MYLIBS%\*.lib

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
echo * This script creates icepack.lib and icepack.dll from all .h and .cxx files
echo * in the current directory.
echo *
echo * In view of the ROOTCINT processing, the following two standard files
echo * are always required :
echo *
echo * ICEHeaders.h : containing an include of all .h files
echo * ICELinkDef.h : containing the #pragma's to define all classes
echo ****************************************************************************
goto end

:export
echo *** Creation of ROOT loadable export libraries
echo.
rem --- Creation of ROOT dictionary ---
rootcint zzzicepackdict.cxx -c -I%ALIROOT%\RALICE ICEHeaders.h ICELinkDef.h
rem --- Compilation step ---
cl %mscomp% *.cxx
rem --- Creation of the export LIB ---
bindexplib icepack *.obj > icepack.def
lib /nologo /machine:IX86 *.obj /def:icepack.def /out:icepack.lib
rem --- Creation of the DLL ---
link /nologo /machine:IX86 /DLL *.obj icepack.exp %mslink% /OUT:icepack.dll
rem --- Move the created libs to the SCRIPTS subdirectory
move icepack.lib .\scripts
move icepack.dll .\scripts
rem --- Delete intermediate files
del icepack.def
del icepack.exp
goto root_clean

:full
echo *** Creation of ROOT loadable full version libraries
echo.
rem --- Creation of ROOT dictionary ---
rootcint zzzicepackdict.cxx -c -I%ALIROOT%\RALICE ICEHeaders.h ICELinkDef.h
rem --- Creation of the DLL ---
cl %msdll% *.cxx /link %mslink% /OUT:icepack.dll
rem --- Creation of the full version LIB ---
lib /nologo /machine:IX86 *.obj /out:icepack.lib
rem --- Move the created libs to the SCRIPTS subdirectory
move icepack.lib .\scripts
move icepack.dll .\scripts
rem --- Delete intermediate files
goto root_clean

:root_clean
rem --- Delete all intermediate files --- 
del .def
del zzzicepackdict.h
del zzzicepackdict.cxx
del *.obj
echo.
echo *** mklibs done.
goto end

:end
rem --- Go back to original directory
cd scripts
rem --- End of script ---
