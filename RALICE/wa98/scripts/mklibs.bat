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
rem * This script creates rwa98.lib and rwa98.dll from all .h and .cxx files
rem * in the current directory.
rem *
rem * In view of the ROOTCINT processing, the following two standard files
rem * are always required :
rem *
rem * RWA98Headers.h : containing an include of all .h files
rem * RWA98LinkDef.h : containing the #pragma's to define all classes
rem *
rem * --- NvE 04-apr-2000 UU-SAP Utrecht
rem ****************************************************************************

echo .
echo === Automatic ROOT library production of files rwa98.lib and rwa98.dll ===
echo .

set alice=c:\nick\cxx\source\alice\AliRoot\RALICE

rem --- The option strings for MSVC++ DLL compilation and linking ---
set mscomp=/nologo /c /TP /Ze /MD /GR /GX /I%ROOTSYS%\include /I%alice%
set msdll=/nologo /TP /Ze /MD /LD /GD /GR /GX /I%ROOTSYS%\include /I%alice%
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
echo * This script creates rwa98.lib and rwa98.dll from all .h and .cxx files
echo * in the current directory.
echo *
echo * In view of the ROOTCINT processing, the following two standard files
echo * are always required :
echo *
echo * RWA98Headers.h : containing an include of all .h files
echo * RWA98LinkDef.h : containing the #pragma's to define all classes
echo ****************************************************************************
goto end

:export
echo *** Creation of ROOT loadable export libraries
echo.
rem --- Creation of ROOT dictionary ---
rootcint zzzrwa98dict.cxx -c -Ic:/nick/cxx/source/alice/AliRoot/RALICE RWA98Headers.h RWA98LinkDef.h
rem --- Compilation step ---
cl %mscomp% *.cxx
rem --- Creation of the export LIB ---
bindexplib rwa98 *.obj > rwa98.def
lib /nologo /machine:IX86 *.obj /def:rwa98.def /out:rwa98.lib
rem --- Creation of the DLL ---
link /nologo /machine:IX86 /DLL *.obj rwa98.exp %mslink% /OUT:rwa98.dll
del rwa98.def
del rwa98.exp
goto root_clean

:full
echo *** Creation of ROOT loadable full version libraries
echo.
rem --- Creation of ROOT dictionary ---
rootcint zzzrwa98dict.cxx -c -Ic:/nick/cxx/source/alice/AliRoot/RALICE RWA98Headers.h RWA98LinkDef.h
rem --- Creation of the DLL ---
cl %msdll% *.cxx /link %mslink% /OUT:rwa98.dll
rem --- Creation of the full version LIB ---
lib /nologo /machine:IX86 *.obj /out:rwa98.lib
goto root_clean

:root_clean
rem --- Delete all intermediate files --- 
del .def
del zzzrwa98dict.h
del zzzrwa98dict.cxx
del *.obj
echo.
echo *** mklibs done.
goto end

:end
rem --- End of script ---
