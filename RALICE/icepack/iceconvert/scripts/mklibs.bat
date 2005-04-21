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
rem * This script creates iceconvert.lib and iceconvert.dll from all .h, .c
rem * and .cxx files in the source directory.
rem *
rem * In view of the ROOTCINT processing, the following two standard files
rem * are always required :
rem *
rem * ICEConvHeaders.h : containing an include of all .h files
rem * ICEConvLinkDef.h : containing the #pragma's to define all classes
rem *
rem * --- NvE 11-mar-2005 Utrecht University
rem ****************************************************************************

echo .
echo === Automatic ROOT library production of files iceconvert.lib and iceconvert.dll ===
echo .
rem --- Set the IceConvert source directory as working directory
cd %ALIROOT%\RALICE\icepack\iceconvert

rem --- The option strings for MSVC++ DLL compilation and linking ---
set msc=/nologo /c /Ze /MD /GR /GX
set mscomp=/nologo /c /TP /Ze /MD /GR /GX /I%ROOTSYS%\include /I%ALIROOT%\RALICE /I%ALIROOT%\RALICE\icepack
set msdll=/nologo /Ze /MD /LD /GD /GR /GX /I%ROOTSYS%\include /I%ALIROOT%\RALICE /I%ALIROOT%\RALICE\icepack
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
echo * This script creates iceconvert.lib and iceconvert.dll from all .h, .c
echo * and .cxx files in the source directory.
echo *
echo * In view of the ROOTCINT processing, the following two standard files
echo * are always required :
echo *
echo * ICEConvHeaders.h : containing an include of all .h files
echo * ICEConvLinkDef.h : containing the #pragma's to define all classes
echo ****************************************************************************
goto end

:export
echo *** Creation of ROOT loadable export libraries
echo.
cl %msc% *.c
rem --- Creation of ROOT dictionary ---
rootcint zzziceconvertdict.cxx -c -p -I%ALIROOT%\RALICE -I%ALIROOT%\RALICE\icepack ICEConvHeaders.h ICEConvLinkDef.h
rem --- Compilation step ---
cl %mscomp% *.cxx
rem --- Creation of the export LIB ---
bindexplib iceconvert *.obj > iceconvert.def
lib /nologo /machine:IX86 *.obj /def:iceconvert.def /out:iceconvert.lib
rem --- Creation of the DLL ---
link /nologo /machine:IX86 /DLL *.obj iceconvert.exp %mslink% /OUT:iceconvert.dll
rem --- Move the created libs to the SCRIPTS subdirectory
move iceconvert.lib .\scripts
move iceconvert.dll .\scripts
rem --- Delete intermediate files
del iceconvert.def
del iceconvert.exp
goto root_clean

:full
echo *** Creation of ROOT loadable full version libraries
echo.
rem --- Creation of ROOT dictionary ---
rootcint zzziceconvertdict.cxx -c -I%ALIROOT%\RALICE -I%ALIROOT%\RALICE\icepack ICEConvHeaders.h ICEConvLinkDef.h
rem --- Creation of the DLL ---
cl %msdll% *.c *.cxx /link %mslink% /OUT:iceconvert.dll
rem --- Creation of the full version LIB ---
lib /nologo /machine:IX86 *.obj /out:iceconvert.lib
rem --- Move the created libs to the SCRIPTS subdirectory
move iceconvert.lib .\scripts
move iceconvert.dll .\scripts
rem --- Delete intermediate files
goto root_clean

:root_clean
rem --- Delete all intermediate files --- 
del .def
del zzziceconvertdict.h
del zzziceconvertdict.cxx
del *.obj
echo.
echo *** mklibs done.
goto end

:end
rem --- Go back to original directory
cd scripts
rem --- End of script ---
