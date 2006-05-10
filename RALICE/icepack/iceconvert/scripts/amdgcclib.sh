#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared library
lib=iceconvert.so
#
### The option string for GCC compilation of the .c code ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gcccomp="-fPIC -c -g0 -Wall -Wno-long-long -I$ROOTSYS/include -I$ALIROOT/RALICE -I$ALIROOT/RALICE/icepack"
#
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gccroot="-fPIC -shared -g0 -ansi -pedantic -Wall -Wno-long-long -Woverloaded-virtual -I$ROOTSYS/include -I$ALIROOT/RALICE -I$ALIROOT/RALICE/icepack -o $lib"
#
echo "lib = " $lib
echo "gcccomp = " $gcccomp 
echo "gccroot = " $gccroot 
#
### Go to the directory with the source files
cd $ALIROOT/RALICE/icepack/iceconvert
#
### Create the dictionary files
rootcint -f zzziceconvertdict.cxx -c -p -I$ALIROOT/RALICE -I$ALIROOT/RALICE/icepack ICEConvHeaders.h ICEConvLinkDef.h
# 
### Compile and create the ROOT loadable shared library
gcc $gcccomp *.c   
g++ $gccroot *.cxx *.o   
# 
rm zzziceconvertdict.*
rm *.o
# 
### Move the created lib to the scripts directory and go there
mv $lib scripts
cd scripts
#
echo ' ' 
echo '*** gcclib done. Result in ' $lib 
