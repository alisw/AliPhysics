#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared library
lib=rwa98.so
# 
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gccroot="-fPIC -shared -g0 -ansi -pedantic -Wall -Wno-long-long -I$ROOTSYS/include -I$ALIROOT/RALICE -o $lib"
#
echo "lib = " $lib
echo "gccroot = " $gccroot 
#
### Go to the directory with the source files
cd $ALIROOT/RALICE/wa98
#
### Create the dictionary files
rootcint -f zzzrwa98dict.cxx -c -I$ALIROOT/RALICE RWA98Headers.h RWA98LinkDef.h
# 
### Compile and create the ROOT loadable shared library
g++ $gccroot *.cxx   
# 
rm zzzrwa98dict.*
rm *.o
# 
### Move the created lib to the scripts directory and go there
mv $lib scripts
cd scripts
#
echo ' ' 
echo '*** gcclib done. Result in ' $lib 
