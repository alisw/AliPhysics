#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared library
lib=rwa98.so
#
### The ALICE source directory
alice=$HOME/cxx/source/alice
# 
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gccroot="-shared -g0 -ansi -pedantic -Wall -Wno-long-long -I$ROOTSYS/include -I$alice -o $lib"
#
echo "lib = " $lib
echo "alice = " $alice
echo "gccroot = " $gccroot 
#
### Create the dictionary files
rootcint -f zzzrwa98dict.cxx -c -I$alice RWA98Headers.h RWA98LinkDef.h
# 
### Compile and create the ROOT loadable shared library
g++ $gccroot *.cxx   
# 
rm zzzwa98dict.*
rm *.o
# 
echo ' ' 
echo '*** gcclib done. Result in ralice.so' 
