#!/bin/sh
### Shell script to create a ROOT loadable ALPHA-CXX shared lib out of .cxx source code
###
### NvE 05-jul-2000 UU-SAP Utrecht

### Name of the produced shared library
lib=rwa98.so

### The ALICE source directory
alice=$HOME/cxx/source/alice

### The option string for ALPHA-CXX shared lib compilation and linking ***
### For the ALPHA-CXX ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
alpharoot="-x cxx -g0 -shared -w1 -I$ROOTSYS/include -I$alice -o $lib"

### Create the dictionary files
rootcint zzzrwa98dict.cxx -c -I$alice RWA98Headers.h RWA98LinkDef.h

### Compile and create the ROOT loadable shared library
cxx $alpharoot *.cxx   

rm zzzrwa98dict.*
rm *.o
rm so_locations

echo ' '
echo '*** alphacxxlib done. Result in ' $lib
