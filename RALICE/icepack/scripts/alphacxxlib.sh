#!/bin/sh
### Shell script to create a ROOT loadable ALPHA-CXX shared lib out of .cxx source code
###
### NvE 05-jul-2000 UU-SAP Utrecht

### Name of the produced shared library
lib=icepack.so

### The option string for ALPHA-CXX shared lib compilation and linking ***
### For the ALPHA-CXX ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
alpharoot="-x cxx -g0 -shared -w1 -I$ROOTSYS/include -I$ALIROOT/RALICE -o $lib"

### Go to the directory with the source files
cd $ALIROOT/RALICE/icepack

### Create the dictionary files
rootcint zzzicepackdict.cxx -c -I$ALIROOT/RALICE ICEHeaders.h ICELinkDef.h

### Compile and create the ROOT loadable shared library
cxx $alpharoot *.cxx   

rm zzzicepackdict.*
rm *.o
rm so_locations

### Move the created lib to the scripts directory and go there
mv $lib scripts
cd scripts

echo ' '
echo '*** alphacxxlib done. Result in ' $lib
