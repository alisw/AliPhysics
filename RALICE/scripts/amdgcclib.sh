#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared library
lib=ralice.so
# 
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gccroot="-fPIC -shared -g0 -ansi -pedantic -Wall -Wno-long-long -Woverloaded-virtual -I$ROOTSYS/include -o $lib"
#
echo "lib = " $lib
echo "gccroot = " $gccroot 
#
### Go to the directory with the source files
cd ..
#
### Create the dictionary files
rootcint -f zzzralicedict.cxx -c RALICEHeaders.h RALICELinkDef.h
# 
### Compile and create the ROOT loadable shared library
g++ $gccroot *.cxx   
# 
rm zzzralicedict.*
rm *.o
# 
### Move the created lib to the scripts directory and go there
mv $lib scripts
cd scripts
#
echo ' ' 
echo '*** gcclib done. Result in ' $lib 
