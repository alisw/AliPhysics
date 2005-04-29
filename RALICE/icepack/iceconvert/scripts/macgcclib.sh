#!/bin/sh
### Shell script to create a ROOT loadable GCC shared lib out of .cxx source code
###
### NvE 23-may-2000 UU-SAP Utrecht
# 
### Name of the produced shared libraries
lib1=iceconvert.so
lib2=iceconvert.dylib
#
### Some MAC specific settings
export MACOSX_DEPLOYMENT_TARGET=10.3
unset LD_PREBIND
# 
### The option string for GCC shared lib compilation and linking ***
### For the GCC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
gcccomp="-c -g0 -Wall -Wno-long-long -I$ROOTSYS/include -I$ALIROOT/RALICE -I$ALIROOT/RALICE/icepack"
gccroot="-c -g0 -ansi -pedantic -Wall -Wno-long-long -I$ROOTSYS/include -I$ALIROOT/RALICE -I$ALIROOT/RALICE/icepack"
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
#
# Compilation phase
gcc $gcccomp *.c   
g++ $gccroot *.cxx   
#
# Creating ralice.so library 
g++ -v -bundle -undefined dynamic_lookup -o $lib1 *.o
#
# Creating ralice.dylib library 
g++ -v -dynamiclib -undefined dynamic_lookup -single_module -o $lib2 *.o
### On some systems the following extra "-read_only_relocs" flag might be needed
# g++ -v -dynamiclib -undefined dynamic_lookup -single_module -read_only_relocs -o $lib2 *.o
#
rm zzziceconvertdict.*
rm *.o
# 
### Move the created lib to the scripts directory and go there
mv $lib1 scripts
mv $lib2 scripts
cd scripts
#
echo ' ' 
echo '*** macgcclib.sh done. Results in iceconvert.so and iceconvert.dylib' 
