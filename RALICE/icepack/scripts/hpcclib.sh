#!/bin/sh
### Shell script to create a ROOT loadable HP-CC shared lib out of .cxx source code
###
### NvE 28-jun-1999 UU-SAP Utrecht
#
### The option strings for HP-CC shared lib compilation and linking ***
### For the HP-CC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***
hpcomp="-c -s -z +z +a1 +w +DAportable -I$ROOTSYS/include -I$ALIROOT/RALICE"
hproot="-c -s -z +z +a1 +DAportable -I$ROOTSYS/include -I$ALIROOT/RALICE"
hplink="-L$ROOTSYS/lib/ -l*.sl -lm"

### Go to the directory with the source files
cd $ALIROOT/RALICE/icepack

rootcint zzzicepackdict.cxx -c -I$ALIROOT/RALICE ICEHeaders.h ICELinkDef.h

CC $hproot *.cxx   

CC -b -o icepack.sl *.o

rm zzzicepackdict.*
rm *.o

### Move the created lib to the scripts directory and go there
mv icepack.sl scripts
cd scripts

echo '*** hpcclib.sh done. Result in icepack.sl'
