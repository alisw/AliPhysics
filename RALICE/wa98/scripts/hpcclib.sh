#!/bin/sh
### Shell script to create a ROOT loadable HP-CC shared lib out of .cxx source code
###
### NvE 28-jun-1999 UU-SAP Utrecht
#
### The option strings for HP-CC shared lib compilation and linking ***
### For the HP-CC ROOT loadable shared lib the strict requirements are ***
### dropped to avoid many warnings from the rootcint generated code ***

hpcomp="-c -s -z +z +a1 +w +DAportable -I$ROOTSYS/include"
hproot="-c -s -z +z +a1 +DAportable -I$ROOTSYS/include"
hplink="-L$ROOTSYS/lib/ -l*.sl -lm"

rootcint zzzrwa98dict.cxx -c RWA98Headers.h RWA98LinkDef.h

CC $hproot *.cxx   

CC -b -o rwa98.sl *.o

rm zzzrwa98dict.*
rm *.o

echo '*** hpcclib done. Result in rwa98.sl'
