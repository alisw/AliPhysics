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

### Go to the directory with the source files
cd ..

rootcint zzzralicedict.cxx -c RALICEHeaders.h RALICELinkDef.h

CC $hproot *.cxx   

CC -b -o ralice.sl *.o

rm zzzralicedict.*
rm *.o

### Move the created lib to the scripts directory and go there
mv ralice.sl scripts
cd scripts

echo '*** hpcclib done. Result in ralice.sl'
