#!/bin/bash
#
# Shell script to do all the LEGO plots 
#
for i in Inner ITS PIPE FMD Nothing ; do 
    aliroot -l -b -q FMD/scripts/MakeLego.C\(\"$i\"\)
done 

root -l -q FMD/scripts/DrawLego.C
