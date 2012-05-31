#!/bin/sh
# $Id$

# Script for installation of cuts.cxx program

CURDIR=`pwd`

echo "... Installing cuts program"

g++ -I$ROOTSYS/include \
    `root-config --glibs` -lGeomPainter -lGeom cuts.cxx \
    -o cuts2

echo "... Installing cuts program finished"

cd $CURDIR
