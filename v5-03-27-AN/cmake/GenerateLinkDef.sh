#!/bin/sh

# 
# Script to generate LinkDef file from a list of classes. 
#
# Usage: GenerateLinkDef.sh [CLASS ...]
#
# Author: Christian Holm Christensen <cholm@nbi.dk>
#

# --- Help message ---------------------------------------------------
usage()
{
    cat<<EOF
Usage: [OPTIONS] $0 [CLASS ...] 

Options:
	-h,--help		This help 
	-S,--no-evolution	Do not enable schema evolution 
	-C,--custom-streamer	Use custom streamers 
EOF
}

# --- Default options ------------------------------------------------
schema="+"
custom=""

# --- Process options ------------------------------------------------
while test $# -gt 0 ; do 
    case $1 in 
	-h|--help) 		usage ; exit 0 ;; 
	-S|--no-evolution)	schema="" ;;
	-C|--custom-streamer)	custom="-" ;; 
	*) break ;; 
    esac
    shift 
done 

# --- Header ---------------------------------------------------------
cat <<EOF
// Auto generated file - do not edit
#ifndef __CINT__
# error Not for compilation
#else 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

EOF

# --- Link definitions -----------------------------------------------
while test $# -gt 0 ; do 
    echo "#pragma link C++ class ${1}${schema}${custom};" 
    shift 
done 

# --- Trailer --------------------------------------------------------
cat <<EOF

#endif // __CINT__
//
// EOF
//
EOF

#
# EOF
#



