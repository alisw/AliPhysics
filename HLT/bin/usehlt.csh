#!/bin/tcsh
###########
# $Id$

if ( "a$1" == "a" )  then
 setenv ALIHLT_USEPACKAGE ALIROOT
else
 setenv ALIHLT_USEPACKAGE $1
endif
echo HLT ALIHLT_USEPACKAGE=$ALIHLT_USEPACKAGE

setenv ALIHLT_LIBDIR ${ALIHLT_TOPDIR}/lib_${ALIHLT_USEPACKAGE}

#setenv ALIHLT_NOLOGGING false
#setenv ALIHLT_DOMC true
#setenv ALIHLT_ALIDETECT true
