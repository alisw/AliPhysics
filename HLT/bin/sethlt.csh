#!/bin/csh
###########
# $Id$

setenv CVS_RSH ssh
setenv CVSROOT ${USER}@kjekspc1.fi.uib.no:/cvs/hltcvs

setenv ALIHLT_USEPACKAGE ALIROOT
#setenv ALIHLT_USEPACKAGE ROOT
#setenv ALIHLT_USEPACKAGE STANDALONE

setenv ALIHLT_BASEDIR $HOME/work/hlt
setenv ALIHLT_TOPDIR ${ALIHLT_BASEDIR}/level3code
setenv ALIHLT_LIBDIR ${ALIHLT_TOPDIR}/lib_${ALIHLT_USEPACKAGE}

setenv ALIHLT_NOLOGGING false
setenv ALIHLT_DOMC true
setenv ALIHLT_ALIDETECT true
setenv ALIHLT_ROWHOUGH false
setenv ALIHLT_MLUCDIR ${ALIHLT_BASEDIR}/kip/MLUC

setenv ALIHLT_DATADIR /data1/head
#setenv ALIHLT_TRANSFORMFILE ${ALIHLT_DATADIR}/l3transform.config
#setenv ALIHLT_GEOPATH ${ALIHLT_DATADIR}

setenv LD_LIBRARY_PATH ${ROOTSYS}/lib\:${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}\:${ALIHLT_MLUCDIR}/lib

