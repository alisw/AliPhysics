#ifndef CERNLIB_ISAJET_ZLINKA_INC
#define CERNLIB_ISAJET_ZLINKA_INC
*
*
* zlinka.inc
*
C
C              general Zebra link area
C   use with utility subroutines GSLINK,GRLINK,RSLINK and RRLINK
      INTEGER NSLINK,NRLINK,LSLINK,LRLINK
      PARAMETER (NSLINK=100)
      PARAMETER (NRLINK=100)
      COMMON/ZLINKA/LSLINK(NSLINK),LRLINK(NRLINK)
      SAVE /ZLINKA/

#endif
