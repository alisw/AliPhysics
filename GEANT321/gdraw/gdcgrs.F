*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:21  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDCGRS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *     Reset the pointers and counters for the Zebra Bank         *
C.    *     of CG, HS and Wire Objects.                                *
C.    *                                                                *
C.    *    ==>Called by : GDRAW,GDCGHI                                 *
C.    *                                                                *
C.    *                                                                *
C.    *       Author : J. Salt ; S. Giani  **********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcgobj.inc"
*
      IST=2
*      DO 10 I=1,NST
*         LPST(I)=0
* 10   CONTINUE
*      LPST(1)=1
*      LPST(2)=4000
*      LPST(3)=8000
      IFCG=3
      ILCG=2
*SG
*      NTCUR=0
*      NFILT=0
*      NTNEX=0
*SG
      KCGST=0
      NCGVOL=0
      IVCLOS=0
      IVFUN=1
      NCLAS1=0
*SG
*      NCLAS2=0
*      NCLAS3=0
*SG
      END
 
