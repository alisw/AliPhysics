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
*
      SUBROUTINE GDCGPR
C.    ******************************************************************
C.    *                                                                *
C.    *    The aim of this subroutine is:                              *
C.    *                                                                *
C.    *    - To control the number of structures used                  *
C.    *    - To test if the maximum size of memory is achieved         *
C.    *                                                                *
C.    *    - To increase the size of the CG-Bank if there              *
C.    *    is not enough space for new CG structures.                  *
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GDCGCL,GDCGHI,GDCGOB                         *
C.    *                                                                *
C.    *       Author : J. Salt   **********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcgobj.inc"
#include "geant321/gchiln.inc"
*
      INCND=2500
*
*  Control on the Maximum Number of CG Structures
*
      IF(IST.GT.IQ(JCOUNT-1))THEN
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
         KCGST=-1
         RETURN
      ENDIF
*
*   Control the Size of the IXCG
*   Bank Extension. Garbage Collection
*
      NSUM=NFILT+NTNEX
      NWFIL=IQ(JCG-1)
      IF(NSUM.GT.NWFIL)THEN
         NRES=NSUM-IQ(JCG-1)
         IF(NRES.GT.INCND)INCND=NRES
**
         IXCG=IXSTOR+1
         CALL MZPUSH(IXCG,JCG,0,INCND,' ')
      ENDIF
      IXCG=IXSTOR+1
      CALL MZGARB(IXCG,0)
*
10000 FORMAT('  Maximum Number of Structures exceeded')
*
      END
 
