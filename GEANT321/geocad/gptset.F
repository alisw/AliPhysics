*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:46  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE GPTSET (IVOLNA, NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print tracking media data structure JTMED     *
C.    *       NUMB     medium number                                   *
C.    *                                                                *
C.    *    ==>Called by : GDEXCA    PRINTS the media of the volumes    *
C.    *       exported to SET in the material file .mat.               *
C.    *       Author: Nils Hoimyr                                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcsetf.inc"
*
      CHARACTER CHVOLN*4, CHNATM*20, CHNAMA*20
*
      IF (JTMED.LE.0) RETURN
      JTM = LQ(JTMED-NUMB)
      IF (JTM.NE.0) THEN
         CALL UHTOC(IVOLNA,4,CHVOLN,4)
         CALL UHTOC(IQ(JTM+1),4,CHNATM,20)
         IMAT = Q(JTM+6)
         JMA = LQ(JMATE-IMAT)
         CALL UHTOC(IQ(JMA+1),4,CHNAMA,20)
         DENS = Q(JMA+8)
         WRITE (NUNIT2,10000) CHVOLN,NUMB,CHNATM,IMAT,CHNAMA,DENS
      ENDIF
*
10000 FORMAT(1X,A4,6X,I6,1X,A20,I4,1X,A20,E15.8)
      END
 
