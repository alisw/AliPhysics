*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GSTPAR(ITMED,CHPAR,PARVAL)
*
************************************************************************
*                                                                      *
*  To change the value of cut  or mechanism "CHPAR"                    *
*      to a new value PARVAL  for tracking medium ITMED                *
*    The  data   structure  JTMED   contains  the   standard  tracking *
*  parameters (CUTS and flags to control the physics processes)  which *
*  are used  by default  for all  tracking media.   It is  possible to *
*  redefine individually  with GSTPAR  any of  these parameters  for a *
*  given tracking medium.                                              *
*  ITMED     tracking medium number                                    *
*  CHPAR     is a character string (variable name)                     *
*  PARVAL    must be given as a floating point.                        *
*  For     example     to     change    CUTGAM     to     0.0001       *
*                                                                      *
*    ==>Called by : <USER>                                             *
*       Author    R.Brun     *********                                 *
*                                                                      *
************************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      DIMENSION CUTS(10),MECA(5,13)
      EQUIVALENCE (CUTS(1),CUTGAM),(MECA(1,1),IPAIR)
      CHARACTER*(*) CHPAR
C.
C.    ------------------------------------------------------------------
C.
      IF(ITMED.LE.0)GO TO 90
      IF(ITMED.GT.NTMED)GO TO 90
      JTM=LQ(JTMED-ITMED)
      IF(JTM.LE.0)GO TO 90
      JTMN=LQ(JTM)
      IF(JTMN.EQ.0)THEN
         CALL MZBOOK(IXCONS,JTMN,JTM,0,'TCUT',0,0,40,3,0)
         IQ(JTMN-5)=ITMED
         DO 10 I=1,10
            Q(JTMN+I)=CUTS(I)
  10     CONTINUE
         DO 20 I=1,13
            Q(JTMN+10+I)=MECA(1,I)
  20     CONTINUE
*
* *** New mechanisms in version 3.16
         Q(JTMN+10+21)=ILABS
         Q(JTMN+10+22)=ISYNC
         Q(JTMN+10+23)=ISTRA
      ENDIF
C
      ITPAR=0
      IF(CHPAR.EQ.'CUTGAM')ITPAR=1
      IF(CHPAR.EQ.'CUTELE')ITPAR=2
      IF(CHPAR.EQ.'CUTNEU')ITPAR=3
      IF(CHPAR.EQ.'CUTHAD')ITPAR=4
      IF(CHPAR.EQ.'CUTMUO')ITPAR=5
      IF(CHPAR.EQ.'BCUTE' )ITPAR=6
      IF(CHPAR.EQ.'BCUTM' )ITPAR=7
      IF(CHPAR.EQ.'DCUTE' )ITPAR=8
      IF(CHPAR.EQ.'DCUTM' )ITPAR=9
      IF(CHPAR.EQ.'PPCUTM')ITPAR=10
      IF(CHPAR.EQ.'PAIR'  )ITPAR=11
      IF(CHPAR.EQ.'COMP'  )ITPAR=12
      IF(CHPAR.EQ.'PHOT'  )ITPAR=13
      IF(CHPAR.EQ.'PFIS'  )ITPAR=14
      IF(CHPAR.EQ.'DRAY'  )ITPAR=15
      IF(CHPAR.EQ.'ANNI'  )ITPAR=16
      IF(CHPAR.EQ.'BREM'  )ITPAR=17
      IF(CHPAR.EQ.'HADR'  )ITPAR=18
      IF(CHPAR.EQ.'MUNU'  )ITPAR=19
      IF(CHPAR.EQ.'DCAY'  )ITPAR=20
      IF(CHPAR.EQ.'LOSS'  )ITPAR=21
      IF(CHPAR.EQ.'MULS'  )ITPAR=22
      IF(CHPAR.EQ.'RAYL'  )ITPAR=23
      IF(CHPAR.EQ.'GHCOR1')ITPAR=26
      IF(CHPAR.EQ.'GHCOR2')ITPAR=27
      IF(CHPAR.EQ.'GHCOR3')ITPAR=28
      IF(CHPAR.EQ.'GHCOR4')ITPAR=29
      IF(CHPAR.EQ.'BIRK1' )ITPAR=27
      IF(CHPAR.EQ.'BIRK2' )ITPAR=28
      IF(CHPAR.EQ.'BIRK3' )ITPAR=29
      IF(CHPAR.EQ.'LABS' ) ITPAR=31
      IF(CHPAR.EQ.'SYNC' ) ITPAR=32
      IF(CHPAR.EQ.'STRA' ) ITPAR=33
      IF(ITPAR.NE.0)THEN
         Q(JTMN+ITPAR)=PARVAL
         IF(ITPAR.EQ.21)THEN
            KLOSS=PARVAL+0.001
            IF(KLOSS.EQ.3.OR.KLOSS.EQ.1)Q(JTMN+15)=1.
         ENDIF
      ELSE
         WRITE(CHMAIL,1000)ITMED,CHPAR
         CALL GMAIL(0,0)
      ENDIF
      GO TO 99
C
  90  WRITE(CHMAIL,2000)ITMED
      CALL GMAIL(0,0)
C
 1000 FORMAT(' ***** GSTPAR error for tracking medium ',
     +  I3,' Tracking parameter ',A,'  not defined ***** ')
 2000 FORMAT(' ***** GSTPAR error. Tracking medium NR ',
     +  I3,' not defined ***** ')
  99  END
 
