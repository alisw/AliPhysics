*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:17  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GSTMED(KTMED,NATMED,NMAT,ISVOL,IFIELD,FIELDM,TMAXFD,
     +        STEMAX,DEEMAX,EPSIL,STMIN,UBUF,NWBUF)
*
***********************************************************************
*                                                                     *
*                                                                     *
*       Store tracking media parameters                               *
*                                                                     *
*    Stores the  parameters of the tracking  medium ITMED in  the data*
*  structure JTMED.                                                   *
*  ITMED     tracking medium number 0<ITMED<100                       *
*  NATMED    tracking medium name (up to 20 characters ended by $)    *
*  NMAT      material number corresponding to ITMED                   *
*  ISVOL     =0 if not a sensitive volume                             *
*  IFIELD    = 0  if no magnetic field                                *
*            = -1  reserved for user decision in GUSWIM               *
*            = 1  tracking performed with GRKUTA                      *
*            = 2  tracking performed with GHELIX                      *
*            = 3  tracking performed with GHELX3                      *
*  FIELDM    maximum field value (in Kilogauss)                       *
*  TMAXFD    maximum  angle due  to field  permitted in  one step  (in*
*            degrees)                                                 *
*  STEMAX    Maximum step allowed                                     *
*            (in cm)                                                  *
*  DEEMAX    maximum fractional energy loss in one step               *
*  EPSIL     tracking precision (in cm)                               *
*  STMIN     minimum step  due to energy  loss or  multiple scattering*
*            (in cm)                                                  *
*  UBUF      array of NWBUF additional parameters                     *
*  NWBUF                                                              *
*                                                                     *
*                                                                     *
*          The Tracking Medium data structure JTMED                   *
*          ----------------------------------------                   *
*                                                                     *
*                                         | JTMED                     *
*    NTMED           ITMED                v                           *
*     ..........................................................      *
*     |               | |                | | Standard Trac.media      *
*     ..........................................................      *
*                      | JT                                           *
*                      v                                              *
*                    ..........................                       *
*                    | 1 |                    |                       *
*                    .....                    |                       *
*                    | 2 |  Tracking medium   |                       *
*                    |...|                    |                       *
*                    | 3 |   Name             |                       *
*                    |...|                    |                       *
*                    | 4 |                    |                       *
*                    |...|                    |                       *
*                    | 5 |                    |                       *
*                    ..........................                       *
*                    | 6 |   NMAT             |                       *
*                    |...|....................|                       *
*                    | 7 |   ISVOL            |                       *
*                    |...|....................|                       *
*                    | 8 |   IFIELD           |                       *
*                    |...|....................|                       *
*                    | 9 |   FIELDM           |                       *
*                    |...|....................|                       *
*                    | 10|   TMAXFD           |                       *
*                    |...|....................|                       *
*                    | 11|   STEMAX           |                       *
*                    |...|....................|                       *
*                    | 12|   DEEMAX           |                       *
*                    |...|....................|                       *
*                    | 13|   EPSIL            |                       *
*                    |...|....................|                       *
*                    | 14|   STMIN            |                       *
*                    |...|....................|                       *
*                    | 15|   User words ....  |                       *
*                    ..........................                       *
*   JT = LQ(JTMED-ITMED) pointer to tracking medium ITMED             *
*                                                                     *
*    ==>Called by : <USER>, UGEOM    ,<GXINT> GINC3                   *
*       Author    R.Brun  *********                                   *
*                                                                     *
***********************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gctrak.inc"
      DIMENSION MECA(5,13)
      EQUIVALENCE (MECA(1,1),IPAIR)
      DIMENSION UBUF(1),CUTVEC(10)
      EQUIVALENCE (CUTVEC,CUTGAM)
      CHARACTER*(*) NATMED
      CHARACTER*20 NAME
C.
C.    ------------------------------------------------------------------
C.
      ITMED=ABS(KTMED)
      IF(JTMED.LE.0)THEN
         CALL MZBOOK(IXCONS,JTMED,JTMED,1,'TMED',NTMED,NTMED,40,3,0)
         CALL UCOPY(CUTVEC,Q(JTMED+1),10)
         IQ(JTMED-5)=0
         DO 10 I=1,13
            Q(JTMED+10+I)=MECA(1,I)
   10    CONTINUE
         Q(JTMED+10+21)=ILABS
         Q(JTMED+10+22)=ISYNC
         Q(JTMED+10+23)=ISTRA
      ENDIF
      IF(ITMED.GT.NTMED)THEN
         CALL MZPUSH(IXCONS,JTMED,ITMED-NTMED,0,'I')
         NTMED=ITMED
         JTM1=0
      ELSE
         JTM1=LQ(JTMED-ITMED)
         IF(JTM1.GT.0) THEN
            WRITE(CHMAIL,10100)
            CALL GMAIL(1,0)
            CALL GPTMED(ITMED)
            CALL MZDROP(IXCONS,LQ(JTMED-ITMED),' ')
         ENDIF
      ENDIF
      NW=NWBUF+14
      CALL MZBOOK(IXCONS,JTM,JTMED,-ITMED,'TMED',10,10,NW,IOTMED,0)
C
      NAME=NATMED
      NCH=LNBLNK(NAME)
      IF(NCH.GT.0)THEN
         IF(NAME(NCH:NCH).EQ.'$')NAME(NCH:NCH)=' '
      ENDIF
      CALL UCTOH(NAME,IQ(JTM+1),4,20)
C
      EPS=EPSIL
      IF(EPSIL.LE.0.0) THEN
         WRITE(CHMAIL,10000) ITMED, EPSIL
         CALL GMAIL(0,0)
         EPS=1.E-4
      END IF
      IF(IFIELD.NE.0.AND.FIELDM.EQ.0.0) THEN
         WRITE(CHMAIL,10200) ITMED, IFIELD
         CALL GMAIL(0,0)
      END IF
      IF(IGAUTO.NE.0.AND.ITMED.GT.0)THEN
         DE=-1.
         ST=-1.
         SM=-1.
      ELSE
         DE=DEEMAX
         ST=STMIN
         SM=STEMAX
      ENDIF
      Q(JTM + 6) = NMAT
      Q(JTM + 7) = ISVOL
      Q(JTM + 8) = IFIELD
      Q(JTM + 9) = FIELDM
      Q(JTM + 10) = TMAXFD
      Q(JTM + 11) = SM
      Q(JTM + 12) = DE
      Q(JTM + 13) = EPS
      Q(JTM + 14) = ST
      IF(NWBUF.GT.0)CALL UCOPY(UBUF,Q(JTM+15),NWBUF)
C
      IF(JTM1.GT.0) THEN
         CALL GPTMED(-ITMED)
      ENDIF
C
10000 FORMAT('0*** GSTMED *** Warning, medium = ',I5,
     +       ', value of EPSIL=',E10.3,' reset to 1 micron')
10100 FORMAT(' *** GSTMED *** Warning, tracking medium redefinition:')
10200 FORMAT('0*** GSTMED *** Warning, medium = ',I5,
     +       ', IFIELD = ',I3,' and FIELDM = 0.0 is illegal')
  999 END
