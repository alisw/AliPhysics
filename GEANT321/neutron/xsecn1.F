*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:00  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 23/02/95  14.46.01  by  S.Giani
*-- Author :
      SUBROUTINE XSECN1(NII,KE,IN,ICOM,IREC,IUNIT,LNUMB,IND,
     +                  BUF,IBUF,INEL)
C       THIS ROUTINE READS THE SECOND RECORD ON INPUT
C       I/O UNIT (MICROS)  (I.E. THE B CONTROL BLOCK)
#include "geant321/minput.inc"
#include "geant321/mmicap.inc"
#include "geant321/mconst.inc"
      DIMENSION BUF(*),IBUF(*),ICOM(*),KE(*),IREC(*),IND(*),IN(*)
      DIMENSION INEL(*),LNUMB(*),IUNIT(*)
      INTEGER NII
C       READ THE B CONTROL BLOCK OFF INPUT I/O UNIT
      LT = LTEMP
      LZ = 1
      IU = 1
   10 CONTINUE
        NU = IQ(LT+NTUNIT)
        NIJ = IQ(LT+NTMPNI)
        LZZ=3*NIJ
        READ(NU,'((8I10))')(IBUF(I),I=LZ,LZZ+LZ)
C       INITIALIZE IND ARRAY AND IREC ARRAY TO ZERO
        DO 20 I=IU,IU+NIJ-1
          IUNIT(I) = NU
   20   CONTINUE
        IU = IU+NIJ
        LZ = LZ + LZZ
        LT = LQ(LT)
      IF(LT.GT.0) GOTO 10
      DO 30 I=1,NII
         INEL(I)=0
         IREC(I)=0
   30 CONTINUE
      DO 40 I=1,NMIX
   40 IND(I)=0
      II=0
      JI=0
   50 II=II+1
CZ      IF(II.GT.NII)GO TO 90
      IF(3*II.GT.LZ)GO TO 90
      NEL=IBUF(3*II-2)
      INEL(II)=NEL
      DO 60 IJ=1,NMIX
C correct element AND the correct unit ?
         IF(NEL.EQ.KE(IJ)) GO TO 70
   60 CONTINUE
      IREC(II)=0
      GO TO 50
   70 I=IN(IJ)
C       ICOM RELATES THE ISOTOPE NUMBER TO THE DICTIONARY NUMBER
      IF(ICOM(I).GT.0) IREC(ICOM(I)) = 0
      ICOM(I)=II
C total length of x-section data in words
      LNUMB(I)  = IBUF(3*II)
      IREC(II)  = IBUF(3*II-1)
C       SET INDICATORS
      DO 80  I=IJ,NMIX
         IF(NEL.NE.KE(I))GO TO 80
         IND(I)=1
         JI=JI+1
   80 CONTINUE
      GO TO 50
   90 RETURN
      END
