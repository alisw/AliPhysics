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
      SUBROUTINE MICSET(MATNO,NKEY,CNAME)
C***********************************************************************
C set a option in MICAP
C
C INPUT:   MATNO  - GEANT material number
C          NKEY   - 0 -> use single isotopes instead of the
C                                 natural composition in material MATNO
C                   1 -> use natural composition
C                   10-> define additional x-section file (MATNO ignored)
C          CNAME  - in case NKEY=10 the file name
C
C************************************************************************
C
#include "geant321/mmicap.inc"
C
      INTEGER MATNO,NKEY
      CHARACTER*(*) CNAME
C
      LOGICAL FIRST,FIRST10
      DATA FIRST/.TRUE./
      DATA FIRST10/.TRUE./
C
      IF(FIRST.AND.NKEY.LT.10) THEN
         FIRST = .FALSE.
         NWW = 100
         CALL CHKZEB(NWW,IXCONS)
         CALL MZLINK(IXCONS,'MICPAR',LMIST,LMIST,LMIST)
         CALL MZBOOK(IXCONS,LMIST,0,2,'MIST',0,0,NWW,0,0)
      ELSE IF(FIRST10) THEN
         FIRST10 = .FALSE.
         NFIL = 101
         CALL CHKZEB(NFIL,IXCONS)
         CALL MZLINK(IXCONS,'MICFIL',LMIFIL,LMIFIL,LMIFIL)
         CALL MZBOOK(IXCONS,LMIFIL,0,2,'MIFL',0,0,NFIL,0,0)
      ELSE IF(NKEY.EQ.10) THEN
C increase the bank for the x-section file name
        NFIL = 101 + IQ(LMIFIL-1)
        CALL CHKZEB(NFIL,IXCONS)
        CALL MZPUSH(IXCONS,LMIFIL,0,101,'I')
      ENDIF
      IF(NKEY.LT.10) THEN
C store MICAP option for material MATNO in bank 'MIST'
   10   CONTINUE
        DO 20 I=1,IQ(LMIST-1),2
           IF(IQ(LMIST+I).EQ.MATNO) THEN
              IQ(LMIST+I+1) = NKEY
              GOTO 30
           ENDIF
           IF(IQ(LMIST+I).EQ.0) THEN
              IQ(LMIST+I) = MATNO
              IQ(LMIST+I+1) = NKEY
              GOTO 30
           ENDIF
   20   CONTINUE
C
C  Bank got to small, increase the size
        NWW = 100 + IQ(LMIST-1)
        CALL CHKZEB(NWW,IXCONS)
        CALL MZPUSH(IXCONS,LMIST,0,100,'I')
        GOTO 10
      ELSE
C store x-section file name in bank 'MIFL'
C find the last free index in the bank
        IF(LNBLNK(CNAME).GT.0) THEN
          I = LMIFIL+IQ(LMIFIL-1)-100+1
          CALL UCTOH(CNAME,IQ(I),4,LNBLNK(CNAME))
          IQ(I-1) = LNBLNK(CNAME)
        ELSE
          PRINT*,' MICSET : invalid file name '
        ENDIF
      ENDIF
   30 RETURN
      END
