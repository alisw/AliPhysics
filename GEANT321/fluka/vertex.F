*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLDNAME)
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*=== vertex ===========================================================*
*
      SUBROUTINE VERTEX(IT,LT,LL,KFA,E0,IV,RE,KFR1,KFR2,AME,IOPT)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
      REAL RNDM(2)
C     CHOICE OF THE VERTEX
      DIMENSION RE(*),KFR1(*),KFR2(*),IV(*)
C     WAEHLT PASSENDEN VERTEX IV AUS
       IF(LT.EQ.1) WRITE(LUNOUT,288)IT,LT,LL,KFA,E0,AME,IOPT
  288 FORMAT (4I5,2E12.4,I5,' VERTEX IT,LT,LL,KFA,E0,AME')
      PDB=0.8D0
      IVA=1
      I=IT
      J=IT-1
      IVX=IVA
      IF(IOPT.EQ.2) IVX=2
      IF(IOPT.EQ.4.AND.KFA.GT.6.AND.LL.EQ.0) IVX=2
      IF(IOPT.EQ.4.AND.KFA.LE.6.AND.LL.EQ.1) IVX=2
      IF(IOPT.EQ.5) IVX=2
      IF(I.GT.1)IVX=IV(J)
      IF(LL.EQ.1.AND.I.GT.1)IVX=IV(J)-5
      CALL GRNDM(RNDM,2)
      Z=RNDM(1)
      Y=RNDM(2)
      GO TO (1,2,1,2,2),IVX
    1 IVX=1
      IF(Z.GE.AME) IVX=2
      GO TO 20
    2 IVX=3
      IF(Z.GT.PDB) IVX=4
      IF(Z.GT.PDB.AND.Y.LE.0.5D0)  IVX=5
      GO TO 20
   20 CONTINUE
      IF(LL.EQ.1)IVX=IVX+5
      IV(I)=IVX
      IF(LT.EQ.0)GO TO 60
C     AUSDRUCKEN
      IF(I.EQ.1)GO TO 30
      WRITE(LUNOUT,50)IV(J),IV(I),LL
      GO TO 40
   30 WRITE(LUNOUT,50)IVA,IV(I),LL
   40 CONTINUE
   50 FORMAT(1H0,15HIV(J),IV(I),LL=,3I5)
   60 CONTINUE
      RETURN
      END
#endif
