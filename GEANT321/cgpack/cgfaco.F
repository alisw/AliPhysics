*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
*
      SUBROUTINE CGFACO(VVX,VVY,IJKLMN,LFULL,APROSC)
*
********************************************************************
*                                                                  *
* Function: Fill fully visible faces with solid colour.            *
*           The colour is determined by the table in GDCOTA.       *
*           The intensity is determined in GDSHAD.                 *
*                                                                  *
* Author: S. Giani                                                 *
*                                                                  *
* I/O parameters:                                                  *
*                                                                  *
* VVX,VVY = set of edge's coordinates                              *
* IJKLMN  = number of edges *2.                                    *
* LFULL   = flag for full visibility                               *
* APROSC  = light intensity                                        *
*                                                                  *
********************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cghpar.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
*
      DIMENSION VVX(500),VVY(500),IZ(500)
      DIMENSION TVVX(500),TVVY(500),IACTU(500)
*
      IACTU(1)=0
      IACT=1
      LINFIL=IBITS(LINATT,13,3)
      CALL ISFAIS(1)
      LINCOL=IBITS(LINATT,16,8)
      MPRECM=9-LINFIL
      IF(MPRECM.EQ.8)THEN
         CALL ISFACI(LINCOL)
         CALL ISPLCI(LINCOL)
      ELSE
         CALL GDSHAD(LINCOL,APROSC)
      ENDIF
*
      DO 10  I=1,IJKLMN
         TVVX(I)=VVX(I)
         TVVY(I)=VVY(I)
   10 CONTINUE
      VVX(1)=TVVX(1)
      VVY(1)=TVVY(1)
      VVX(2)=TVVX(2)
      VVY(2)=TVVY(2)
      IZ(1)=1
      IZ(2)=2
      KZ=3
      IGOR=IJKLMN/2
      DO 90  II=3,IGOR
         DO 40  JJ=3,IJKLMN
            DO 30 I=1,KZ-1
               IF(JJ.EQ.IZ(I))GOTO 40
   30       CONTINUE
            C11=ABS(TVVX(JJ)-VVX(II-1))
            C12=ABS(TVVY(JJ)-VVY(II-1))
            IF(C11.LT..001.AND.C12.LT..001)THEN
               PDJ=JJ*.5
               IPDJ=PDJ
               ARG=PDJ-IPDJ
               IF(ARG.GT..49)THEN
                  VVX(II)=TVVX(JJ+1)
                  VVY(II)=TVVY(JJ+1)
                  IZ(KZ)=JJ
                  IZ(KZ+1)=JJ+1
                  KZ=KZ+2
               ELSE
                  VVX(II)=TVVX(JJ-1)
                  VVY(II)=TVVY(JJ-1)
                  IZ(KZ)=JJ
                  IZ(KZ+1)=JJ-1
                  KZ=KZ+2
               ENDIF
               GOTO 90
            ENDIF
   40    CONTINUE
         IACT=IACT+1
         IACTU(IACT)=II-2
         NUMPON=IACTU(IACT)-IACTU(IACT-1)
         NUPO=IACTU(IACT-1)+1
         IF(NUMPON.GT.2)THEN
*            IF(LEP.LT.0)THEN
*               DO 50  IIJ=NUPO,NUPO+NUMPON-1
*                  IF(VVY(IIJ).LT.1.)VVY(IIJ)=1.
*   50          CONTINUE
*            ENDIF
            IF(IDVIEW.EQ.0)THEN
               CALL IFA(NUMPON,VVX(NUPO),VVY(NUPO))
            ELSE
               LLLINA=LINATT
               CALL MVBITS(LINCOL,0,8,LINATT,16)
               CALL GVIEWF(VVX(NUPO),VVY(NUPO),NUMPON)
               LINATT=LLLINA
            ENDIF
         ENDIF
         DO 70 JJA=3,IJKLMN
            DO 60 I=1,KZ-1
               IF(JJA.EQ.IZ(I))GOTO 70
   60       CONTINUE
            JJAO=JJA
            GOTO 80
   70    CONTINUE
         PRINT *,'Error in CGFACO'
   80    CONTINUE
         PDJJA=JJAO*.5
         IPDJJA=PDJJA
         ARG1=PDJJA-IPDJJA
         IF(ARG1.LT..01)PRINT *,'Error in CGFACO'
         VVX(II-1)=TVVX(JJAO)
         VVY(II-1)=TVVY(JJAO)
         VVX(II)=TVVX(JJAO+1)
         VVY(II)=TVVY(JJAO+1)
         IZ(KZ)=JJAO
         IZ(KZ+1)=JJAO+1
         KZ=KZ+2
   90 CONTINUE
*
      IACT=IACT+1
      IACTU(IACT)=IGOR
      NUMPON=IACTU(IACT)-IACTU(IACT-1)
      NUPO=IACTU(IACT-1)+1
      IF(NUMPON.GT.2)THEN
*         IF(LEP.LT.0)THEN
*            DO 100 IIJ=NUPO,NUPO+NUMPON-1
*               IF(VVY(IIJ).LT.1.)VVY(IIJ)=1.
*  100       CONTINUE
*         ENDIF
         IF(IDVIEW.EQ.0)THEN
            CALL IFA(NUMPON,VVX(NUPO),VVY(NUPO))
         ELSE
            LLLINA=LINATT
            CALL MVBITS(LINCOL,0,8,LINATT,16)
            CALL GVIEWF(VVX(NUPO),VVY(NUPO),NUMPON)
            LINATT=LLLINA
         ENDIF
      ENDIF
*
*      IF(IREP.EQ.0)PRINT *,IREP
      DO 120 I=1,KZ
         IZ(I)=0
         IACTU(I)=0
  120 CONTINUE
*
      END
