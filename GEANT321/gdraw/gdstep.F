*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GDSTEP
*
************************************************************************
*                                                                      *
*       RAY-TRACING                                                    *
*       Light processing at the end of each step                       *
*                                                                      *
*       Author: S.Giani.                                               *
*                                                                      *
************************************************************************
*
#include "geant321/gctmed.inc"
#include "geant321/gcking.inc"
#include "geant321/gcflag.inc"
#include "geant321/gctrak.inc"
********************************************************************************
#include "geant321/gcbank.inc"
#include "geant321/gcvdma.inc"
#include "geant321/gcpixe.inc"
#include "geant321/gcrayt.inc"
 
      DIMENSION PROTU(1),PROTV(1)
      DIMENSION OUT(10)
      DIMENSION CKDIFF(8),CKSPRE(8),NEXPO(8)
      DATA CKDIFF/0.40,0.80,0.50,0.50,0.25,0.25,0.00,0.00/
      DATA CKSPRE/0.00,0.00,0.50,0.50,0.75,0.75,1.00,1.00/
      DATA NEXPO/0,0,2,7,2,7,1,10/
*
*  LSTY keys:
*   0=darkmatt; 1=brightmatt; 2=plastic; 3=ceramic; 4=roughmetal;
*   5=shinymetal; 6=glass; 7=mirror.
*
      SAVE MYMED,APROSH,MYOMED,SPONOR,XYZMOD
********************************************************************************
      BOFLAG=0.
      IF(ISSEEN.EQ.1.OR.ISSEEN.EQ.-2)THEN
      IF(INWVOL.EQ.1.OR.(INWVOL.EQ.3.AND.IOMBRA.EQ.9))THEN
        BOFLAG=1.
                   PROTU(1)=UUU
                   PROTV(1)=VVV
       IF(IOMBRA.NE.9)THEN
       IF(NOFLAG.EQ.0)THEN
         out(1)=veccos(1)*vect(4)+veccos(2)*vect(5)+veccos(3)*vect(6)
         out(2)=veccos(1)
         out(3)=veccos(2)
         out(4)=veccos(3)
         OBSNOR=ABS(OUT(1))
         APROSC=ABS(OUT(1))*CKDIFF(ISLSTY+1)
         if(aprosc.lt.0.)then
            print *,'warning: aprosc.lt.0'
         endif
         IF(ISLSTY.NE.0.AND.ISLSTY.NE.1)THEN
           BPROSC=(2.*(OUT(1)**2))-1.
           IF(BPROSC.GT.0.)THEN
             BPROSC=CKSPRE(ISLSTY+1)*(BPROSC**NEXPO(ISLSTY+1))
           ELSE
             BPROSC=0.
           ENDIF
           APROSC=(APROSC+BPROSC)
         ENDIF
       ELSE
         OBSNOR=APFLAG
         APROSC=APFLAG*CKDIFF(ISLSTY+1)
         if(aprosc.lt.0.)then
           print *,'warning: aprosc.lt.0',aprosc
         endif
         IF(ISLSTY.NE.0.AND.ISLSTY.NE.1)THEN
           BPROSC=(2.*(APFLAG**2))-1.
           IF(BPROSC.GT.0.)THEN
             BPROSC=CKSPRE(ISLSTY+1)*(BPROSC**NEXPO(ISLSTY+1))
           ELSE
             BPROSC=0.
           ENDIF
           APROSC=(APROSC+BPROSC)
         ENDIF
       ENDIF
       ENDIF
       IF(IOMBRA.EQ.10)THEN
         XYZMOD=1./SQRT((XLPOS-VECT(1))**2+(YLPOS-VECT(2))**2+
     +                  (ZLPOS-VECT(3))**2)
         IF(NOFLAG.EQ.0)THEN
*         call ggperp(vect,veccos,ierr)
*         if(ierr.eq.1)then
*           veccos(1)=1.
*           veccos(2)=0.
*           veccos(3)=0.
*         endif
         out(1)=veccos(1)*vect(4)+veccos(2)*vect(5)+veccos(3)*vect(6)
         out(2)=veccos(1)
         out(3)=veccos(2)
         out(4)=veccos(3)
          XLDIR=(XLPOS-VECT(1))*XYZMOD
          YLDIR=(YLPOS-VECT(2))*XYZMOD
          ZLDIR=(ZLPOS-VECT(3))*XYZMOD
          APROSH=((OUT(2)*XLDIR)+(OUT(3)*YLDIR)+(OUT(4)*ZLDIR))*
     +    CKDIFF(ISLSTY+1)
          SPONOR=APROSH/CKDIFF(ISLSTY+1)
          IF(ISLSTY.NE.0.AND.ISLSTY.NE.1)THEN
          IF(SPONOR.GT.0.)THEN
            B1=OUT(4)*YLDIR-OUT(3)*ZLDIR
            B2=OUT(2)*ZLDIR-OUT(4)*XLDIR
            B3=OUT(3)*XLDIR-OUT(2)*YLDIR
            VRX=OUT(4)*B2-OUT(3)*B3+OUT(2)*SPONOR
            VRY=OUT(2)*B3-OUT(4)*B1+OUT(3)*SPONOR
            VRZ=OUT(3)*B1-OUT(2)*B2+OUT(4)*SPONOR
            COSB1=VRX*XCOSXS+VRY*YCOSYS+VRZ*ZCOSZS
            IF(COSB1.GT.0.)THEN
              BPROSH=CKSPRE(ISLSTY+1)*(COSB1**NEXPO(ISLSTY+1))
*              print *,COSB1**ISLSTY
            ELSE
              BPROSH=0.
            ENDIF
            APROSH=(APROSH+BPROSH)
          ENDIF
          ENDIF
         ELSE
          SPONOR=APFLAH
          APROSH=APFLAH*CKDIFF(ISLSTY+1)
          IF(ISLSTY.NE.0.AND.ISLSTY.NE.1)THEN
          IF(SPONOR.GT.0.)THEN
            IF(IXYFLA.EQ.1)THEN
              OUT(2)=CCXX(2)
              OUT(3)=CCXX(3)
              OUT(4)=CCXX(4)
            ELSEIF(IXYFLA.EQ.2)THEN
              OUT(2)=CCYY(2)
              OUT(3)=CCYY(3)
              OUT(4)=CCYY(4)
            ELSEIF(IXYFLA.EQ.3)THEN
              OUT(2)=CCZZ(2)
              OUT(3)=CCZZ(3)
              OUT(4)=CCZZ(4)
            ELSE
              print *,'problem for ixyfla'
            ENDIF
            B1=OUT(4)*YLDIR-OUT(3)*ZLDIR
            B2=OUT(2)*ZLDIR-OUT(4)*XLDIR
            B3=OUT(3)*XLDIR-OUT(2)*YLDIR
            VRX=OUT(4)*B2-OUT(3)*B3+OUT(2)*SPONOR
            VRY=OUT(2)*B3-OUT(4)*B1+OUT(3)*SPONOR
            VRZ=OUT(3)*B1-OUT(2)*B2+OUT(4)*SPONOR
            COSB1=VRX*XCOSXS+VRY*YCOSYS+VRZ*ZCOSZS
            IF(COSB1.GT.0.)THEN
              BPROSH=CKSPRE(ISLSTY+1)*(COSB1**NEXPO(ISLSTY+1))
*              print *,COSB1**ISLSTY
            ELSE
              BPROSH=0.
            ENDIF
            APROSH=(APROSH+BPROSH)
          ENDIF
          ENDIF
         ENDIF
       ENDIF
       SAVMED=NUMED
       IF(IOMBRA.EQ.10)MYOMED=ISCOLO
       IF(IOMBRA.NE.9)THEN
           CALL GDLUMI(ISCOLO,APROSC)
       ENDIF
       IF(IOMBRA.NE.9.AND.IOMBRA.NE.10)THEN
        IF(IMAP.EQ.1)THEN
         JON=1
        ELSE
**         IF(ISWIT(9).NE.54321)THEN
          CALL IPM(1,PROTU,PROTV)
**         ELSE
**          CALL IXPIXEL(IXXX,IYYY,ICOLOR)
**         ENDIF
        ENDIF
       ELSEIF(IOMBRA.EQ.10)THEN
         MYMED=ISCOLO
       ELSEIF(IOMBRA.EQ.9)THEN
         IF((INWVOL.EQ.3.AND.SPONOR.GT.0.).OR.(INWVOL.EQ.1.AND.
     +       SPONOR.GT.0..AND.SLENG.GE.(1./XYZMOD)))THEN
           CALL GDLUMI(MYOMED,APROSH)
           SINTEN=INTEN*.1
           FINMED=((1.-SINTEN)*MYMED)+(SINTEN*MYOMED)
           MEDFIN=FINMED
           IF((FINMED-MEDFIN).GT..5)MEDFIN=MEDFIN+1
           CALL ISPMCI(MEDFIN)
           CALL IPM(1,PROTU,PROTV)
           IF(INWVOL.EQ.3.AND.SPONOR.GT.0.)IOMBRA=10
         ELSE
          SINTEN=INTEN*.1
           IF(MYOMED.EQ.2)THEN
            MYOMED=17
           ELSEIF(MYOMED.EQ.3)THEN
            MYOMED=67
           ELSEIF(MYOMED.EQ.4)THEN
            MYOMED=117
           ELSEIF(MYOMED.EQ.5)THEN
            MYOMED=42
           ELSEIF(MYOMED.EQ.6)THEN
            MYOMED=142
           ELSEIF(MYOMED.EQ.7)THEN
            MYOMED=92
           ENDIF
           FINMED=((1.-SINTEN)*MYMED)+(SINTEN*MYOMED)
           MEDFIN=FINMED
           IF((FINMED-MEDFIN).GT..5)MEDFIN=MEDFIN+1
           CALL ISPMCI(MEDFIN)
           CALL IPM(1,PROTU,PROTV)
           IF(INWVOL.EQ.3.AND.SPONOR.LE.0.)IOMBRA=10
         ENDIF
       ENDIF
       NUMED=SAVMED
      ENDIF
      ENDIF
********************************************************************************
*
      END
