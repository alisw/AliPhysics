*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 03/07/94  19.26.33  by  S.Giani
*-- Author :
      SUBROUTINE GDRAYT(NAME,UTHET,UPHI,UPSI,UU0,UV0,SU,SV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       RAY-TRACING                                              *
C.    *       This routine draws the objects, in the Jvolum data       *
C.    *       data structure, that are visible. The objects are drawn  *
C.    *       at the screen point (UU0,UV0), with the                  *
C.    *       screen factors SU and SV acting on the U and V           *
C.    *       dimensions respectively;                                 *
C.    *       the object is rotated by an angle UTHET along Y-axis     *
C.    *       and UPHI along Z-axis and the resulting 2-D picture      *
C.    *       is also rotated by an angle UPSI along the line of       *
C.    *       projection (i.e. the normal to the 2-D view plane).      *
C.    *                                                                *
C.    *       Author:  S.Giani.                                        *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gconst.inc"
#include "geant321/gcmutr.inc"
********************************************************************************
#include "geant321/gcflag.inc"
#include "geant321/gcstak.inc"
#include "geant321/gcrayt.inc"
#include "geant321/gcpixe.inc"
      COMMON/INIFIR/NFIRST
********************************************************************************
*
      CHARACTER*4 NAME
 
      DIMENSION V(3),T(4,3)
      SAVE IFIRST
      DATA IFIRST/0/
*
      CALL UCTOH('PERS',IPERS,4,4)
*
         IF(NFIRST.EQ.0) THEN
            CALL GDCOTA
            NFIRST = 1
         ENDIF
C
      IF (IFIRST.NE.0) GO TO 40
C
      IFIRST=1
      DPHI=PI/20.
      PHI=0.
C
      DO 30 I=1,40
         GSIN(I)=SIN(PHI)
         GCOS(I)=COS(PHI)
         PHI=PHI+DPHI
   30 CONTINUE
C
      GSIN(41)=GSIN(1)
      GCOS(41)=GCOS(1)
C
   40 CONTINUE
C
C             Theta, phi and psi angles are normalized in [0-360[ range
C
 
      GTHETA=MOD(ABS(UTHET),360.)
      GPHI=MOD(ABS(UPHI),360.)
      GPSI=MOD(ABS(UPSI),360.)
      IMOD=0
*
*              Set Transformation Matrix T as for CG Package
*
         V(1)=GTHETA
         V(2)=GPHI
         V(3)=GPSI
         CALL GDCGVW(V,T)
********************************************************************************
            CALL ISFACI(1)
            CALL ISFAIS(1)
            if(iswit(9).ne.67890)CALL IGBOX(0.,20.,20.,0.)
                           DO 101 IHH=1,3
                              DO 102 JHH=1,4
                                 AROTS(IHH,JHH)=T(JHH,IHH)
  102                         CONTINUE
  101                      CONTINUE
                           AROTS(4,1)=0.
                           AROTS(4,2)=0.
                           AROTS(4,3)=0.
                           AROTS(4,4)=1.
                           DO 105 IHH=1,4
                              DO 106 JHH=1,4
                                 ZROTS(IHH,JHH)=AROTS(IHH,JHH)
  106                         CONTINUE
  105                      CONTINUE
                           CALL RINV(4,ZROTS,4,RRR,IFAIL)
         IFLAPE=0
         IF(IPERS.EQ.IPRJ)THEN
           IFLAPE=1
           IF(DPERS.GT.5000.)DPERS=5000.
           FZUV=8.4+DPERS
                                 FPINTX=ZROTS(1,4)+ZROTS(1,1)*
     +                           10.+ZROTS(1,2)*10.+ZROTS(1,3)*
     +                           FZUV
                                 FPINTY=ZROTS(2,4)+ZROTS(2,1)*
     +                           10.+ZROTS(2,2)*10.+ZROTS(2,3)*
     +                           FZUV
                                 FPINTZ=ZROTS(3,4)+ZROTS(3,1)*
     +                           10.+ZROTS(3,2)*10.+ZROTS(3,3)*
     +                           FZUV
         ENDIF
         ISFILL=0
         CALL GFIND('*','FILL',ISFILL)
         IF(ISFILL.EQ.0)ISFILL=10
*         CALL ISMKSC(0.)
         CALL ISMK(1)
         IMAP=0
         IF(NMAP.GT.4)NMAP=4
         IF(NMAP.GT.0)THEN
          NNIM=2
         ELSE
          NNIM=1
         ENDIF
        DO 211 IM=1,NNIM
         IF(NMAP.GT.0)IMAP=IMAP+1
         IF(IMAP.EQ.1)THEN
           IF(ISFILL.LT.NMAP)ISFILL=NMAP
           ISSAVE=ISFILL
           IOSAVE=IOMBRA
           ISFILL=NMAP
           IOMBRA=0
         ELSEIF(IMAP.EQ.2)THEN
           ZNMAP1=1./NMAP
           ISFILL=ISSAVE
           IOMBRA=IOSAVE
         ENDIF
         UUU=1./ISFILL
         VVV=1./ISFILL
         IXXX=1
         IYYY=1
         IF(IPERS.EQ.IPRJ)THEN
           ZUV=DPERS
         ELSE
           ZUV=5000.
         ENDIF
         LIMPRE=20*ISFILL
*
                                 XPINTS=ZROTS(1,4)+ZROTS(1,1)*
     +                           UUU+ZROTS(1,2)*VVV+ZROTS(1,3)*
     +                           ZUV
                                 YPINTS=ZROTS(2,4)+ZROTS(2,1)*
     +                           UUU+ZROTS(2,2)*VVV+ZROTS(2,3)*
     +                           ZUV
                                 ZPINTS=ZROTS(3,4)+ZROTS(3,1)*
     +                           UUU+ZROTS(3,2)*VVV+ZROTS(3,3)*
     +                           ZUV
*            print *,xpints,ypints,zpints,'vertex from gdraw'
            ISSEEN=0
            CALL GTRIGC
*            IEVENT=IEVENT+1
            CALL GTRIGI
            CALL GTRIG
 211    CONTINUE
*
********************************************************************************
*
 999      END
