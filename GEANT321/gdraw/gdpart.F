*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDPART(IT,ISEL,SIZE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draws attributes of particle generating track IT         *
C.    *       (for which space points have been stored in banks JXYZ). *
C.    *       IT=0 means all tracks.                                   *
C.    *                                                                *
C.    *       If ISEL=xxxxxx1   then draw track number;                *
C.    *       if ISEL=xxxxx1x   then draw particle name.               *
C.    *                                                                *
C.    *       SIZE is character size used.                             *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
      CHARACTER*4 IPTAB(100),IPTAB1(100),IPTAB2(100)
      CHARACTER*5 CHPART
      CHARACTER*12 CHTEXT
      SAVE IPTAB,IPTAB1,IPTAB2
C
      DATA IPTAB /'<[G$','<E>[','<E>[','<[N$','<[M>','<[M>','<[P>',
     +            '<[P>','<[P>','K?L&','K[^+','K[^-','<N$ ','<P$ ',
     +            '<-P$','K?<S','<[C$','[L$ ','[S^+','[S^0','[S^-',
     +            '[X^0','[X^-','[W^-','<-N$','<->[','<->[','<->[',
     +            '<->[','<->[','<->[','<->[','<[T>','<[T>','D[^+',
     +            'D[^-','D[^0','<->D','F[^+','F[^-','[L]?','W[^+',
     +            'W[^-','Z[^0','D$  ','T$  ','<[A$','G$  ',
     +                52*'"4$ '/
      DATA IPTAB1/'    ','^+$ ','^-$ ','    ','^+$ ','^-$ ','^0$ ',
     +            '^+$ ','^-$ ','[^0$','$   ','$   ','    ','    ',
     +            '    ','>&[^','    ','    ','$   ','$   ','$   ',
     +            '$   ','$   ','$   ','    ','L$  ','S^-$','S^0$',
     +            'S^+$','X^0$','X^+$','W^+$','^+$ ','^-$ ','$   ',
     +            '$   ','$   ','[^0$','$   ','$   ','C&[^','$   ',
     +            '$   ','$   ','    ','    ','    ','    ',
     +                52*'    '/
      DATA IPTAB2/'    ','    ','    ','    ','    ','    ','    ',
     +            '    ','    ','L   ','    ','    ','    ','    ',
     +            '    ','0$  ','    ','    ','    ','    ','    ',
     +            '    ','    ','    ','    ','    ','    ','    ',
     +            '    ','    ','    ','    ','    ','    ','    ',
     +            '    ','    ','    ','    ','    ','+$  ','    ',
     +            '    ','    ','    ','    ','    ','    ',
     +                52*'    '/
C.
C.    ------------------------------------------------------------------
C.
#if defined(CERNLIB_DEBUG)
      IF (IT.LT.0) THEN
         IPART=0
         DO 10 IV=19,1,-2
            DO 10 IU=1,19,2
               IPART=IPART+1
               CHTEXT=IPTAB(IPART)//IPTAB1(IPART)//IPTAB2(IPART)
               U=IU
               V=IV
               CALL GDRAWT(U,V,CHTEXT,SIZE,0.,1,0)
   10    CONTINUE
      ENDIF
#endif
C
      CALL UCTOH('OFF ',IOFF,4,4)
      CALL UCTOH('$   ',IDOL,4,4)
C
C            Set IOBJ to TRACK
C
      IOBJ=2
C
      ISEL1=ISEL/10
      ISEL1=ISEL1-ISEL1/2*2
      ISEL0=ISEL/1
      ISEL0=ISEL0-ISEL0/2*2
C
      IF(JXYZ.LE.0)GO TO 30
      NXYZ=IQ(JXYZ-2)
      N1  =1
      N2  =NXYZ
      IF(IT.NE.0)N1=IT
      IF(IT.NE.0)N2=IT
      IF(N1.LE.0)GO TO 30
      IF(N2.GT.NXYZ)GO TO 30
C
      DO 20 I=N1,N2
C
         JX=LQ(JXYZ-I)
         IF(JX.LE.0)GO TO 20
         IPOINT=JX+3
         IF(IPOINT.GT.JX+IQ(JX-1))GO TO 20
         NPOINT=Q(IPOINT)
         IF(NPOINT.LE.0)GO TO 20
         IPART=Q(IPOINT+1)
C
C            Is THRZ option set on (R-Z projection)
C            and the cut lateral (vertical or horizontal) ?
C
         IF (ITHRZ.NE.IOFF) THEN
            IF (ICUT.EQ.1.OR.ICUT.EQ.2) ITR3D=-I
         ELSE
            ITR3D=0
         ENDIF
C
         CALL GDFR3D(Q(IPOINT+3*NPOINT-1),1,U0,V0)
C
         IF (ISEL0.NE.0) THEN
C
C             Draw track number I
C
            WRITE (CHPART,'(I4,''$'')') I
            U=U0-0.5*SIZE
            V=V0-1.5*SIZE
            CALL GDRAWT(U,V,CHPART,SIZE,0.,1,0)
         ENDIF
C
         IF (ISEL1.NE.0) THEN
C
C             Draw particle name IPTAB...(IPART)
C
            CHTEXT=IPTAB(IPART)//IPTAB1(IPART)//IPTAB2(IPART)
            U=U0-0.5*SIZE
            V=V0+0.5*SIZE
            CALL GDRAWT(U,V,CHTEXT,SIZE,0.,1,0)
         ENDIF
C
   20 CONTINUE
C
   30 CONTINUE
C
C             Reset ITR3D to 'standard projection'
C
      ITR3D=0
      IOBJ=0
C
      END
