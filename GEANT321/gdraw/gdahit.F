*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:19  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      SUBROUTINE GDAHIT (X,Y,Z,ISYMB,SIZSYM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw an HIT point at (X,Y,Z)                             *
C.    *                                                                *
C.    *       The hit belongs to track -ITR3D (in common GCDRAW)       *
C.    *       if GDAHIT is called from GDHITS or GDCHIT, otherwise     *
C.    *       it belongs to track ITRA (in common GCKINE)              *
C.    *       if GDAHIT is called directly from user hit routines.     *
C.    *                                                                *
C.    *       The following characters may be chosen :                 *
C.    *          ISYMB = -1   hardware character    (fast),            *
C.    *                =  0   crosses X by TVDRAW   (default),         *
C.    *                       (size = SIZSYM in cm, default = 0.1),    *
C.    *                = 840  characters by HPLSOF                     *
C.    *                   .   according to the table p. 25             *
C.    *                   .   of the HPLOT manual                      *
C.    *                  853  (size = SIZSYM in cm)                    *
C.    *             (844,854  white, black asterisk * in addition).    *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GDCHIT, GDHITS                       *
C.    *       Authors : R.Brun, P.Zanarini   **********                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gckine.inc"
      CHARACTER*4 LSYMB(10),CHTEXT
      DIMENSION        HITS(3),UU(2),VV(2),JSYMB(10)
      SAVE JSYMB,LSYMB,IFIRST,IOFF
      DATA             JSYMB/  840 ,  841 ,  842 ,  843 ,  844
     +,                        850 ,  851 ,  852 ,  853 ,  854 /
      DATA             LSYMB/'"<0$','"<1$','"<2$','"<3$','"<4$'
     +,                      '">0$','">1$','">2$','">3$','">4$'/
      DATA IFIRST/0/
C.
C.    ------------------------------------------------------------------
C.
      IF (IFIRST.EQ.0) THEN
         IFIRST=1
         CALL UCTOH('OFF ',IOFF,4,4)
      ENDIF
C
C            Set IOBJ to HIT
C
      IOBJ=3
C
      HITS(1)= X
      HITS(2)= Y
      HITS(3)= Z
C
C            Find mode of drawing
C
      KSYMB      = ISYMB
      IF (ISYMB.NE.-1)   KSYMB=IUCOMP(ISYMB,JSYMB(1),10)
      IF (KSYMB.GT. 0)   CALL UCTOH(LSYMB(KSYMB),MSYMB,4,4)
      SSYMB=SIZSYM
      IF (SSYMB.EQ.0.)   SSYMB=0.1
C
C            If THRZ option is set on (R-Z projection) and
C            the cut is lateral (vertical or horizontal) then :
C            if ITR3D is not set by GDHITS or GDCHIT
C            then use -ITRA from common /GCKINE/
C            else do not modify ITR3D already set
C
      IF (ITHRZ.NE.IOFF) THEN
         IF (ICUT.EQ.1.OR.ICUT.EQ.2) THEN
            IF (ITR3D.EQ.0) ITR3D=-ITRA
         ENDIF
      ELSE
         ITR3D=0
      ENDIF
C
C            Draw the hit
C
      CALL GDFR3D(HITS,1,U,V)
C
CM             Short circuit in case of GMR
C
      IF ( IGMR .NE. 0 ) GO TO 999
CM
C
      IF (KSYMB.LT.0) THEN
C
C            ' ' for hardware points, or 'x' for hardware crosses
C
         CALL IPM(1,U,V)
C
      ELSE IF (KSYMB.GT.0) THEN
C
C            soft-characters 840...854 (see above)
C
         V=V-SSYMB*0.5
         CALL UHTOC(MSYMB,4,CHTEXT,4)
         CALL GDRAWT(U,V,CHTEXT,SSYMB,0.,1,0)
C
      ELSE
C
C            default crosses of size SSYMB cm
C
         DP = 0.5 * SSYMB
         UU(1) = U-DP
         UU(2) = U+DP
         VV(1) = V+DP
         VV(2) = V-DP
         CALL GDRAWV(UU,VV,2)
         UU(1) = U+DP
         UU(2) = U-DP
         VV(1) = V+DP
         VV(2) = V-DP
         CALL GDRAWV(UU,VV,2)
C
      ENDIF
C
  999 RETURN
      END
