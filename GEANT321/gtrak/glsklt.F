*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GLSKLT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GLSKLT                                                *
C.    *                                                                *
C.    *   Prepares the Skeleton banks for parallel tracking            *
C.    *   Also lifts the stack bank JTRACK                             *
C.    *                                                                *
C.    *   Called by : GTREVE                                           *
C.    *   Authors   : S.Banerjee, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
C.
      COMMON /GCSKLT/ LOCAL(2), JSK, JSKL, JVOLX
      CHARACTER*12    CFORM
C.
C.    ------------------------------------------------------------------
*
      IF (JSKLT.EQ.0) THEN
*
*  **    Initialize a temporary link area
*
         CALL MZLINT (IXSTOR, '/GCSKLT/', LOCAL, JSK, JVOLX)
         CALL MZFORM ('SKIN', '2I -F', IOSKIN)
         JVOLX = LQ(JVOLUM)
         NLVT  = IQ(JVOLX-1)
*
*  **    Lift the top level bank
*
         CALL MZBOOK (IXCONS, JSKLT, JSKLT, 1, 'SKLT', NLVT, NLVT, 0,
     +                2, 0)
         CALL MZBOOK (IXCONS, JSKL, JSKLT, -1, 'SKLV', 1, 1, 1, 2, 0)
         CALL MZBOOK (IXCONS, JSK, JSKL, -1, 'SKIN', 1, 0, 2, 2, -1)
         JVO  = LQ(JVOLUM-1)
         LQ(JSK-1) = JVO + 6
         IQ(JSK+1) = Q(JVO+5)
         IQ(JSK+2) = 1
*
*  **    Loop over the remaining levels
*
         DO 15 ILEV = 2, NLVT
            NINSK = IQ(JVOLX+ILEV)
            ND    = NINSK + ILEV - 1
            CALL MZBOOK (IXCONS, JSKL, JSKLT,-ILEV,'SKLV', NINSK, NINSK,
     +                   ND, 2, 0)
            DO 10 IN = 1, NINSK
               CALL MZBOOK (IXCONS, JSK, JSKL, -IN, 'SKIN', 1, 0, 15,
     +                      IOSKIN, 1)
   10       CONTINUE
   15    CONTINUE
*
*  **    Now create the Stack bank JTRACK
*
         NWR    = NWTRAC - NWINT
         WRITE (CFORM, 1001) NWINT, NWR
         CALL MZFORM ('TRAC', CFORM, IOTRAC)
         ND     = NWTRAC * NJTMAX
         CALL MZBOOK (IXCONS, JTRACK, JTRACK, 1, 'TRAC', 0, 0, ND,
     +                IOTRAC, -1)
*
         LOCAL(1) = 0
*
      ELSE
*
*  **    Clear the pointers in the skeleton
*
         DO 25 ILEV = 1, NLEVMX
            JSKL = LQ(JSKLT-ILEV)
            DO 20 I = 1, IQ(JSKL-3)
               IQ(JSKL+I) = 0
   20       CONTINUE
   25    CONTINUE
      ENDIF
*
* *** Fill up the skeleton upto NLEVEL
*
      IF (NLEVEL.GT.1) THEN
         DO 60 ILEV = 2, NLEVEL
            JSKL  = LQ(JSKLT-ILEV)
            NINSK = LINMX(ILEV)
            JOFF  = JSKL + IQ(JSKL-3)
            DO 40 IL = 1, ILEV-1
               IF (IQ(JOFF+IL).EQ.LINDEX(IL)) GO TO 40
               DO 30 I = IL, ILEV-1
                  IQ(JOFF+I) = LINDEX(I)
   30          CONTINUE
               DO 35 I = 1, NINSK
                 JSK  = LQ(JSKL-I)
                 IQ(JSK+1) = 0
   35          CONTINUE
               GO TO 45
   40       CONTINUE
   45       JSK  = LQ(JSKL-LINDEX(ILEV))
            IF (IQ(JSK+1).LE.0) THEN
               LQ(JSK-1) = LQ(JGPAR-ILEV)
               IQ(JSK+1) = IQ(JGPAR+ILEV)
               IQ(JSK+2) = LVOLUM(ILEV)
               DO 50 I = 1, 3
                  Q(JSK+2+I) = GTRAN(I,ILEV)
   50          CONTINUE
               DO 55 I = 1, 10
                  Q(JSK+5+I) = GRMAT(I,ILEV)
   55          CONTINUE
            ENDIF
   60    CONTINUE
      ENDIF
*
* *** Initialize pointers
*
      NJFREE = 1
      NJGARB = 0
      NJINVO = 0
      NLDOWN = 1
*
 1001 FORMAT ('/ ',I3,'I ',I3,'F ')
*                                                             END GLSKLT
      END
