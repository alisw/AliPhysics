*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CCOSUB (CHIN,NNIN,CHEX,JLP,JRP, TOKEN,SUBST)
C
C CERN PROGLIB# M432    CCOSUB          .VERSION KERNFOR  4.41  950228
C ORIG. 20/01/94, JZ
C
C-    Copy CHIN(1:NNIN) to CHEX(JL:JR) substituting each occurrence
C-    of TOKEN by SUBST.
C-    Return  NDSLAT   number of characters stored
C-            NESLAT   first character after the last stored
C-            NFSLAT non-zero:  CHEX too small to receive the copy
C-            NGSLAT = 0  if the copy is identical to the original
C-                     1  if changed

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHIN*512, CHEX*512, TOKEN*(*), SUBST*(*)


      NTO = LEN(TOKEN)
      NSU = LEN(SUBST)
      LIM = NNIN+1 - NTO

      NEX  = JLP - 1
      NNEX = JRP

      NOSP = 0
      MODI = 0

C--           look for the start of TOKEN

      NIN = 0
   21 JS = NIN
   22 IF (JS.LT.LIM) THEN
          JS = ICFIND (TOKEN(1:1), CHIN,JS+1,LIM)
          IF (NGSLAT.NE.0)         GO TO 23
        ENDIF

      JS = NNIN + 1
      NGSLAT = 0

C--           copy leading text

   23 N  = JS-1 - NIN
      IF (N.GT.0)  THEN
          NN = MIN (N, NNEX-NEX)
          IF (NN.LE.0)             GO TO 28
          CHEX(NEX+1:NEX+NN) = CHIN(NIN+1:NIN+NN)
          NEX = NEX + NN
          IF (NN.NE.N)             GO TO 28
        ENDIF

      IF (NGSLAT.EQ.0)             GO TO 29

      NIN = JS - 1
      IF (CHIN(NIN+1:NIN+NTO).NE.TOKEN(1:NTO))  GO TO 22

C--           copy substitution

      MODI = 1
      NN   = MIN (NSU, NNEX-NEX)
      IF (NN.GT.0)  CHEX(NEX+1:NEX+NN) = SUBST(1:NN)
      NEX = NEX + NN
      NIN = NIN + NTO
      IF (NN.EQ.NSU)               GO TO 21

C--       CHEX overflows

   28 NOSP = 1

   29 NESLAT = NEX + 1
      NDSLAT = NESLAT - JLP
      NFSLAT = NOSP
      NGSLAT = MODI
      RETURN
      END
