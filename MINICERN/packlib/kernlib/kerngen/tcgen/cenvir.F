*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:42  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE CENVIR (CHIN,NINP, CHEX,JLP,JRP, IFLEX)
C
C CERN PROGLIB# M432    CENVIR          .VERSION KERNFOR  4.41  950228
C ORIG. 07/12/93, JZ
C
C-    Copy CHIN(1:NINP) to CHEX(JLP:JRP) substituting environment
C-    variables designated by ${name}; on machines running UNIX
C-    the form $name is also recognised.
C-    For undefined variables the original text is either skipped
C-    or copied if IFLEX is zero or not zero.
C-    Return  NFSLAT  bit 1  undefined env.v. encountered
C-                        2  syntax error
C-                        3  not enough space to translate the complete
C-                                                             line
C-            NGSLAT = 0  if the copy is identical to the original
C-                     1  if changed

      COMMON /SLATE/ NDSLAT,NESLAT,NFSLAT,NGSLAT, DUMMY(36)
      CHARACTER    CHIN*512, CHEX*512


      NNIN = NINP
      NEX  = JLP - 1
      NNEX = JRP
#if !defined(CERNLIB_QS_UNIX)
      LIM  = NNIN - 3
#endif
#if defined(CERNLIB_QS_UNIX)
      LIM  = NNIN - 1
#endif

      MISS = 0
      NFAI = 0
      NOSP = 0
      MODI = 0

C--           look for the start of ${SYMBOL}

      NIN = 0
   21 JS = NIN
   22 IF (JS.LT.LIM) THEN
          JS = ICFIND ('$', CHIN,JS+1,LIM)
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

C--           get env. variable name

      JA  = JS + 1
      IF (CHIN(JA:JA).NE.'{')  THEN
#if !defined(CERNLIB_QS_UNIX)
          GO TO 22
#endif
#if defined(CERNLIB_QS_UNIX)
          JB = ICFIND (' ', CHIN,JA,NNIN)
          IF (JB.EQ.JA)            GO TO 22
          JE = ICNUMU (CHIN,JA,JB-1) - 1
          N  = NDSLAT
          IF (N.EQ.0)              GO TO 22
#endif
        ELSE
          JA = JA + 1
          JE = ICFIND ('}', CHIN,JA,NNIN)
          IF (NGSLAT.EQ.0)         GO TO 27
          N  = JE - JA
          IF (N.EQ.0)              GO TO 27
        ENDIF

C--           copy env. variable value

      CALL GETENVF (CHIN(JA:JA+N-1), CHEX(NEX+1:NNEX))
      IF (NDSLAT.EQ.0)             GO TO 26
      NEX  = NEX + NDSLAT
      MODI = 1
      NIN  = JE
      GO TO 21

C--       environment variable not defined

   26 MISS = 1
      IF (IFLEX.NE.0)              GO TO 22
      MODI = 1
      NIN  = JE
      GO TO 21

C--       syntax problems

   27 NFAI = 1
      GO TO 22

C--       CHEX overflows

   28 NOSP = 1

   29 NESLAT = NEX + 1
      NDSLAT = NESLAT - JLP
      NFSLAT = MISS + 2*NFAI + 4*NOSP
      NGSLAT = MODI
      RETURN
      END
