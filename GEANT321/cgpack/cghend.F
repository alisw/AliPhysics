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
      SUBROUTINE CGHEND(H,JERR,RSHD)
************************************************************************
*                                                                      *
*     Name: CGHEND                                                     *
*     Author: E. Chernyaev                       Date:    07.08.88     *
*     Revised by: S.Giani                        Date:    25.04.91     *
*                                                                      *
*     Function: Terminate creation of HIDE-structure                   *
*                                                                      *
*     References: CGHSHI, CGHDFA, CGHTRE                               *
*                                                                      *
*     Input : H(*) - HIDE-structure                                    *
*     Output:        H(1) .EQ. 0 - error in input data (see CGHCRE)    *
*                    H(1) .LT. 0 - no space                            *
*                                                                      *
*             IHSERR : Flag for indicating the Memory Size Overflow    *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
*SG
#include "geant321/cghpar.inc"
*
*  I4SIZE is now = 1 to compute the right pointers after that all the
*  integer*2 variables have been eliminated.
*SG
      REAL H(*)
*JS
      COMMON /CGOVER/ IHSERR
      JERR=0
*JS
*-
      IF (H(KHSIZE) .EQ. 0.)    GOTO 999
      NHSIZE = H(KHSIZE)
      NT     = H(KHNT)
      NXYZ   = H(KHNXYZ)
      NFACE  = H(KHNFAC)
      NEDGE  = H(KHNEDG)
*           S P A C E   F O R   N O D E S
      NX     = 3*NXYZ
*           S P A C E   F O R   F A C E S
      NF     = (NFACE+NEDGE*2)/I4SIZE + 1
      NP     = NFACE
      NDF    = 6*NFACE
*           S P A C E   F O R   T R E E
      NTR    = 4*NFACE/I4SIZE
      NST    = NFACE/I4SIZE + 1
      NAL    = NFACE
      NAR    = NFACE
*            C H E C K   S P A C E
      NTOTAL = LHHEAD + NX + NF + NP + NDF + NTR + NST + NAL + NAR
      IF (NTOTAL .LT. NHSIZE)   GOTO 100
**SG
*  If size of Hide Structure is too big, tell GDRAW the new number
*    of words needed to build it.
*
      H(KHSIZE) = NHSIZE - NTOTAL
      JERR=NTOTAL+1
      GOTO 999
**SG
*
**          S E T   P O I N T E R S
*
  100 H(KHSIZE) = NTOTAL
*SG
      RSHD=NTOTAL+1
*SG
      JXYZ1      = LHHEAD + 1
      JFA       = JXYZ1 + NX
      JPFA      = JFA  + NF
      JDFA      = JPFA + NP
      JTRE      = JDFA + NDF
      JSTA      = JTRE + NTR
      JALE      = JSTA + NST
      JARI      = JALE + NAL
      H(KHJXYZ) = JXYZ1
      H(KHJFA)  = JFA
      H(KHJPFA) = JPFA
      H(KHJDFA) = JDFA
      H(KHJTRE) = JTRE
      H(KHJSTA) = JSTA
      H(KHJALE) = JALE
      H(KHJARI) = JARI
*
**          S H I F T   F A C E S
**          C O M P U T E   F A C E S   M I N - M A X
**          B U I L D   T R E E
*
      NH        = NHSIZE*I4SIZE
      CALL CGHSHI(NH,NFACE,NEDGE,H(1),H(JFA),H(JPFA))
      H(KHNFAC) = NFACE
      H(KHNEDG) = NEDGE
      CALL CGHDFA(NFACE,H(JXYZ1),H(JFA),H(JPFA),H(JDFA))
      JERR=IHSERR
      IF(IHSERR.NE.0) RETURN
      CALL CGHTRE(NFACE,H(JDFA),H(JSTA),H(JTRE),H(JALE),H(JARI))
*
  999 RETURN
      END
