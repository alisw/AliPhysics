*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_OLD)
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GMEDIA (X, NUMED)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Finds in which volume/medium the point X is, and updates the *
C.    *    common /GCVOLU/ and the structure JGPAR accordingly.        *
C.    *                                                                *
C.    *   NUMED returns the tracking medium number, or 0 if point is   *
C.    *         outside the experimental setup.                        *
C.    *                                                                *
C.    *   Called by :  GTREVE, GLTRAC, 'User'                          *
C.    *   Authors   : R.Brun, F.Bruyant, A.McPherson                   *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
C.
      DIMENSION  X(*)
      REAL       XC(3)
      LOGICAL    BTEST
C.
C.    ------------------------------------------------------------------
*
      IF (NLEVEL.EQ.0) CALL GMEDIN
*
* SECTION I: The /GCVOLU/ table contains the initial guess for a path
*            in the geometry tree on which X may be found.  Look along this
*            path until X is found inside.  This is the starting position.
*            If this is an ONLY volume with no daughters, we are done;
*            otherwise reset search record variables, proceed to section II.
*
*            The information contained in INFROM has to be invalidated
*            because it has no meaning for the subsequent tracking. INFR
*            is a local variable used to optimise the search in the
*            geometry tree.
*
      INFROM = 0
*
* *** Check if point is in current volume
*
      INFR   = 0
      JVIN   = 0
C*****  Code Expanded From Routine:  GTRNSF
C
  100 IF (GRMAT(10,NLEVEL) .EQ. 0.) THEN
         XC(1) = X(1) - GTRAN(1,NLEVEL)
         XC(2) = X(2) - GTRAN(2,NLEVEL)
         XC(3) = X(3) - GTRAN(3,NLEVEL)
*
      ELSE
         XL1 = X(1) - GTRAN(1,NLEVEL)
         XL2 = X(2) - GTRAN(2,NLEVEL)
         XL3 = X(3) - GTRAN(3,NLEVEL)
         XC(1) = XL1*GRMAT(1,NLEVEL) + XL2*GRMAT(2,NLEVEL) + XL3*GRMAT(3
     +      ,NLEVEL)
         XC(2) = XL1*GRMAT(4,NLEVEL) + XL2*GRMAT(5,NLEVEL) + XL3*GRMAT(6
     +      ,NLEVEL)
         XC(3) = XL1*GRMAT(7,NLEVEL) + XL2*GRMAT(8,NLEVEL) + XL3*GRMAT(9
     +      ,NLEVEL)
 
      ENDIF
C*****  End of Code Expanded From Routine:  GTRNSF
*
      JVO  = LQ(JVOLUM-LVOLUM(NLEVEL))
      JPAR = LQ(JGPAR-NLEVEL)
      CALL GINME (XC, Q(JVO+2), Q(JPAR+1), IYES)
      IF (IYES.EQ.0) THEN
*
*  **   Point not in current volume, go up the tree
*
         IF (NLEVEL.GT.1) THEN
            NLEVEL = NLEVEL -1
            JVO  = LQ(JVOLUM-LVOLUM(NLEVEL))
            NIN = Q(JVO+3)
            IF(NIN.GT.0) THEN
*
*       Do not set INFR whne going up the tree. GMEDIA can be called
*       by the user and it  should not assume  that the previous
*       position has something to do with the current search. INFR
*       is otherwise useful when searching in a 'MANY' volume
*       configuration. This statement is commented for the above reason.
*
*              INFR  =LINDEX(NLEVEL+1)
            ELSE
               INFR  =0
            ENDIF
            GO TO 100
         ELSE
*
*   *      Point is outside setup
*
            NUMED = 0
            GO TO 999
         ENDIF
      ENDIF
*
*  **   Point is in current volume
*
      IF(INFR  .GT.0) THEN
         JIN=LQ(JVO-INFR  )
         IQ(JIN) = IBSET(IQ(JIN),4)
         JVIN = JIN
      ENDIF
      NLMIN = NLEVEL
      NLMANY = 0
*
* SECTION II: X is found inside current node at NLEVEL in /GCVOLU/.
*             Search all contents recursively for any containing X.
*             Take the first one found, if any, and continue at that
*             level, incrementing NLEVEL and extending /GCVOLU/ tables.
*             This is continued until a level is reached where X is not
*             found in any of the contents, or there are no contents.
* Note: Since Section II is re-entered from Section III, a blocking word
* is used to mark those contents already checked.  Upon exit from Section
* II, these blocking words are cleared at NLEVEL, but may remain set in
* levels between NLEVEL-1 and NLMIN, if any.  They must be cleared at exit.
*
*  **  Check contents, if any
*
  200 JVO = LQ(JVOLUM-LVOLUM(NLEVEL))
      NIN = Q(JVO+3)
*
*   *   Case with no contents
*
      IF (NIN.EQ.0) THEN
         GO TO 300
*
*   *   Case with contents defined by division
*
      ELSEIF (NIN.LT.0) THEN
         CALL GMEDIV (JVO, IN, XC, 1)
         IF (IN.GT.0) THEN
            INFR   = 0
            GO TO 200
         ENDIF
*
*   *  Case with contents positioned
*
      ELSE
         JCONT  = LQ(JVO-NIN-1)+1
         NCONT  = IQ(JCONT)
         ISEARC = Q(JVO+1)
         IF (ISEARC.LT.0) THEN
*
*       Prepare access to contents, when ordered by GSORD
*
            JSB = LQ(LQ(JVO-NIN-1))
            IAX = Q(JSB+1)
            NSB = Q(JSB+2)
            IF (IAX.LE.3) THEN
               IDIV = LOCATF (Q(JSB+3), NSB, XC(IAX))
            ELSE
               CALL GFCOOR (XC, IAX, CX)
               IDIV = LOCATF (Q(JSB+3), NSB, CX)
            ENDIF
            IF (IDIV.LT.0) IDIV = -IDIV
            IF (IDIV.EQ.0) THEN
               IF (IAX.NE.6) GO TO 260
               IDIV = NSB
            ELSEIF (IDIV.EQ.NSB) THEN
               IF (IAX.NE.6) GO TO 260
            ENDIF
            JSC0  = LQ(JVO-NIN-2)
            NCONT = IQ(JSC0+IDIV)
            JCONT = LQ(JSC0-IDIV)
         ELSE
*
*       otherwise, scan contents (possibly a user selection of them)
*
            JNEAR = LQ(JVO-NIN-1)
            IF (ISEARC.GT.0) THEN
#if !defined(CERNLIB_USRJMP)
               CALL GUNEAR (ISEARC, 1, XC, JNEAR)
#endif
#if defined(CERNLIB_USRJMP)
               CALL JUMPT4(JUNEAR,ISEARC, 1, XC, JNEAR)
#endif
            ELSEIF (INFR  .GT.0) THEN
               JNUP = LQ(LQ(JVO-INFR  )-1)
               IF (JNUP.GT.0) THEN
                  JNEAR = JNUP
               ENDIF
            ENDIF
            JCONT = JNEAR +1
            NCONT = IQ(JCONT)
         ENDIF
*
*     For each selected content in turn, check if point is inside
*
         DO 259 ICONT=1,NCONT
            IN = IQ(JCONT+ICONT)
            IF(IN.EQ.0) THEN
*
*     If the value IQ(JCONT+ICONT)=0 then we are back in the mother.
*     So jump to 260, the search is finished. Clean-up should be done
*     only up to ICONT-1, so we set:
*
               NCONT=ICONT-1
               GOTO 260
            ELSE
            JIN = LQ(JVO-IN)
            IF (.NOT.BTEST(IQ(JIN),4)) THEN
               CALL GMEPOS (JVO, IN, XC, 1)
               IF (IN.GT.0) THEN
                  IF (GONLY(NLEVEL).NE.0.) NLMANY = 0
                  INFR   = 0
                  GO TO 200
               ELSE
                  IQ(JIN) = IBSET(IQ(JIN),4)
               ENDIF
            ENDIF
            ENDIF
  259    CONTINUE
*
  260    IF(NCONT.EQ.NIN) THEN
         DO 268 IN=1,NIN
            JIN = LQ(JVO-IN)
            IQ(JIN) = IBCLR(IQ(JIN),4)
  268    CONTINUE
         ELSE
         DO 269 ICONT=1,NCONT
            IN  = IQ(JCONT+ICONT)
            JIN = LQ(JVO-IN)
            IQ(JIN) = IBCLR(IQ(JIN),4)
  269    CONTINUE
         IF(INFR  .GT.0) THEN
            JIN = LQ(JVO-INFR  )
            IQ(JIN) = IBCLR(IQ(JIN),4)
         ENDIF
         ENDIF
*
      ENDIF
*
* SECTION III: X is found at current node (NLEVEL in /GCVOLU) but not in
*              any of its contents, if any.  If this is a MANY volume,
*              save it as a candidate best-choice, and continue the search
*              by backing up the tree one node and proceed to Section II.
*              If this is an ONLY volume, proceed to Section IV.
*
* *** Point is in current volume/medium, and not in any content
*
  300 IF (GONLY(NLEVEL).EQ.0.) THEN
*
*  **   Lowest level is 'NOT ONLY'
*
         IF (NLEVEL.GT.NLMANY) THEN
            CALL GSCVOL
            NLMANY = NLEVEL
         ENDIF
*
*   *   Go up the tree up to a volume with positioned contents
*
  310    INFR   = LINDEX(NLEVEL)
         NLEVEL = NLEVEL -1
         JVO    = LQ(JVOLUM-LVOLUM(NLEVEL))
         NIN    = Q(JVO+3)
         IF (NIN.LT.0) GO TO 310
*
C*****  Code Expanded From Routine:  GTRNSF
C
         IF (GRMAT(10,NLEVEL) .EQ. 0.) THEN
            XC(1) = X(1) - GTRAN(1,NLEVEL)
            XC(2) = X(2) - GTRAN(2,NLEVEL)
            XC(3) = X(3) - GTRAN(3,NLEVEL)
*
         ELSE
            XL1 = X(1) - GTRAN(1,NLEVEL)
            XL2 = X(2) - GTRAN(2,NLEVEL)
            XL3 = X(3) - GTRAN(3,NLEVEL)
            XC(1) = XL1*GRMAT(1,NLEVEL) + XL2*GRMAT(2,NLEVEL) + XL3*
     +      GRMAT(3,NLEVEL)
            XC(2) = XL1*GRMAT(4,NLEVEL) + XL2*GRMAT(5,NLEVEL) + XL3*
     +      GRMAT(6,NLEVEL)
            XC(3) = XL1*GRMAT(7,NLEVEL) + XL2*GRMAT(8,NLEVEL) + XL3*
     +      GRMAT(9,NLEVEL)
 
         ENDIF
C*****  End of Code Expanded From Routine:  GTRNSF
*
         JIN = LQ(JVO-INFR  )
         IQ(JIN) = IBSET(IQ(JIN),4)
         NLMIN = MIN(NLEVEL,NLMIN)
         GO TO 200
      ENDIF
*
* SECTION IV: This is the end of the search.  The current node (NLEVEL
*             in /GCVOLU/) is the lowest ONLY volume in which X is found.
*             If X was also found in any of its contents, they are MANY
*             volumes: the best-choice is the one among them at the greatest
*             level in the tree, and it is stored.  Otherwise the current
*             volume is the solution.  Before exit, all of the blocking
*             words leftover in the tree must be reset to zero.
* Note: A valid structure is assumed, in which no ONLY volumes overlap.
* If this rule is violated, or if a daughter is not entirely contained
* within the mother volume, the results are unpredictable.
*
      DO 419 NL=NLMIN,NLEVEL-1
         JVO = LQ(JVOLUM-LVOLUM(NL))
         NIN = Q(JVO+3)
         DO 418 IN=1,NIN
            JIN = LQ(JVO-IN)
            IQ(JIN) = IBCLR(IQ(JIN),4)
  418    CONTINUE
  419 CONTINUE
*
      IF (NLMANY.GT.0) CALL GFCVOL
      JVO = LQ(JVOLUM-LVOLUM(NLEVEL))
      IF(JVIN.NE.0) IQ(JVIN) = IBCLR(IQ(JVIN),4)
      NUMED = Q(JVO+4)
*                                                             END GMEDIA
  999 IF(JGSTAT.NE.0) CALL GFSTAT(2)
      END
#endif
