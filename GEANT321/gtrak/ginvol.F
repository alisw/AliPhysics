*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if !defined(CERNLIB_OLD)
*CMZ :  3.21/02 03/07/94  17.14.15  by  S.Giani
*-- Author :
      SUBROUTINE GINVOL (X, ISAME)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GINVOL (X, ISAME*)                                    *
C.    *                                                                *
C.    *   Checks if particle at point X has left current volume/medium *
C.    *   If so, returns ISAME = 0 and prepares information useful to  *
C.    *    identify the new volume entered.                            *
C.    *   Otherwise, returns ISAME = 1                                 *
C.    *                                                                *
C.    *   Note : INGOTO is set by GTNEXT, to transmit the information  *
C.    *       on the one volume which has limited the step SNEXT,      *
C.    *       >0 : INth content                                        *
C.    *       =0 : current volume                                      *
C.    *       <0 : -NLONLY, with NLONLY defined as the first 'ONLY'    *
C.    *           level up in the tree for the 'NOT-ONLY' volume       *
C.    *           where the point X is found to be.                    *
C.    *                                                                *
C.    *   Called by : GNEXT, GTELEC, GTHADR, GTMUON, GTNEXT            *
C.    *   Authors   : S.Banerjee, R.Brun, F.Bruyant                    *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#include "geant321/gchvir.inc"
C.
      DIMENSION  X(*)
      REAL       XC(6), XT(3)
      LOGICAL    BTEST
C.
C.    ------------------------------------------------------------------
*
* SECTION I: The /GCVOLU/ table contains the presumed location of X in the
*            geometry tree, at level NLEVEL.  The suggestion is that INGOTO
*            is the index of a content at NLEVEL which may also contain X.
*            If this is so, ISAME=0 and return.  INGOTO is left unchanged.
*            If this is not so, have we left the volume at NLEVEL altogether?
*            If so, ISAME=0 and INGOTO=0, return.  Otherwise, this is the
*            starting position for a search.  Reset search record variables
*            and proceed to section II.
*
* *** Check if point is in current volume
*
      INGT = 0
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
         XC(1) = XL1*GRMAT(1,NLEVEL) + XL2*GRMAT(2,NLEVEL) + XL3*
     +      GRMAT(3,NLEVEL)
         XC(2) = XL1*GRMAT(4,NLEVEL) + XL2*GRMAT(5,NLEVEL) + XL3*
     +      GRMAT(6,NLEVEL)
         XC(3) = XL1*GRMAT(7,NLEVEL) + XL2*GRMAT(8,NLEVEL) + XL3*
     +      GRMAT(9,NLEVEL)
*
      ENDIF
      xc(4)=0.
      xc(5)=0.
      xc(6)=0.
C*****  End of Code Expanded From Routine:  GTRNSF
*
      JVO  = LQ(JVOLUM-LVOLUM(NLEVEL))
*
* Note: At entry the variable INGOTO may contain the index of a volume
* contained within the current one at NLEVEL.  If so, begin by checking
* if X lies inside.  This improves the search speed over that of GMEDIA.
*
      NIN = Q(JVO+3)
      IF ((INGOTO.LE.0).OR.(INGOTO.GT.NIN)) THEN
         INGOTO = 0
      ELSE
*
* ***   Entrance in content INGOTO predicted by GTNEXT
*
         JIN  = LQ(JVO-INGOTO)
         IVOT = Q(JIN+2)
         JVOT = LQ(JVOLUM-IVOT)
         JPAR = LQ(JGPAR-NLEVEL-1)
*
         IROTT = Q(JIN+4)
C*****  Code Expanded From Routine:  GITRAN
C.
C.    ------------------------------------------------------------------
C.
         IF (IROTT .EQ. 0) THEN
            XT(1) = XC(1) - Q(5+JIN)
            XT(2) = XC(2) - Q(6+JIN)
            XT(3) = XC(3) - Q(7+JIN)
*
         ELSE
            XL1 = XC(1) - Q(5+JIN)
            XL2 = XC(2) - Q(6+JIN)
            XL3 = XC(3) - Q(7+JIN)
            JR = LQ(JROTM-IROTT)
            XT(1) = XL1*Q(JR+1) + XL2*Q(JR+2) + XL3*Q(JR+3)
            XT(2) = XL1*Q(JR+4) + XL2*Q(JR+5) + XL3*Q(JR+6)
            XT(3) = XL1*Q(JR+7) + XL2*Q(JR+8) + XL3*Q(JR+9)
*
         ENDIF
C*****  End of Code Expanded From Routine:  GITRAN
*
*   *   Check if point is in content
*
         CALL GINME (XT, Q(JVOT+2), Q(JPAR+1), IYES)
         IF (IYES.NE.0) THEN
*
*          If so, prepare information for volume retrieval, and return
*
            NLEVIN = NLEVEL +1
            LVOLUM(NLEVIN) = IVOT
            NAMES(NLEVIN)  = IQ(JVOLUM+IVOT)
            NUMBER(NLEVIN) = Q(JIN+3)
            LINDEX(NLEVIN) = INGOTO
            LINMX(NLEVIN)  = Q(JVO+3)
            GONLY(NLEVIN)  = Q(JIN+8)
            IF (LQ(LQ(JVOLUM-IVOT)).EQ.0) THEN
               NLDEV(NLEVIN) = NLDEV(NLEVEL)
            ELSE
               NLDEV(NLEVIN) = NLEVIN
            ENDIF
            CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), Q(JIN+5),
     +                   IROTT, GTRAN(1,NLEVIN), GRMAT(1,NLEVIN))
            ISAME = 0
            GO TO 999
         ENDIF
      ENDIF
*
* End of INGOTO processing
*
      JPAR = LQ(JGPAR-NLEVEL)
      CALL GINME (XC, Q(JVO+2), Q(JPAR+1), IYES)
      IF (IYES.EQ.0) THEN
         ISAME  = 0
         INGOTO = 0
         GO TO 999
      ENDIF
*
*  **   Point is in current volume
*
      NLEVIN = NLEVEL
      NLMIN = NLEVEL
      IF ((INFROM.LE.0).OR.(INFROM.GT.NIN)) THEN
         INFROM = 0
      ENDIF
      INFR = INFROM
      NLMANY = 0
      IF (INGOTO.GT.0) THEN
         INGT = INGOTO
         JIN = LQ(JVO-INGOTO)
         IQ(JIN) = IBSET(IQ(JIN),4)
      ENDIF
*
* SECTION II: X is found inside current node at NLEVEL in /GCVOLU/.
*             Search all contents for any containing X.  Take the
*             first one found, incrementing NLEVEL and extending the
*             /GCVOLU/ tables.  Otherwise if the list of contents is
*             exhausted without finding X inside, proceed to Section III.
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
            IF ((GONLY(NLEVEL).EQ.0).AND.
     +          (NLEVEL.LE.NLEVIN)) THEN
                INFR = 0
                INGT = 0
                GO TO 200
             ELSE
                GO TO 450
             ENDIF
         ENDIF
*
*   *  Case with contents positioned
*
      ELSE
       if(nin.gt.1)then
        clmoth=q(jvirt+4*(LVOLUM(NLEVEL)-1)+3)
        chmoth=q(jvirt+4*(LVOLUM(NLEVEL)-1)+4)
        ndivto=q(jvirt+4*(LVOLUM(NLEVEL)-1)+2)
        iaxis =q(jvirt+4*(LVOLUM(NLEVEL)-1)+1)
        if(iaxis.le.3)then
          ivdiv=((xc(iaxis)-clmoth)*ndivto/(chmoth-clmoth))+1
            if(ivdiv.lt.1)then
              ivdiv=1
            elseif(ivdiv.gt.ndivto)then
              ivdiv=ndivto
            endif
        else
          call gfcoor(xc,iaxis,cx)
          if(iaxis.eq.6)then
            if((cx-clmoth).lt.-1.)then
              cx=cx+360.
            elseif((cx-chmoth).gt.1.)then
              cx=cx-360.
            endif
            if(cx.gt.chmoth)then
              cx=chmoth
            elseif(cx.lt.clmoth)then
              cx=clmoth
            endif
          endif
          ivdiv=((cx-clmoth)*ndivto/(chmoth-clmoth))+1
            if(ivdiv.lt.1)then
              ivdiv=1
            elseif(ivdiv.gt.ndivto)then
              ivdiv=ndivto
            endif
        endif
        jvdiv=lq(jvirt-LVOLUM(NLEVEL))
        iofset=iq(jvdiv+ivdiv)
        ncont=iq(jvdiv+iofset+1)
        jcont=jvdiv+iofset+1
        if(ncont.eq.0)goto 260
       else
         JCONT = LQ(JVO-NIN-1)+1
         NCONT = 1
       endif
*
*     For each selected content in turn, check if point is inside
*
         DO 259 ICONT=1,NCONT
           if(nin.eq.1)then
            in=1
           else
            IN = IQ(JCONT+ICONT)
           endif
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
                  IF ((GONLY(NLEVEL).EQ.0).AND.
     +                (NLEVEL.LE.NLEVIN)) THEN
                     INFR = 0
                     INGT = 0
                     GO TO 200
                  ELSE
                     GO TO 450
                  ENDIF
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
           if(nin.eq.1)then
            in=1
           else
            IN  = IQ(JCONT+ICONT)
           endif
            JIN = LQ(JVO-IN)
            IQ(JIN) = IBCLR(IQ(JIN),4)
  269    CONTINUE
         IF(INFR.NE.0) THEN
            JIN = LQ(JVO-INFR)
            IQ(JIN) = IBCLR(IQ(JIN),4)
         ENDIF
         IF(INGT.NE.0) THEN
            JIN = LQ(JVO-INGT)
            IQ(JIN) = IBCLR(IQ(JIN),4)
         ENDIF
         ENDIF
         ingt=0
*
      ENDIF
*
* SECTION III: X is found at current node (NLEVEL in /GCVOLU/) but not in
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
         IF (NLMANY.EQ.0) THEN
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
            XC(1) = XL1*GRMAT(1,NLEVEL) + XL2*GRMAT(2,NLEVEL) +
     +      XL3* GRMAT(3,NLEVEL)
            XC(2) = XL1*GRMAT(4,NLEVEL) + XL2*GRMAT(5,NLEVEL) +
     +      XL3* GRMAT(6,NLEVEL)
            XC(3) = XL1*GRMAT(7,NLEVEL) + XL2*GRMAT(8,NLEVEL) +
     +      XL3* GRMAT(9,NLEVEL)
 
         ENDIF
C*****  End of Code Expanded From Routine:  GTRNSF
*
         JIN = LQ(JVO-INFR)
         IQ(JIN) = IBSET(IQ(JIN),4)
         NLMIN = MIN(NLEVEL,NLMIN)
         GO TO 200
      ENDIF
*
* SECTION IV: This is the end of the search.
*             (1) Entry at 400:  ISAME = 1     The current node (NLEVEL
*             in /GCVOLU/) is an ONLY volume and there were no contents
*             in the tree below it which could claim X.
*             (2) Entry at 450:  ISAME = 0    Section II has just found
*             another volume which has more claim to X than the current
*             one: either another ONLY or a deeper MANY was found.
* Note: A valid structure is assumed, in which no ONLY volumes overlap.
* If this rule is violated, or if a daughter is not entirely contained
* within the mother volume, the results are unpredictable.
*
  400 ISAME = 1
      GOTO 480
 
  450 ISAME = 0
 
  480 DO 489 NL=NLMIN,NLEVEL-1
         JVO = LQ(JVOLUM-LVOLUM(NL))
         NIN = Q(JVO+3)
         DO 488 IN=1,NIN
            JIN = LQ(JVO-IN)
            IQ(JIN) = IBCLR(IQ(JIN),4)
  488    CONTINUE
  489 CONTINUE
*
      IF (NLMANY.GT.0) THEN
         CALL GFCVOL
         NLEVIN = NLEVEL
      ELSEIF (NLEVEL.GT.NLEVIN) THEN
         INGOTO = LINDEX(NLEVEL)
         NL = NLEVIN
         NLEVIN = NLEVEL
         NLEVEL = NL
      ENDIF
*                                                             END GINVOL
  999 IF(JGSTAT.NE.0) CALL GFSTAT(ISAME)
      END
#endif
