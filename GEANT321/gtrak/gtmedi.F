*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if !defined(CERNLIB_OLD)
*CMZ :  3.21/02 29/03/94  15.22.08  by  S.Giani
*-- Author :
      SUBROUTINE GTMEDI (X, NUMED)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Finds in which volume/medium the point X is, and updates the *
C.    *    common /GCVOLU/ and the structure JGPAR accordingly.        *
C.    *                                                                *
C.    *   NUMED returns the tracking medium number, or 0 if point is   *
C.    *         outside the experimental setup.                        *
C.    *                                                                *
C.    *   Note : For INWVOL = 2, INFROM set to a positive number is    *
C.    *      interpreted by GTMEDI as the number IN of the content     *
C.    *      just left by the current track within the mother volume   *
C.    *      where the point X is assumed to be.                       *
C.    *                                                                *
C.    *   Note : INFROM is set correctly by this routine but it is     *
C.    *      used on entrance only in the case GSNEXT has been called  *
C.    *      by the user. In other words the value of INFROM received  *
C.    *      on entrance is not considered necessarily valid. This     *
C.    *      assumption has been made for safety. A wrong value of     *
C.    *      INFROM can cause wrong tracking.                          *
C.    *                                                                *
C.    *   Called by : GTRACK                                           *
C.    *   Authors   : S.Banerjee, R.Brun, F.Bruyant, A.McPherson       *
C.    *               S.Giani.                                         *
C.    *                                                                *
C.    *   Modified by S.Giani (1993) to perform the search according   *
C.    *    to the new 'virtual divisions' algorithm and to build the   *
C.    *    stack of the 'possible overlapping volumes' in the case of  *
C.    *    MANY volumes. Any kind of boolean operation is now possible.*
C.    *    Divisions along arbitrary axis are now possible.            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#include "geant321/gchvir.inc"
#include "geant321/gcvdma.inc"
      COMMON/GCCHAN/LSAMVL
      LOGICAL LSAMVL
C.
      CHARACTER*4 NAME
      DIMENSION  X(*)
      REAL       XC(6), XT(3)
      LOGICAL    BTEST
C.
C.    ------------------------------------------------------------------
*
      nvmany=0
      nfmany=0
      new2fl=0
      neufla=0
      if(raytra.eq.1.)then
        JVO  = LQ(JVOLUM-LVOLUM(NLEVEL))
        if(ingoto.eq.-1.and.q(jvo+3).lt.0)then
          neufla=1
        elseif(ingoto.eq.0)then
              call ggperp(x,veccos,ierr)
              veccos(1)=-veccos(1)
              veccos(2)=-veccos(2)
              veccos(3)=-veccos(3)
              if(ierr.eq.1)then
                veccos(1)=1.
                veccos(2)=0.
                veccos(3)=0.
              endif
        endif
      endif
*
* SECTION I: The /GCVOLU/ table contains the initial guess for a path
*            in the geometry tree on which X may be found.  Look along this
*            path until X is found inside.  This is the starting position.
*            If this is an ONLY volume with no daughters, we are done;
*            otherwise reset search record variables, proceed to section II.
*
* *** Check if point is in current volume
*
      INFR = 0
      INGT = 0
      JVIN = 0
*
* *** LSAMVL is a logical variable that indicates whether we are still
* *** in the current volume or not. It is used in GTRACK to detect
* *** precision problems.
      LSAMVL = .TRUE.
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
      if(raytra.eq.1..and.imyse.eq.1)then
            CALL UHTOC(NAMES(NLEVEL),4,NAME,4)
            CALL GFIND(NAME,'SEEN',ISSEEN)
            if(isseen.eq.-2.or.isseen.eq.-1)goto 189
      endif
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
            LSAMVL = .FALSE.
            NL1 = NLEVEL +1
            LVOLUM(NL1) = IVOT
            NAMES(NL1)  = IQ(JVOLUM+IVOT)
            NUMBER(NL1) = Q(JIN+3)
            LINDEX(NL1) = INGOTO
            LINMX(NL1)  = Q(JVO+3)
            GONLY(NL1)  = Q(JIN+8)
            IF (LQ(LQ(JVOLUM-IVOT)).EQ.0) THEN
               NLDEV(NL1) = NLDEV(NLEVEL)
            ELSE
               NLDEV(NL1) = NL1
            ENDIF
            CALL GTRMUL (GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), Q(JIN+5),
     +                   IROTT, GTRAN(1,NL1), GRMAT(1,NL1))
            NLEVEL = NL1
            XC(1) = XT(1)
            XC(2) = XT(2)
            XC(3) = XT(3)
            JVO = JVOT
            INFROM = 0
            if(raytra.eq.1.)then
              call ggperp(x,veccos,ierr)
              if(ierr.eq.1)then
                veccos(1)=1.
                veccos(2)=0.
                veccos(3)=0.
              endif
            endif
            GO TO 190
         ENDIF
      ENDIF
*
* End of INGOTO processing
*
 189  JPAR = LQ(JGPAR-NLEVEL)
      CALL GINME (XC, Q(JVO+2), Q(JPAR+1), IYES)
      IF (IYES.EQ.0) THEN
*
*  **   Point not in current volume, go up the tree
*
         LSAMVL = .FALSE.
         INGOTO = 0
         IF (NLEVEL.GT.1) THEN
            NLEVEL = NLEVEL -1
            JVO = LQ(JVOLUM-LVOLUM(NLEVEL))
            NIN = Q(JVO+3)
            IF(NIN.GT.0) THEN
               INFROM=LINDEX(NLEVEL+1)
            ELSE
               INFROM=0
            ENDIF
            INFR = INFROM
            GO TO 100
         ELSE
*
*   *      Point is outside setup
*
            NUMED = 0
            GO TO 999
         ENDIF
      ELSE
*
*   *      Point in current volume but not in INGOTO. We block the
*   *      corresponding volume
*
         IF (INGOTO.GT.0) THEN
            INGT = INGOTO
            JIN = LQ(JVO-INGOTO)
            IQ(JIN) = IBSET(IQ(JIN),4)
         ENDIF
      ENDIF
*
*   *      Found a volume up the tree which contains our point. We block
*   *      the branch we came up from.
*
      IF(INFR.GT.0) THEN
         JIN=LQ(JVO-INFR)
         IQ(JIN) = IBSET(IQ(JIN),4)
         JVIN = JIN
      ENDIF
*
*  **   Point is in current volume
*
  190 INGOTO = 0
      NLMIN = NLEVEL
      IF (INWVOL.NE.2) INFROM = 0
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
      if(raytra.eq.1..and.imyse.eq.1)then
            CALL UHTOC(NAMES(NLEVEL),4,NAME,4)
            CALL GFIND(NAME,'SEEN',ISSEEN)
            if(isseen.eq.-2.or.isseen.eq.-1)goto 300
      endif
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
            if(raytra.eq.1..and.neufla.eq.1)then
              neufla=0
              call ggperp(x,veccos,ierr)
              if(ierr.eq.1)then
                veccos(1)=1.
                veccos(2)=0.
                veccos(3)=0.
              endif
            endif
            INFROM = 0
            INFR   = 0
            INGT   = 0
            LSAMVL = .FALSE.
            GO TO 200
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
                  new2fl=0
                  IF (GONLY(NLEVEL).NE.0.) THEN
                    NLMANY = 0
                    nvmany = 0
                    nfmany = 0
                  ENDIF
                  INFROM = 0
                  INGT   = 0
                  INFR   = 0
                  LSAMVL = .FALSE.
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
            INFR = 0
         ENDIF
         IF(INGT.NE.0) THEN
            JIN = LQ(JVO-INGT)
            IQ(JIN) = IBCLR(IQ(JIN),4)
            INGT = 0
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
            nfmany=nvmany+1
         ENDIF
         if(new2fl.eq.0)then
            nvmany=nvmany+1
            manyle(nvmany)=nlevel
            do 401 i = 1,nlevel
              manyna(nvmany,i)=names(i)
              manynu(nvmany,i)=number(i)
 401        continue
         endif
*
*   *   Go up the tree up to a volume with positioned contents
*
         new2fl=-1
  310    INFROM = LINDEX(NLEVEL)
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
*
         ENDIF
C*****  End of Code Expanded From Routine:  GTRNSF
*
         INFR = INFROM
         JIN = LQ(JVO-INFROM)
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
      if(nlmany.eq.0)then
        nvmany=0
        nfmany=0
      endif
      IF (NLMANY.GT.0) CALL GFCVOL
      JVO = LQ(JVOLUM-LVOLUM(NLEVEL))
      IF(JVIN.NE.0) IQ(JVIN) = IBCLR(IQ(JVIN),4)
      NUMED = Q(JVO+4)
*                                                             END GTMEDI
  999 IF(JGSTAT.NE.0) CALL GFSTAT(4)
      END
#endif
