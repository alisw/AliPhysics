*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if !defined(CERNLIB_OLD)
*CMZ :  3.21/02 29/03/94  15.06.41  by  S.Giani
*-- Author :
      SUBROUTINE GNEXT (X, SNEXT, SAFETY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   SUBR. GNEXT (X, SNEXT, SAFETY)                               *
C.    *                                                                *
C.    *   Computes SNEXT and SAFETY                                    *
C.    *     SNEXT  (output) : distance to closest boundary             *
C.    *                      from point X(1-3) along X(4-6)            *
C.    *     SAFETY (output) : shortest distance to any boundary        *
C.    *                                                                *
C.    *   Called by : User                                             *
C.    *   Author   : S.Giani (1993)                                    *
C.    *                                                                *
C.    *   This routine is now based on the new 'virtual divisions'     *
C.    *    algorithm to speed up the tracking.                         *
C.    *   The tracking for MANY volumes is not anymore based on a step *
C.    *    search: it is now based on a search through the list of     *
C.    *    'possible overlapping volumes' built by GTMEDI.             *
C.    *    Boolean operations and divisions along arbitrary axis are   *
C.    *     now supported.                                             *
C.    *    Important notice: in case of MANY volumes, it might happen  *
C.    *     that at the end of GNEXT the volume where X was found to   *
C.    *     be is not anymore the current volume in GCVOLU (being the  *
C.    *     tree ready for the volume where the particles is entering).*
C.    *     When this happens, the variable MYCOUN in the common block *
C.    *     PHYCOU is greater than one; the user can take the proper   *
C.    *     actions (before calling GINVOL, for example) like restoring*
C.    *     the old tree: this is possible because GMEDIA has saved    *
C.    *     the needed informations in MANYLE(NFMANY) (for nlevel), in *
C.    *     MANYNA(NFMANY,i) (for names(i)) and in MANYNU(NFMANY,i)    *
C.    *     (for number(i)); this is sufficient to restore the old tree*
C.    *     by calling GLVOLU.                                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcshno.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#include "geant321/gchvir.inc"
#include "geant321/gcvdma.inc"
      DIMENSION NUMTMP(15),NAMTMP(15)
      dimension iarrin(500),cxm(3),xxm(6)
      REAL    X(6), X0(6), XC(6), XT(6)
      INTEGER IDTYP(3,12)
      LOGICAL BTEST
      SAVE IDTYP
C.
      DATA  IDTYP / 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 2, 3, 1,
     +              2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 4, 3, 1, 1, 1,
     +              2, 3, 1, 2, 3, 1/
C.
C.    ------------------------------------------------------------------
*
* *** Transform current point and direction into local reference system
*
      mycoun=0
      myinfr=0
      newfl=0
      manyfl=0
      tsafet=big
      tsnext=big
 401  IF (GRMAT(10,NLEVEL).EQ.0.) THEN
         DO 19 I = 1,3
            XC(I)   = X(I) -GTRAN(I,NLEVEL)
            XC(I+3) = X(I+3)
   19    CONTINUE
      ELSE
*       (later, code in line)
         CALL GTRNSF (X, GTRAN(1,NLEVEL), GRMAT(1,NLEVEL), XC)
         CALL GROT (X(4), GRMAT(1,NLEVEL), XC(4))
      ENDIF
*
* *** Compute distance to boundaries
*
      SNEXT  = BIG
      SAFETY = BIG
      JVO    = LQ(JVOLUM-LVOLUM(NLEVEL))
      ISH    = Q(JVO+2)
      IF (Q(JVO+3).EQ.0.) GO TO 300
      NIN = Q(JVO+3)
      IF (NIN.LT.0) GO TO 200
*
* *** Case with contents positioned
*
      sneold=SNEXT
      nnn=0
      nflag=0
      mmm=0
      snxtot=0.
 111  if(nin.gt.1)then
        if(nnn.gt.0)goto 112
        clmoth=q(jvirt+4*(LVOLUM(NLEVEL)-1)+3)
        chmoth=q(jvirt+4*(LVOLUM(NLEVEL)-1)+4)
        ndivto=q(jvirt+4*(LVOLUM(NLEVEL)-1)+2)
        iaxis =q(jvirt+4*(LVOLUM(NLEVEL)-1)+1)
        if(iaxis.eq.4)then
         do 1 i=1,6
          xxm(i)=xc(i)
 1       continue
        endif
        divthi=(chmoth-clmoth)/ndivto
        if(iaxis.le.3)then
          cx=xc(iaxis)
          if(xc(iaxis+3).ge.0.)then
            inc=1
          else
            inc=-1
          endif
          xvdiv=((cx-clmoth)*ndivto/(chmoth-clmoth))+1
          ivdiv=xvdiv
          if((xvdiv-ivdiv).lt.0.0001.and.inc.eq.-1)ivdiv=ivdiv-1
          if(ivdiv.lt.1)then
            ivdiv=1
          elseif(ivdiv.gt.ndivto)then
            ivdiv=ndivto
          endif
        else
          call gfcoor(xc,iaxis,cx)
          if(iaxis.eq.4)then
            dr= xc(1)*xc(4)+xc(2)*xc(5)
*            if(dr.eq.0.)print *,'dr.eq.0.'
            if(dr.ge.0.)then
              inc=1
            else
              inc=-1
            endif
          elseif(iaxis.eq.6)then
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
            dfi=xc(1)*xc(5)-xc(2)*xc(4)
            if(dfi.ge.0)then
              inc=1
            else
              inc=-1
            endif
          endif
          xvdiv=((cx-clmoth)*ndivto/(chmoth-clmoth))+1
          ivdiv=xvdiv
          if((xvdiv-ivdiv).lt.0.0001.and.inc.eq.-1)ivdiv=ivdiv-1
          if(ivdiv.lt.1)then
            ivdiv=1
          elseif(ivdiv.gt.ndivto)then
            ivdiv=ndivto
          endif
        endif
        jvdiv=lq(jvirt-LVOLUM(NLEVEL))
 112    iofset=iq(jvdiv+ivdiv)
        jcont2=jvdiv+iofset+1
        ncont=iq(jcont2)
        if(ncont.eq.0)then
          idmi=iq(jcont2+1)
          idma=iq(jcont2+2)
          llflag=0
        elseif(ncont.eq.1)then
          idmi=iq(jcont2+2)
          idma=iq(jcont2+3)
          in=iq(jcont2+1)
        else
          idmi=iq(jcont2+ncont+1)
          idma=iq(jcont2+ncont+2)
          iii=1
          in=iq(jcont2+iii)
        endif
        if(nnn.eq.0)then
         cxold=cx
         if(inc.gt.0)then
          cmin=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
          if(iaxis.ne.6)then
           safety=min(safety,(cxold-cmin))
          else
           safefi=min(90.,(cxold-cmin))
           saferr=sqrt(xc(1)**2+xc(2)**2)
           safe22=saferr*sin(safefi)
           safety=min(safety,safe22)
          endif
         else
          cmax=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
          if(iaxis.ne.6)then
           safety=min(safety,(cmax-cxold))
          else
           safefi=min(90.,(cmax-cxold))
           saferr=sqrt(xc(1)**2+xc(2)**2)
           safe22=saferr*sin(safefi)
           safety=min(safety,safe22)
          endif
         endif
        endif
        if(ncont.eq.0)goto 181
      elseif(nin.eq.1)then
        in=1
      endif
*
  150 if(nin.gt.1.and.ncont.gt.1)then
        in=iq(jcont2+iii)
      endif
      if(nin.gt.0)then
        jin=lq(jvo-in)
        if(.NOT.BTEST(iq(jin),4))then
        else
          goto 171
        endif
      endif
      if(nin.gt.1)then
        llflag=0
        if(mmm.le.500)then
         do 151 ll=1,mmm
          if(iarrin(ll).eq.in)then
            llflag=1
            goto 171
          endif
 151     continue
        endif
        if(llflag.eq.0)then
          mmm=mmm+1
          if(mmm.le.500)then
           iarrin(mmm)=in
          endif
        endif
      endif
      IF(IN.LE.0)GO TO 300
      JIN   = LQ(JVO-IN)
      IVOT  = Q(JIN+2)
      JVOT  = LQ(JVOLUM-IVOT)
      IROTT = Q(JIN+4)
*
      IF (NLEVEL.GE.NLDEV(NLEVEL)) THEN
*       (case with JVOLUM structure locally developed)
         JPAR = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         DO 169 ILEV = NLDEV(NLEVEL), NLEVEL
            IF (IQ(JPAR+1).EQ.0) THEN
               IF (ILEV.EQ.NLEVEL) THEN
                  JPAR = LQ(JPAR-IN)
               ELSE
                  JPAR = LQ(JPAR-LINDEX(ILEV+1))
               ENDIF
               IF (JPAR.EQ.0) GO TO 170
            ELSE IF (IQ(JPAR-3).GT.1) THEN
               JPAR = LQ(JPAR-LINDEX(ILEV+1))
            ELSE
               JPAR = LQ(JPAR-1)
            ENDIF
  169    CONTINUE
         JPAR = JPAR + 5
         GO TO 179
      ENDIF
*     (normal case)
  170 NPAR = Q(JVOT+5)
      IF (NPAR.EQ.0) THEN
         JPAR = JIN +9
      ELSE
         JPAR = JVOT +6
      ENDIF
 179  if((nin.eq.1).or.(nin.gt.1.and.llflag.eq.0))then
*
*   * Compute distance to boundary of current content
*
  180 IF (IROTT.EQ.0) THEN
         DO 189 I = 1,3
            XT(I)   = XC(I) -Q(JIN+4+I)
            XT(I+3) = XC(I+3)
  189    CONTINUE
      ELSE
*       (later, code in line)
         CALL GITRAN (XC, Q(JIN+5), IROTT, XT)
         CALL GRMTD (XC(4), IROTT, XT(4))
      ENDIF
*
      IACT = 2
      ISHT = Q(JVOT+2)
      IF (ISHT.LT.5) THEN
         IF (ISHT.EQ.1) THEN
            CALL GNOBOX (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.2) THEN
            CALL GNOTRA(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.3) THEN
            CALL GNOTRA(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE
            CALL GNOTRP (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ENDIF
      ELSE IF (ISHT.LE.10) THEN
         IF (ISHT.EQ.5) THEN
            CALL GNOTUB(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.6) THEN
            CALL GNOTUB(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.7) THEN
            CALL GNOCON(XT,Q(JPAR+1),IACT,1,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.8) THEN
            CALL GNOCON(XT,Q(JPAR+1),IACT,2,SNEXT,SNXT,SAFE)
         ELSE IF (ISHT.EQ.9) THEN
            CALL GNOSPH (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ELSE
            CALL GNOPAR (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
         ENDIF
      ELSE IF (ISHT.EQ.11) THEN
         CALL GNOPGO (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.12) THEN
         CALL GNOPCO (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.13) THEN
         CALL GNOELT (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.14) THEN
         CALL GNOHYP (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE IF (ISHT.EQ.28) THEN
         CALL GSNGTR (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE,0)
      ELSE IF (ISHT.EQ.NSCTUB) THEN
         CALL GNOCTU (XT,Q(JPAR+1),IACT,SNEXT,SNXT,SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISHT
         STOP
      ENDIF
*
      safe=max(safe,0.)
      if(snxt.le.-prec)snxt=big
      snxt=max(snxt,0.)
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT) THEN
         SNEXT = SNXT
      ENDIF
      endif
 171  if(nin.eq.1)then
        goto 300
      elseif(nin.ge.1.and.ncont.gt.1)then
           iii=iii+1
           if(iii.le.ncont)goto 150
      endif
*
*   *         Compute distance to boundary of current volume
*
 181  if(nnn.eq.0)then
               JPAR = LQ(JGPAR-NLEVEL)
               IACT = 2
               ISH  = Q(JVO+2)
               IF (ISH.LT.5) THEN
                  IF (ISH.EQ.1) THEN
                     CALL GNBOX (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
                  ELSE IF (ISH.EQ.2) THEN
                     CALL GNTRAP (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.3) THEN
                     CALL GNTRAP (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE
                     CALL GNTRP (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
                  ENDIF
               ELSE IF (ISH.LE.10) THEN
                  IF (ISH.EQ.5) THEN
                     CALL GNTUBE (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.6) THEN
                     CALL GNTUBE (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.7) THEN
                     CALL GNCONE (XC, Q(JPAR+1),IACT,1, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.8) THEN
                     CALL GNCONE (XC, Q(JPAR+1),IACT,2, SNEXT,SNXT,SAFE)
                  ELSE IF (ISH.EQ.9) THEN
                     CALL GNSPHR (XC, Q(JPAR+1),IACT, SNEXT, SNXT, SAFE)
                  ELSE
                     CALL GNPARA (XC, Q(JPAR+1),IACT, SNEXT, SNXT, SAFE)
                  ENDIF
               ELSE IF (ISH.EQ.12) THEN
                  CALL GNPCON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.11) THEN
                  CALL GNPGON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.13) THEN
                  CALL GNELTU (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.14) THEN
                  CALL GNHYPE (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE IF (ISH.EQ.28) THEN
                  CALL GSNGTR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,1)
               ELSE IF (ISH.EQ.NSCTUB) THEN
                  CALL GNCTUB (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
               ELSE
                  PRINT *, ' GNEXT : No code for shape ', ISH
                  STOP
               ENDIF
*
               safe=max(safe,0.)
               if(snxt.le.-prec)snxt=big
               snxt=max(snxt,0.)
               IF (SAFE.LT.SAFETY) SAFETY = SAFE
               IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
      endif
      if(iaxis.eq.4)then
       if(idma.eq.ndivto.and.inc.gt.0)goto 400
        cxm(1)=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
        if(idmi.eq.idma)then
          cxm(2)=cxm(1)+divthi
        else
          cxm(2)=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
        endif
        cxm(3)=20000.
        call gntube(xxm,cxm,3,1,SNEXT,snxnew,safe)
        snxnew=snxnew+.004
        snxtot=snxtot+snxnew
        if(snxtot.lt.SNEXT)then
          xxm(1)=xxm(1)+snxnew*xxm(4)
          xxm(2)=xxm(2)+snxnew*xxm(5)
          xxm(3)=xxm(3)+snxnew*xxm(6)
          call gfcoor(xxm,iaxis,cxnew)
          xevdiv=((cxnew-clmoth)*ndivto/(chmoth-clmoth))+1
          ivdiv=xevdiv
          dr= xxm(1)*xxm(4)+xxm(2)*xxm(5)
*          if(dr.eq.0.)print *,'dr.eq.0.'
          if(dr.ge.0.)then
              inc=1
          else
              inc=-1
          endif
          if((xevdiv-ivdiv).lt.0.0001.and.inc.eq.-1)ivdiv=ivdiv-1
          if(ivdiv.lt.1)then
              ivdiv=1
          elseif(ivdiv.gt.ndivto)then
              ivdiv=ndivto
          endif
          nnn=nnn+1
          goto 111
        else
          if(inc.gt.0)then
           cmax=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
           safety=min(safety,(cmax-cxold))
          else
           cmin=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
           safety=min(safety,(cxold-cmin))
          endif
          goto 400
        endif
      endif
          if(nnn.ne.0.and.SNEXT.eq.sneold)goto 199
               x0(1) = xc(1) + SNEXT*xc(4)
               x0(2) = xc(2) + SNEXT*xc(5)
               x0(3) = xc(3) + SNEXT*xc(6)
               x0(4) = xc(4)
               x0(5) = xc(5)
               x0(6) = xc(6)
          if(iaxis.le.3)then
            cx=x0(iaxis)
            xevdiv=((cx-clmoth)*ndivto/(chmoth-clmoth))+1
            ievdiv=xevdiv
            if((xevdiv-ievdiv).lt.0.0001.and.inc.eq.-1)ievdiv=ievdiv-1
            if(ievdiv.lt.1)then
              ievdiv=1
            elseif(ievdiv.gt.ndivto)then
              ievdiv=ndivto
            endif
          else
            call gfcoor(x0,iaxis,cx)
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
            xevdiv=((cx-clmoth)*ndivto/(chmoth-clmoth))+1
            ievdiv=xevdiv
            if((xevdiv-ievdiv).lt.0.0001.and.inc.eq.-1)ievdiv=ievdiv-1
            if(ievdiv.lt.1)then
              ievdiv=1
            elseif(ievdiv.gt.ndivto)then
              ievdiv=ndivto
            endif
          endif
 199      if(ievdiv.ge.idmi.and.ievdiv.le.idma)then
            if(inc.gt.0)then
             cmax=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
             if(iaxis.ne.6)then
              safety=min(safety,(cmax-cxold))
             else
              safefi=min(90.,(cmax-cxold))
              safe22=saferr*sin(safefi)
              safety=min(safety,safe22)
             endif
            else
             cmin=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
             if(iaxis.ne.6)then
              safety=min(safety,(cxold-cmin))
             else
              safefi=min(90.,(cxold-cmin))
              safe22=saferr*sin(safefi)
              safety=min(safety,safe22)
             endif
            endif
            goto 400
          endif
          if(iaxis.eq.6.or.iaxis.le.3)then
           if(ievdiv.lt.idmi.and.inc.gt.0)then
            if(nnn.eq.0.and.iaxis.eq.6)nflag=1
            if(nflag.eq.0)then
*             print *,'ievdiv=',ievdiv,' ;idmi=',idmi,' inc.gt.0'
*             print *,isht,'=isht; ',iaxis,'=iaxis; ',ish,'=ish;'
             if(iaxis.le.3)then
               cmax=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
               safety=min(safety,abs(cmax-cxold))
             elseif(iaxis.eq.6)then
               cmax=(clmoth+(idma-1)*(chmoth-clmoth)/ndivto)+divthi
               safefi=min(90.,(cmax-cxold))
               safe22=saferr*sin(safefi)
               safety=min(safety,safe22)
             endif
             goto 400
            endif
           elseif(ievdiv.gt.idma.and.inc.lt.0)then
            if(nnn.eq.0.and.iaxis.eq.6)nflag=1
            if(nflag.eq.0)then
*             print *,'ievdiv=',ievdiv,' ;idma=',idma,' inc.lt.0'
*             print *,isht,'=isht; ',iaxis,'=iaxis; ',ish,'=ish;'
             if(iaxis.le.3)then
               cmin=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
               safety=min(safety,abs(cxold-cmin))
             elseif(iaxis.eq.6)then
               cmin=clmoth+(idmi-1)*(chmoth-clmoth)/ndivto
               safefi=min(90.,(cxold-cmin))
               safe22=saferr*sin(safefi)
               safety=min(safety,safe22)
             endif
             goto 400
            endif
           endif
          endif
          nnn=nnn+1
          sneold=SNEXT
          if(inc.gt.0)then
            if(iaxis.eq.6)then
             if(idma.eq.ndivto.and.(chmoth-clmoth).eq.360.)then
               ivdiv=1
             else
               ivdiv=idma+1
             endif
            else
             ivdiv=idma+1
            endif
          else
            if(iaxis.eq.6)then
             if(idmi.eq.1.and.(chmoth-clmoth).eq.360.)then
               ivdiv=ndivto
             else
               ivdiv=idmi-1
             endif
            else
             ivdiv=idmi-1
            endif
          endif
          goto 111
*
* ***    Case of volume incompletely divided
*
  200 JDIV  = LQ(JVO-1)
      IAXIS = Q(JDIV+1)
      IVOT  = Q(JDIV+2)
      JVOT  = LQ(JVOLUM-IVOT)
      ISHT  = Q(JVOT+2)
*
*  ** Get the division parameters
*
      IF (NLEVEL.LT.NLDEV(NLEVEL)) THEN
         JPARM = 0
      ELSE
*        (case with JVOLUM structure locally developed)
         JPARM = LQ(LQ(JVOLUM-LVOLUM(NLDEV(NLEVEL))))
         IF (NLEVEL.EQ.NLDEV(NLEVEL)) GO TO 215
         DO 210 ILEV = NLDEV(NLEVEL), NLEVEL-1
            IF (IQ(JPARM+1).EQ.0) THEN
               JPARM = LQ(JPARM-LINDEX(ILEV+1))
               IF (JPARM.EQ.0) GO TO 215
            ELSE IF (IQ(JPARM-3).GT.1) THEN
               JPARM = LQ(JPARM-LINDEX(ILEV+1))
            ELSE
               JPARM = LQ(JPARM-1)
            ENDIF
            IF (ILEV.EQ.NLEVEL-1) THEN
               NDIV = IQ(JPARM+1)
               ORIG =  Q(JPARM+2)
               SDIV =  Q(JPARM+3)
            ENDIF
  210    CONTINUE
         GO TO 220
      ENDIF
*     (normal case)
  215 NDIV = Q(JDIV+3)
      ORIG = Q(JDIV+4)
      SDIV = Q(JDIV+5)
*
*  ** Look at the first and the last divisions only
*
  220 IDT  = IDTYP(IAXIS, ISH)
      IF (IDT.EQ.1) THEN
         IN2 = 0
         IF (XC(IAXIS).LT.ORIG) THEN
            IN  = 1
         ELSE
            IN  = NDIV
         ENDIF
      ELSE IF (IDT.EQ.2) THEN
         R   = XC(1)**2 + XC(2)**2
         IF (ISH.EQ.9) R = R + XC(3)**2
         R   = SQRT(R)
         IN2 = 0
         IF (ISH.EQ.5.OR.ISH.EQ.6.OR.ISH.EQ.9) THEN
            IF (R.LT.ORIG) THEN
               IN  = 1
            ELSE
               IN  = NDIV
            ENDIF
         ELSE
            PRINT *, ' GNEXT : Partially divided ',ISH,IAXIS
            IN  = 1
            IF (NDIV.GT.1) IN2 = NDIV
         ENDIF
      ELSE IF (IDT.EQ.4) THEN
         IN2 = 0
         RXY = XC(1)**2 + XC(2)**2
         RXY = SQRT(RXY)
         IF (XC(3).NE.0.0) THEN
            THET = RADDEG * ATAN (RXY/XC(3))
            IF (THET.LT.0.0) THET = THET + 180.0
         ELSE
            THET = 90.
         ENDIF
         IF (THET.LE.ORIG) THEN
            IN  = 1
         ELSE
            IN  = NDIV
         ENDIF
      ELSE
         PRINT *, ' GNEXT : Partially divided ',ISH,IAXIS
         IN2 = 0
         IF (ISH.EQ.5.OR.ISH.EQ.7) THEN
            IN  = 1
            IF (NDIV.GT.1) IN2 = NDIV
         ELSE
            IF (XC(1).NE.0.0.OR.XC(2).NE.0.0) THEN
               PHI = RADDEG * ATAN2 (XC(2), XC(1))
            ELSE
               PHI = 0.0
            ENDIF
            IF (ISH.EQ.6.OR.ISH.EQ.8) THEN
               IF (PHI.LT.ORIG) THEN
                  IN  = 1
               ELSE
                  IN  = NDIV
               ENDIF
            ELSE
               IN  = 1
               IF (NDIV.GT.1) IN2 = NDIV
            ENDIF
         ENDIF
      ENDIF
*
  225 IF (IDT.EQ.1) THEN
         DO 231 I = 1, 3
            X0(I) = 0.0
  231    CONTINUE
         X0(IAXIS) = ORIG + (IN - 0.5) * SDIV
         IF (ISH.EQ.4.OR.(ISH.EQ.10.AND.IAXIS.NE.1)) THEN
            CALL GCENT (IAXIS, X0)
         ENDIF
         DO 232 I = 1, 3
            XT(I)   = XC(I) - X0(I)
            XT(I+3) = XC(I+3)
  232    CONTINUE
      ELSE IF (IDT.EQ.3) THEN
         PH0  = DEGRAD * (ORIG + (IN - 0.5) * SDIV)
         CPHR = COS(PH0)
         SPHR = SIN(PH0)
         DO 233 I = 1, 4, 3
            XT(I)   = XC(I)*CPHR + XC(I+1)*SPHR
            XT(I+1) = XC(I+1)*CPHR - XC(I)*SPHR
            XT(I+2) = XC(I+2)
  233    CONTINUE
      ELSE
         DO 234 I = 1, 6
            XT(I) = XC(I)
  234    CONTINUE
      ENDIF
*
      IF (JPARM.NE.0) THEN
         IF (IQ(JPARM-3).GT.1) THEN
            JPAR = LQ(JPARM-IN)
         ELSE
            JPAR = LQ(JPARM-1)
         ENDIF
         JPAR = JPAR + 5
      ELSE
         JPAR = JVOT + 6
      ENDIF
*
      IACT = 2
      IF (ISHT.LT.5) THEN
         IF (ISHT.EQ.1) THEN
            CALL GNOBOX (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.2) THEN
            CALL GNOTRA (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.3) THEN
            CALL GNOTRA (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNOTRP (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISHT.LE.10) THEN
         IF (ISHT.EQ.5) THEN
            CALL GNOTUB (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.6) THEN
            CALL GNOTUB (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.7) THEN
            CALL GNOCON (XT, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.8) THEN
            CALL GNOCON (XT, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISHT.EQ.9) THEN
            CALL GNOSPH (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNOPAR (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISHT.EQ.11) THEN
         CALL GNOPGO (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.12) THEN
         CALL GNOPCO (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.13) THEN
         CALL GNOELT (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.14) THEN
         CALL GNOHYP (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISHT.EQ.28) THEN
         CALL GSNGTR (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,0)
      ELSE IF (ISHT.EQ.NSCTUB) THEN
         CALL GNOCTU (XT, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISHT
         STOP
      ENDIF
*
      safe=max(safe,0.)
      if(snxt.le.-prec)snxt=big
      snxt=max(snxt,0.)
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
*
      IF (IN2.NE.0) THEN
         IF (IN2.NE.IN) THEN
            IN  = IN2
            GO TO 225
         ENDIF
      ENDIF
*
* ***  Calculate SNEXT and SAFETY with respect to the Mother
* ***            SAFETY only for concave volumes if finite SNEXT
* ***            has been found with respect to one of its contents
*
  300 IACT = 2
      IF (SNEXT.LT.0.9*BIG) THEN
         IF (.NOT.BTEST(IQ(JVO),2)) IACT = 0
      ENDIF
      if(nin.eq.1.and.SNEXT.LT.0.9*BIG)then
        if(q(jin+8).eq.0.)iact=2
      endif
      JPAR = LQ(JGPAR-NLEVEL)
      IF (ISH.LT.5) THEN
         IF (ISH.EQ.1) THEN
            CALL GNBOX (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE )
         ELSE IF (ISH.EQ.2) THEN
            CALL GNTRAP (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.3) THEN
            CALL GNTRAP (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNTRP (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISH.LE.10) THEN
         IF (ISH.EQ.5) THEN
            CALL GNTUBE (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.6) THEN
            CALL GNTUBE (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.7) THEN
            CALL GNCONE (XC, Q(JPAR+1), IACT, 1, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.8) THEN
            CALL GNCONE (XC, Q(JPAR+1), IACT, 2, SNEXT, SNXT, SAFE)
         ELSE IF (ISH.EQ.9) THEN
            CALL GNSPHR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ELSE
            CALL GNPARA (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
         ENDIF
      ELSE IF (ISH.EQ.12) THEN
         CALL GNPCON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.11) THEN
         CALL GNPGON (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.13) THEN
         CALL GNELTU (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.14) THEN
         CALL GNHYPE (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE IF (ISH.EQ.28) THEN
         CALL GSNGTR (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE,1)
      ELSE IF (ISH.EQ.NSCTUB) THEN
         CALL GNCTUB (XC, Q(JPAR+1), IACT, SNEXT, SNXT, SAFE)
      ELSE
         PRINT *, ' GNEXT : No code for shape ', ISH
         STOP
      ENDIF
*
      safe=max(safe,0.)
      if(snxt.le.-prec)snxt=big
      snxt=max(snxt,0.)
      IF (SAFE.LT.SAFETY) SAFETY = SAFE
      IF (SNXT.LT.SNEXT)  SNEXT  = SNXT
*
 400  if(myinfr.gt.0)then
        jin=lq(jvo-myinfr)
        iq(jin)=ibclr(iq(jin),4)
        myinfr=0
      endif
      if(gonly(nlevel).eq.0..or.nvmany.ne.0) THEN
         if(safety.lt.tsafet)tsafet=safety
         if(snext.lt.tsnext)then
          mycoun=mycoun+1
          tsnext=snext
          call gscvol
         endif
         if(gonly(nlevel).eq.0.)then
 404       continue
           if(gonly(nlevel-1).eq.0..or.newfl.eq.0)then
             if(gonly(nlevel-1).ne.0.)newfl=1
             nlevel=nlevel-1
             jvo=lq(jvolum-lvolum(nlevel))
             nin=q(jvo+3)
             if(nin.lt.0)goto 404
             myinfr=lindex(nlevel+1)
             jin=lq(jvo-myinfr)
             iq(jin)=ibset(iq(jin),4)
             ignext=0
             goto 401
           endif
         endif
 403   continue
       if(manyfl.lt.nvmany)then
         manyfl=manyfl+1
         if(manyfl.eq.nfmany)goto 403
         levtmp=manyle(manyfl)
         do 402 i=1,levtmp
          namtmp(i)=manyna(manyfl,i)
          numtmp(i)=manynu(manyfl,i)
 402     continue
         call glvolu(levtmp,namtmp,numtmp,ier)
         if(ier.ne.0)print *,'Fatal error in GLVOLU'
         ignext=0
         goto 401
       endif
       if(tsafet.le.safety)safety=tsafet
       if(tsnext.le.snext)then
         snext=tsnext
         call gfcvol
       endif
      endif
*                                                              END GNEXT
  999 END
 
#endif
