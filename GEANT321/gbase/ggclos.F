*
* $Id$
*
* $Log$
* Revision 1.2  1997/11/14 17:44:00  mclareni
* Make sure the maximum angle is greater than the minimun
*
* Revision 1.1.1.1  1995/10/24 10:20:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if !defined(CERNLIB_OLD)
*CMZ :  3.21/04 13/12/94  15.29.27  by  S.Giani
*-- Author :
      SUBROUTINE GGCLOS
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    Closes off the geometry setting.                            *
C.    *    Initializes the search list for the contents of each        *
C.    *    volume following the order they have been positioned, and   *
C.    *    inserting the content '0' when a call to GSNEXT (-1) has    *
C.    *    been required by the user.                                  *
C.    *    Performs the development of the JVOLUM structure for all    *
C.    *    volumes with variable parameters, by calling GGDVLP.        *
C.    *    Interprets the user calls to GSORD, through GGORD.          *
C.    *    Computes and stores in a bank (next to JVOLUM mother bank)  *
C.    *    the number of levels in the geometrical tree and the        *
C.    *    maximum number of contents per level, by calling GGNLEV.    *
C.    *    Sets status bit for CONCAVE volumes, through GGCAVE.        *
C.    *    Completes the JSET structure with the list of volume names  *
C.    *    which identify uniquely a given physical detector, the      *
C.    *    list of bit numbers to pack the corresponding volume copy   *
C.    *    numbers, and the generic path(s) in the JVOLUM tree,        *
C.    *    through the routine GHCLOS.                                 *
C.    *                                                                *
C.    *    Called by : <USER>                                          *
C.    *    Authors   : R.Brun, F.Bruyant, S.Giani  *********           *
C.    *                                                                *
C.    *    Modified by S.Giani for automatic initialization of the new *
C.    *    tracking based on virtual divisions (1993).                 *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gclist.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcopti.inc"
#include "geant321/gchvir.inc"
      CHARACTER*4 NAME
      LOGICAL BTEST
C.
C.    ------------------------------------------------------------------
      dimension dx(3),tmpmax(7),ndivto(7),qualit(7),ivoaxi(7)
      data jfirst/0/
      save jfirst
      COMMON /QUEST/ IQUEST(100)
      COMMON/GCDINA/jphi2,jclow,jchig,jbuff
*
* *** Stop the run in case of serious anomaly during initialization
*
      IF (IEORUN.NE.0) THEN
         WRITE (CHMAIL, 1001)
         CALL GMAIL (0, 0)
         STOP
      ENDIF
*
      IF (NVOLUM.LE.0) THEN
         WRITE (CHMAIL, 1002) NVOLUM
         CALL GMAIL (0, 0)
         GO TO 999
      ENDIF
*
      NPUSH = NVOLUM -IQ(JVOLUM-2)
      CALL MZPUSH (IXCONS, JVOLUM, NPUSH, NPUSH,'I')
*
* *** Loop over volumes, create default JNear banks as relevant,
*      and release unused bank space
*
      IDO = 0
      DO 80 IVO = 1,NVOLUM
         JVO = LQ(JVOLUM-IVO)
*
* *** Check if Tracking medium has been defined
*
         NMED=Q(JVO+4)
         IF(NMED.LE.0.OR.NMED.GT.IQ(JTMED-2))THEN
            WRITE(CHMAIL,1003)IQ(JVOLUM+IVO)
            CALL GMAIL (0, 0)
         ELSE
            IF(LQ(JTMED-NMED).EQ.0)THEN
               WRITE(CHMAIL,1003)IQ(JVOLUM+IVO)
               CALL GMAIL (0, 0)
            ENDIF
         ENDIF
         IF (BTEST(IQ(JVO),0)) GO TO 80
         IDO = 1
         IQ(JVO) = IBSET(IQ(JVO),0)
         NINL  = IQ(JVO-2)
         NIN   = Q(JVO+3)
         NUSED = IABS(NIN)
         IF (NIN.GT.0) THEN
*           reserve enough additional space for sorted volumes
            IF(NIN.LE.1.OR.NIN.GT.500.OR.IOPTIM.LT.0)THEN
              NUSED=NUSED+1
            ELSE
              NUSED=NUSED+2
            ENDIF
         ENDIF
*
         NPUSH = NUSED -NINL
         DO 90 IN=NINL,NUSED+1,-1
            JIN = LQ(JVO-IN)
            IF(JIN.GT.0) THEN
               CALL MZDROP(IXCONS,JIN,'L')
            ENDIF
  90     CONTINUE
         CALL MZPUSH (IXCONS, JVO, NPUSH, 0, 'I')
         IF (NIN.LE.0) GO TO 80
*
         IF(BTEST(IQ(JVO),3)) THEN
            IZERO=1
         ELSE
            IZERO=0
         ENDIF
         NEL = NIN +IZERO
         JN = LQ(JVO-NIN-1)
         IF(JN.EQ.0) THEN
            CALL MZBOOK (IXCONS,JN,JVO,-NIN-1,'VONE',0,0,NEL+1,2,0)
         ENDIF
         IQ(JN-5) = IVO
         IQ(JN+1) = NEL
         JN = JN +1
         DO 29 I = 1,NIN
            IQ(JN+IZERO+I) = I
   29    CONTINUE
         IF (IZERO.NE.0) IQ(JN+1) = 0
*
   80 CONTINUE
*
      IF (IDO.NE.0) THEN
*
* ***    Perform development of JVOLUM structure where necessary
*
         CALL GGDVLP
*
* ***    Fill GSORD ordering banks if required
*
* Modified by S.Egli to allow GGORDQ to find the optimum sorting for
* all volumes
*
         IF(IOPTIM.GE.1)THEN
            WRITE(6,'(A)')' GGCLOS: Start automatic volume ordering:'
         ENDIF
         DO 91 IVO = 1,NVOLUM
            JVO = LQ(JVOLUM-IVO)
            NIN = Q(JVO+3)
            ISEARC=Q(JVO+1)
            IF(ISEARC.GT.0) GO TO 91
*           check if sorting not possible or not wanted
            IF(NIN.LE.1.OR.NIN.GT.500.OR.IOPTIM.LT.0)THEN
               Q(JVO+1)=0.
               IF(NIN.GT.500.AND.IOPTIM.GE.1)THEN
                 CALL UHTOC(IQ(JVOLUM+IVO),4,NAME,4)
                 WRITE (CHMAIL,1004) NAME,NIN
                 CALL  GMAIL (0, 0)
               ENDIF
            ELSEIF(IOPTIM.EQ.0)THEN
               IF(ISEARC.LT.0)CALL GGORD (IVO)
            ELSEIF(IOPTIM.EQ.1)THEN
               IF(ISEARC.EQ.0) THEN
                  CALL GGORDQ(IVO)
               ELSE
                  CALL GGORD (IVO)
               END IF
            ELSE
               CALL GGORDQ(IVO)
            ENDIF
   91    CONTINUE
*
* ***    Set status bit for concave volumes
*
         CALL GGCAVE
*
* ***    Compute maximum number of levels and of contents per level
*
         CALL GGNLEV
*
      ENDIF
*
********************************************************************************
*
      if(jfirst.eq.0)then
        jfirst=1
        call mzlink(ixcons,'/GCHVIR/',jvirt,jvdiv,jcont)
        call mzlink(ixstor,'/GCDINA/',jphi2,jbuff,jphi2)
      endif
      jflag=0
      nwjvdi=0
      jphi2=0
      jclow=0
      jchig=0
      jbuff=0
      if(jvirt.ne.0)call mzdrop(ixcons,jvirt,' ')
      nwjvir=5*nvolum+20
      call mzneed(ixcons,nwjvir,'G')
      if(iquest(11).lt.0)then
          print *,'No space for jvirt bank'
      else
          call mzbook(ixcons,jvirt,jvirt,1,'VIRT',nvolum,nvolum,
     +              4*nvolum+20,0,0)
      endif
      dx(1)=0.
      dx(2)=0.
      dx(3)=0.
      ndivst=0
      ndioff=0
      ninmax=0
      do 101 ivo=1,nvolum
        jvo=lq(jvolum-ivo)
        call uhtoc(iq(jvolum+ivo),4,NAME,4)
*         print *,'VOLUME ',NAME
*         print *,' '
        nin=q(jvo+3)
        isearc=q(jvo+1)
*        if(nin.eq.0)then
*          print *,'No daughters.'
*        elseif(nin.lt.0)then
*          print *,'Divided volume.'
*        elseif(nin.le.1)then
*          print *,'Only 1 daughter.'
*        endif
 1      continue
        if(nin.gt.1)then
         if(jflag.eq.0)then
          if(iswit(9).eq.12345)then
            print *,'VOLUME ',NAME
            print *,' '
          endif
         endif
         if(jflag.eq.1)then
           q(jvirt+4*(ivo-1)+1)=itmpq
           iaxlo=itmpq
           iaxhi=itmpq
         else
           iaxlo=1
           iaxhi=7
         endif
          if(nin.gt.ninmax)then
           if(jphi2.ne.0)call mzdrop(ixstor,jphi2,' ')
           if(jclow.ne.0)call mzdrop(ixstor,jclow,' ')
           if(jchig.ne.0)call mzdrop(ixstor,jchig,' ')
           call mzbook(ixstor,jphi2,jphi2,2,'PHI2',0,0,
     +                 nin+20,2,-1)
           call mzbook(ixstor,jclow,jclow,2,'CLOW',0,0,
     +                 nin+20,3,-1)
           call mzbook(ixstor,jchig,jchig,2,'CHIG',0,0,
     +                 nin+20,3,-1)
           if(jflag.eq.1)then
             if(jbuff.ne.0)call mzdrop(ixstor,jbuff,' ')
             call mzbook(ixstor,jbuff,jbuff,2,'BUFF',0,0,
     +                   nin+20,2,-1)
           endif
          endif
         do 110 iaxis=iaxlo,iaxhi
          myphif=0
*          print *,'Quality search for axis ',iaxis
          ish=q(jvo+2)
          if(iaxis.le.3)then
            call gvdcar(iaxis,ish,0,q(jvo+7),clmoth,chmoth,ierr)
            if(ierr.eq.1.or.(chmoth.le.clmoth))then
*              print *,'Not convenient: abandoned!',iaxis
*              print *,' '
              qualit(iaxis)=10000
              goto 110
            endif
          elseif(iaxis.le.5)then
            call gvdrad(iaxis,ish,0,dx,q(jvo+7),clmoth,chmoth,ierr)
             if(iaxis.eq.5)ierr=1
            if(ierr.eq.1.or.(chmoth.le.clmoth))then
*              print *,'Not convenient: abandoned!',iaxis
*              print *,' '
              qualit(iaxis)=10000
              goto 110
            endif
          elseif(iaxis.eq.6)then
            call gvdphi(ish,0,dx,q(jvo+7),clmoth,chmoth,ierr)
            if(ierr.eq.1.or.(chmoth.le.clmoth))then
*              print *,'Not convenient: abandoned!',iaxis
*              print *,' '
              qualit(iaxis)=10000
              goto 110
            elseif((chmoth-clmoth).gt.360..or.chmoth.gt.360)then
              print *,'(chmoth-clmoth).gt.360.or.chmoth.gt.360'
            elseif((chmoth-clmoth).eq.360.)then
              myphif=1
            endif
          elseif(iaxis.eq.7)then
            call gvdthe(ish,0,dx,q(jvo+7),clmoth,chmoth,ierr)
            ierr=1
            if(ierr.eq.1.or.(chmoth.le.clmoth))then
*              print *,'Not convenient: abandoned!',iaxis
*              print *,' '
              qualit(iaxis)=10000
              goto 110
            endif
          endif
          if(jflag.eq.1)then
            q(jvirt+4*(ivo-1)+3)=clmoth
            q(jvirt+4*(ivo-1)+4)=chmoth
          endif
          thimot=abs(chmoth-clmoth)
          thimin=100000.
          do 102 in=1,nin
            iq(jphi2+in)=0
            jin=lq(jvo-in)
            call gvdlim(jvo,in,iaxis,clow,chigh,ierr)
            if(ierr.eq.1.or.(chigh.le.clow))then
*              if(ierr.eq.0)print *,'Error in gvdlim: corrected',iaxis
              clow=clmoth
              chigh=chmoth
            elseif(myphif.eq.1)then
              clowm=clow
              chighm=chigh
              sg=sign(1.0,clow)
              clow=mod(abs(clow),360.0)
              if(chigh.ne.360.0)then
               if(sg.le.0.0)clow=360.-clow
               sg=sign(1.0,chigh)
               chigh=mod(abs(chigh),360.0)
               if(sg.le.0.0)chigh=360.-chigh
              endif
              if(chigh.lt.clow)then
                chightf = clow
                clow = chigh
                chigh = chightf
                iq(jphi2+in)=1
              endif
            elseif(iaxis.eq.6.and.myphif.eq.0)then
             if((chigh-chmoth).gt..01.or.(clmoth-clow).gt..01)then
               if(clmoth.lt.0..and.clow.gt.0.)then
                clow=clow-360.
                chigh=chigh-360.
                if((chigh-chmoth).gt..01)then
                  chigh=chmoth
                  if(chigh.le.clow)clow=clmoth
                elseif((clmoth-clow).gt..01)then
                  clow=clmoth
                  if(clow.ge.chigh)chigh=chmoth
                endif
               elseif(chigh.lt.0..and.chmoth.gt.0.)then
                clow=clow+360.
                chigh=chigh+360.
                if((chigh-chmoth).gt..01)then
                  chigh=chmoth
                  if(chigh.le.clow)clow=clmoth
                elseif((clmoth-clow).gt..01)then
                  clow=clmoth
                  if(clow.ge.chigh)chigh=chmoth
                endif
               endif
             endif
            endif
            if((chigh-chmoth).gt..01)then
*** ONLY FOR DEBUG
*              print *,'iaxis =',iaxis,'protuding daughter',in
**              print *,'myphif =',myphif,'myphi2 =',iq(jphi2+in)
*              print *,'Dhigh=',chigh-chmoth
*              print *,'Dlow=',clmoth-clow
***
              chigh=chmoth
              if(chigh.le.clow)clow=clmoth
            elseif((clmoth-clow).gt..01)then
*** ONLY FOR DEBUG
*              print *,'iaxis =',iaxis,'protuding daughter',in
**              print *,'myphif =',myphif,'myphi2 =',iq(jphi2+in)
*              print *,'Dhigh=',chigh-chmoth
*              print *,'Dlow=',clmoth-clow
***
              clow=clmoth
              if(clow.ge.chigh)chigh=chmoth
            endif
            q(jclow+in)=clow
            q(jchig+in)=chigh
            if(iq(jphi2+in).eq.0)then
             tmpthi=abs(chigh-clow)
            else
             tmpthi=abs(chighm-clowm)
            endif
            if(thimin.gt.tmpthi)thimin=tmpthi
 102      continue
          if((thimin-thimot).gt.1)then
*            print *,'thimin.gt.thimot',thimin-thimot,'iax=',iaxis
            qualit(iaxis)=10000
            goto 110
          endif
          if(thimin.lt.0.04)thimin=0.04
          tmpndi=2.*thimot/thimin
          nditmp=tmpndi+1
*****          print *,nditmp,' divisions asked for ',nin,' daughters.'
*****           if(nditmp.lt.nin)then
*****             nditmp=nin
*****             print *,'Number of divisions corrected to be = ',nin
*****           endif
*****          if(nditmp.gt.1000.)print *,'1000 divisions are enough.'
          ndivto(iaxis)=min(nditmp,1000)
          if(jflag.eq.1)then
            q(jvirt+4*(ivo-1)+2)=ndivto(iaxis)
           jvdiv=lq(jvirt-ivo)
           if(jvdiv.ne.0)call mzdrop(ixcons,jvdiv,' ')
           nwvili=ndivto(iaxis)+ivoaxi(itmpq)+11
           nwjvdi=nwjvdi+nwvili
           call mzneed(ixcons,nwvili,'G')
           if(iquest(11).lt.0)then
            print *,'No space for jvdiv bank',ivo
           else
            call mzbook(ixcons,jvdiv,jvirt,-ivo,'VLIST',0,0,
     +                  nwvili,2,0)
           endif
          endif
          thisli=thimot/ndivto(iaxis)
          clslic=clmoth
          chslic=clmoth+thisli
          avelis=0.
          aveave=0.
          avesta=0.
          ii=0
          tmpmax(iaxis)=0.
          import=0
          if(jflag.eq.1)ioff=ndivto(iaxis)
          do 103 i=1,ndivto(iaxis)
            j=1
            do 104 in=1,nin
             if(iq(jphi2+in).eq.0)then
              if(q(jchig+in).ge.clslic.and.
     +           q(jclow+in).le.chslic)then
                j=j+1
                if(jflag.eq.1)then
                  iq(jbuff+j)=in
                endif
              endif
             else
              if(q(jchig+in).ge.clslic.or.
     +           q(jclow+in).le.chslic)then
                j=j+1
                if(jflag.eq.1)then
                  iq(jbuff+j)=in
                endif
              endif
             endif
 104        continue
            inbuf1=j-1
            if(jflag.eq.1)then
             if(i.gt.1.and.iq(jbuff+1).eq.(j-1))then
               if(j-1.eq.0)then
                 import=1
               elseif(j-1.eq.1)then
                 if(iq(jbuff+2).eq.iq(jvdiv+ioff-nposti+2))then
                   import=1
                 else
                   import=0
                 endif
               else
                import=1
                do 234 ijk=2,nposti-2
                 do 432 kji=2,nposti-2
                  if(iq(jbuff+ijk).eq.iq(jvdiv+ioff-nposti+kji))then
                    goto 234
                  endif
 432             continue
                 import=0
                 goto 235
 234            continue
 235            continue
               endif
               if(import.eq.1)then
                 iq(jvdiv+ioff-nposti+nposti)=i
                 iq(jvdiv+i)=ioff-nposti
                 goto 145
               endif
             else
               import=0
             endif
             iq(jbuff+1)=j-1
             nposti=j+2
             iq(jbuff+j+1)=i
             iq(jbuff+j+2)=i
             iq(jvdiv+i)=ioff
             do 144 m=1,nposti
               iq(jvdiv+ioff+m)=iq(jbuff+m)
 144         continue
             ioff=ioff+nposti
            else
             aveinc=j+2
             avesta=avesta+aveinc
            endif
 145        continue
            if(inbuf1.gt.tmpmax(iaxis))then
             tmpmax(iaxis)=inbuf1
            endif
            if(inbuf1.ne.0.)ii=ii+1
            avelis=avelis+inbuf1
            clslic=chslic
            chslic=clslic+thisli
 103      continue
          if(jflag.eq.1)then
             ndioff=ndioff+ioff
             if(iswit(9).eq.12345)then
              print *,'words booked =',nwvili,'; words used =',ioff
              print *,' '
          endif
*** ONLY FOR DEBUG
**             mymyof=0
**             do 2 mm=1,ndivto(iaxis)
**               myoff=iq(jvdiv+mm)
**              if(myoff.ne.mymyof)then
**               if(iq(jvdiv+myoff+1).eq.0)then
**                 print *,'Lower div =',iq(jvdiv+myoff+2)
**                 print *,'Upper div =',iq(jvdiv+myoff+3)
**               elseif(iq(jvdiv+myoff+1).eq.1)then
**                 print *,'Lower div =',iq(jvdiv+myoff+3)
**                 print *,'Upper div =',iq(jvdiv+myoff+4)
**               endif
**             endif
**               mymyof=iq(jvdiv+mm)
** 2           continue
***
          endif
          if(ii.eq.0)then
            print *,iaxis,'=iax: not filled divisions: error!'
            print *,' '
            aveave=10000
            avelis=10000
            goto 105
          endif
          if(jflag.eq.0)then
            ivoaxi(iaxis)=avesta
          endif
          aveave=avelis/ndivto(iaxis)
          avelis=avelis/ii
 105      continue
          qualit(iaxis)=avelis
*** ONLY FOR DEBUG
**      print *,'Max n. of objects per div = ',tmpmax(iaxis)
**      print *,'Aver. n. of obj. per not-empty div = ',avelis
**      print *,'Average n. of objects per div = ',aveave
**            print *,' '
***
 110     continue
        if(jflag.eq.0)then
         tmpq=10000
         tmpm=10000
         itmpq=0
         itmpm=0
         do 111 iaxis=1,7
          if(qualit(iaxis).lt.tmpq)then
           tmpq=qualit(iaxis)
           itmpq=iaxis
          endif
          if(tmpmax(iaxis).lt.tmpm)then
           tmpqm=tmpmax(iaxis)
           itmpm=iaxis
          endif
 111     continue
         if(iswit(9).eq.12345)then
          print *,'nin=',nin,' iax=',itmpq,' ndiv=',ndivto(itmpq)
          print *,'Max n. of objects per div = ',tmpmax(itmpq)
          print *,'Average n. of objects per div = ',tmpq
         endif
*** ONLY FOR DEBUG
**         if(isearc.lt.0)then
**           jsb=lq(lq(jvo-nin-1))
**           iaxor=q(jsb+1)
**           ndivor=q(jsb+2)-1
**           jsco=lq(jvo-nin-2)
**           tmpqor=0.
**           tmpmor=0.
**           do 133 idivor=1,ndivor
**             if(iq(jsco+idivor).gt.tmpmor)tmpmor=iq(jsco+idivor)
**             tmpqor=tmpqor+iq(jsco+idivor)
** 133       continue
**           tmpqor=tmpqor/ndivor
**           print *,'Gsord: iax=',iaxor,' ndiv=',ndivor
**           print *,'Gsord: Max n. of obj. per div = ',tmpmor
**           print *,'Gsord: Aver. n. of obj. per div = ',tmpqor
**         endif
***
         ndivst=ndivst+(ndivto(itmpq)+ndivto(itmpq)*(3.+tmpq)+10.)
         jflag=1
         goto 1
        else
          jflag=0
*** ONLY FOR DEBUG
**          print *,'nin=',nin,' iax=',q(jvirt+4*(ivo-1)+1),' ndiv=',
**     +q(jvirt+4*(ivo-1)+2)
**          ittmp=0
**          iind=q(jvirt+4*(ivo-1)+2)
**          do 155 n=1,iind
**           jvdiv=lq(jvirt-ivo)
**           iofset=iq(jvdiv+n)
**           nnobj=iq(jvdiv+iofset+1)
**           if(nnobj.gt.ittmp)ittmp=nnobj
** 155      continue
**          print *,'Max n. of objects per div = ',ittmp
**          print *,' '
**          print *,' '
***
        endif
        endif
        if(nin.gt.ninmax)ninmax=nin
 101  continue
       nwtota=ndivst+nvolum*5+10.
       if(iswit(9).eq.12345)then
        print *,'Computed number of words foreseen = ',nwtota
       endif
       nwreal=nwjvir+nwjvdi
       if(iswit(9).eq.12345)then
        print *,'Computed number of words booked = ',nwreal
       endif
       nwneed=nwjvir+ndioff
       if(iswit(9).eq.12345)then
        print *,'Computed number of words needed = ',nwneed
       endif
       if(jphi2.ne.0)call mzdrop(ixstor,jphi2,' ')
       if(jclow.ne.0)call mzdrop(ixstor,jclow,' ')
       if(jchig.ne.0)call mzdrop(ixstor,jchig,' ')
       if(jbuff.ne.0)call mzdrop(ixstor,jbuff,' ')
*
********************************************************************************
*
* *** Scan the volume structure to retrieve the path through
*      the physical tree for all sensitive detectors
*
       CALL GHCLOS
*
* *** Books STAT banks if data card STAT is submitted
*
      IF (NSTAT.GT.0)  CALL GBSTAT
*
      CALL MZGARB (IXCONS, 0)
*
 1001 FORMAT (' Severe diagnostic in initialization phase. STOP')
 1002 FORMAT (' GGCLOS : NVOLUM =',I5,' *****')
 1003 FORMAT (' Illegal tracking medium number in volume : ',A4)
 1004 FORMAT (' GGORDQ : Volume ',A4,' has more than 500 (',
     +        I3,') daughters ; volume sorting not possible !')
*                                                             END GGCLOS
  999 END
 
#endif
