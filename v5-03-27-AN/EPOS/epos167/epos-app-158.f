c 10.04.2003 Main program and random number generator of epos
c (Do not compile it for CONEX or Corsika)

c-----------------------------------------------------------------------
      program aamain
c-----------------------------------------------------------------------

      include 'epos.inc'

      save nopeno
      call aaset(0)
      call atitle
      call xiniall

9999  continue
      call aread
      if(nopen.eq.-1) then !after second round aread
        nopen=nopeno
        iecho=1
        call xiniall
        goto 9999
      endif
      call utpri('aamain',ish,ishini,4)
      if(model.ne.1)call IniModel(model)

      if(iappl.ne.0)then
        do nrebin= 1,noebin
          call ainit
                                ! if (nrebin.eq.1) call xini
          if(nevent.gt.0)then
            call aseed(2)
            write(ifmt,'(a,i10,a,f10.2,a)')'generate',nevent
     &       ,'  events for engy =',engy,'  ...'
            if(icotabr.eq.1)nevent=1
            if(icotabr.eq.1)nrevt=1
            if(nfull.gt.0)nevent=nfreeze*nfull
            if(mod(nevent,nfreeze).ne.0)
     &        stop'nevent must be a multiple of nfreeze!!!!!!!!!      '
            if(istore.ne.0) call bstora
            do  n=1,nevent
              call evgen(n)
              if(istore.ge.1.and.istore.le.4) call bstore
              if(istore.eq.5) call ustore
            enddo
            call astati
          else
            call xana
          endif
        enddo
        call bfinal
      endif
      if(istore.eq.3) write(ifdt,*) ' 0 0 '
 99   write(6,'(a)')'rewind copy-file'
      rewind (ifcp)
      nopeno=nopen
      nopen=-1
      iecho=0
      call utprix('aamain',ish,ishini,4)
      goto 9999

      end

      subroutine evgen(n)
      include 'epos.inc'

      if(irewch.eq.1)rewind(ifch)
      nfr=mod(n-1,nfreeze)*ispherio
      if(nfr.ne.0)goto77
      do nin=1,iabs(ninicon)
        if(icotabr.eq.0)call aepos(isign(1,-ninicon)*nin)
        if(icocore.ne.0)call IniCon(nin)
      enddo
  77  if(ispherio.ne.0)then
      if(nfr.eq.0)write(ifmt,'(a)')
     & 'spherio evolution + hadronization ...'
      if(mod(nfr+1,50).eq.0)
     & write(ifmt,*)'hadronization ',nfr+1,' / ',nfreeze
c                  call spherio2(nrevt,nfr)
      if(ish.ge.2)call alist('list after spherio&',1,nptl)
      call decayall(1)
      if(ish.ge.2)call alist('list after decay&',1,nptl)
      endif
      call aafinal
      if(iurqmd.eq.1)call urqmd
      call afinal
      call xana
      if(nfr+1.eq.nfreeze.or.ispherio.eq.0)call aseed(1)

      end

      subroutine setinp(ifname, nifname)
      include 'epos.inc'
      character*255 ifname

      ifop=19
      open(UNIT=ifop,FILE=ifname(1:nifname),STATUS='UNKNOWN')
      write(*, '(a,a)')'Connecting input to file: ', ifname(1:nifname)

      end
c-----------------------------------------------------------------------
      function rangen()
c-----------------------------------------------------------------------
c     generates a random number
c-----------------------------------------------------------------------
      include 'epos.inc'
1     rangen=ranf()
      if(rangen.le.0.)goto1
      if(rangen.ge.1.)goto1
      if(irandm.eq.1)write(ifch,*)'rangen()= ',rangen

      return
      end

c-----------------------------------------------------------------------
      double precision function drangen(dummy)
c-----------------------------------------------------------------------
c     generates a random number
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision dummy,dranf
1     drangen=dranf()
      if(drangen.le.0.d0)goto1
      if(drangen.ge.1.d0)goto1
      if(irandm.eq.1)write(ifch,*)'rangen()= ',drangen

      return
      end

c-----------------------------------------------------------------------
      function cxrangen(dummy)
c-----------------------------------------------------------------------
c     generates a random number
c-----------------------------------------------------------------------
      include 'epos.inc'
      double precision dummy
1     cxrangen=ranf()
      if(cxrangen.le.0.)goto1
      if(cxrangen.ge.1.)goto1
      if(irandm.eq.1)write(ifch,*)'rangen()= ',cxrangen

      return
      end

c-----------------------------------------------------------------------
      real function ranf()
c-----------------------------------------------------------------------
c     uniform random number generator from cern library
c-----------------------------------------------------------------------
      double precision    dranf,    g900gt,   g900st
      double precision    ds(2),    dm(2),    dseed
      double precision    dx24,     dx48
      double precision    dl,       dc,       du,       dr
      logical             single
      data      ds     /  1665 1885.d0, 286 8876.d0  /
      data      dm     /  1518 4245.d0, 265 1554.d0  /
      data      dx24   /  1677 7216.d0  /
      data      dx48   /  281 4749 7671 0656.d0  /
      single  =  .true.
      goto 10
      entry dranf()
      single  =  .false.
  10  dl  =  ds(1) * dm(1)
      dc  =  dint(dl/dx24)
      dl  =  dl - dc*dx24
      du  =  ds(1)*dm(2) + ds(2)*dm(1) + dc
      ds(2)  =  du - dint(du/dx24)*dx24
      ds(1)  =  dl
      dr     =  (ds(2)*dx24 + ds(1)) / dx48
      if(single)  then
         ranf  =  sngl(dr)
      else
         dranf  =  dr
      endif
      return
      entry g900gt()
      g900gt  =  ds(2)*dx24 + ds(1)
      return
      entry g900st(dseed)
      ds(2)  =  dint(dseed/dx24)
      ds(1)  =  dseed - ds(2)*dx24
      g900st =  ds(1)
      return
      end

c-----------------------------------------------------------------------
      subroutine ranfgt(seed)
c-----------------------------------------------------------------------
      double precision    seed,     g900gt,   g900st,   dummy
      seed  =  g900gt()
      return
      entry ranfst(seed)
      dummy  =  g900st(seed)
      return
      end
