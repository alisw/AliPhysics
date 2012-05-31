c  reshuffled from sem, sto, sha

c    contains DIS, and unused 3P stuff
c             ---      ---------------




c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################
c
c                              DIS
c
c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################










c-----------------------------------------------------------------------
      subroutine lepexp(rxbj,rqsq)
c-----------------------------------------------------------------------
c     generates x_bjorken and q**2 according to an experimental
c     distribution ( given in array xq(nxbj,nqsq) ).
c-----------------------------------------------------------------------
      parameter (nxbj=10,nqsq=10)
      parameter (xbjmin=0.,qsqmin=4.)
      parameter (xbjwid=0.025, qsqwid=4.)
      dimension xq(nxbj,nqsq),vxq(nxbj*nqsq)
      equivalence (xq(1,1),vxq(1))
 
      data (vxq(i),i=1,50)/
     &         1304.02,   366.40,    19.84,    10.79,     6.42,
     &            4.54,     4.15,     3.38,     2.03,     1.56,
     &          241.63,  1637.26,   427.36,   164.51,    73.72,
     &           43.07,    20.73,    12.78,     9.34,     5.83,
     &            0.01,   724.66,   563.79,   275.08,   176.13,
     &          106.44,    85.82,    54.52,    37.12,    28.65,
     &            0.01,   202.40,   491.10,   245.13,   157.07,
     &          104.43,    61.05,    49.42,    37.84,    26.79,
     &            0.01,     3.77,   316.38,   226.92,   133.45,
     &           90.30,    63.67,    48.42,    35.73,    28.04/
      data (vxq(i),i=51,100)/
     &            0.01,     0.01,   153.74,   213.09,   114.14,
     &           76.26,    60.02,    43.15,    43.47,    25.60,
     &            0.01,     0.01,    39.31,   185.74,   108.56,
     &           88.40,    47.29,    39.35,    31.80,    22.91,
     &            0.01,     0.01,     0.01,   104.61,   107.01,
     &           66.24,    45.34,    37.45,    33.44,    23.78,
     &            0.01,     0.01,     0.01,    56.58,    99.39,
     &           67.78,    43.28,    35.98,    34.63,    18.31,
     &            0.01,     0.01,     0.01,    13.56,    76.25,
     &           64.30,    42.80,    28.56,    21.19,    20.75 /
 
      data init/0/
      init=init+1
      if(init.eq.1) then
      n=nxbj*nqsq
      sum=0.
      do 1 i=1,n
      sum=sum+vxq(i)
1     continue
      do 2 i=2,n
2     vxq(i)=vxq(i)+vxq(i-1)
      do 3 i=1,n
3     vxq(i)=vxq(i)/sum
      endif
 
      n=nxbj*nqsq
      r=rangen()
      call utloc(vxq,n,r,iloc)
      if(iloc.ge.n) iloc=iloc-1
      i=mod(iloc,nxbj)+1
      if(i.eq.0) i=nxbj
      j=iloc/nxbj + 1
      dxint=vxq(1)
      if(iloc.gt.0) dxint=vxq(iloc+1)-vxq(iloc)
      dxbj=xbjwid*abs(r-vxq(iloc+1))/dxint
      dy=qsqwid*rangen()
      rxbj=xbjmin+xbjwid*float(i-1)+dxbj
      rqsq=qsqmin+qsqwid*float(j-1)+dy
      return
      end

c-----------------------------------------------------------------------
      subroutine fremny(wp1,wm1,pnx,pny,sm,ic1,ic2,ic3,ic4,coord,ey0)
c-----------------------------------------------------------------------
c  treats remnant from deep inelastic process;
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incsem'
      dimension coord(6),ic(2),ep(4),ey(3),ey0(3),ep3(4)
      double precision  ept(4),ept1(4)

      call utpri('fremny',ish,ishini,5)
      if(ish.ge.5)write(ifch,*)'writing remnant'
      if(ish.ge.5)write(ifch,*)
     *'wp1,wm1,pnx,pny,sm,ic1,ic2,ic3,ic4,coord,ey0:'
      if(ish.ge.5)write(ifch,*)
     *wp1,wm1,pnx,pny,sm,ic1,ic2,ic3,ic4,coord,ey0

        if(ic3.eq.0.and.ic4.eq.0)then

      ic(1)=ic1
      ic(2)=ic2
      nptl=nptl+1
      ep3(3)=pnx
      ep3(4)=pny
      ep3(2)=(wp1-wm1)/2
      ep3(1)=(wp1+wm1)/2
      call pstrans(ep3,ey0,-1)
      pptl(1,nptl)=ep3(3)
      pptl(2,nptl)=ep3(4)
      pptl(3,nptl)=ep3(2)
      pptl(4,nptl)=ep3(1)
      pptl(5,nptl)=sqrt(sm)
      idptl(nptl)=idtra(ic,0,0,3)
      iorptl(nptl)=1
      istptl(nptl)=0
      jorptl(nptl)=0
      do i=1,4
      xorptl(i,nptl)=coord(i)
      enddo
      tivptl(1,nptl)=coord(5)
      tivptl(2,nptl)=coord(6)
      ityptl(nptl)=40

       if(ish.ge.7)then
      write(ifch,*)'proj: nptl, mass**2,id',nptl,sm,idptl(nptl)
      write(ifch,*)'ept',
     *pptl(1,nptl),pptl(2,nptl),pptl(3,nptl),pptl(4,nptl)
       endif

        else

      ic(1)=ic1
      ic(2)=ic2
      nptl=nptl+1
      idptl(nptl)=idtra(ic,0,0,3)
      istptl(nptl)=20
      iorptl(nptl)=1
      jorptl(nptl)=0
      do i=1,4
      xorptl(i,nptl)=coord(i)
      enddo
      tivptl(1,nptl)=coord(5)
      tivptl(2,nptl)=coord(6)
      ityptl(nptl)=40

      ic(1)=ic3
      ic(2)=ic4
      idptl(nptl+1)=idtra(ic,0,0,3)
      istptl(nptl+1)=20
      iorptl(nptl+1)=1
      jorptl(nptl+1)=0
      do i=1,4
      xorptl(i,nptl+1)=coord(i)
      enddo
      tivptl(1,nptl+1)=coord(5)
      tivptl(2,nptl+1)=coord(6)
      ityptl(nptl+1)=40

      ep3(3)=pnx
      ep3(4)=pny
      ep3(2)=(wp1-wm1)/2
      ep3(1)=(wp1+wm1)/2
      call pstrans(ep3,ey0,-1)   !boost to hadronic c.m.s.
      ept(1)=ep3(3)
      ept(2)=ep3(4)
      ept(3)=ep3(2)
      ept(4)=ep3(1)
      do i=1,4
        ept1(i)=ep3(i)
      enddo

      sww=sqrt(sm)
      call psdeftr(sm,ept1,ey)
      ep(1)=.5*sww
      ep(2)=.5*sww
      ep(3)=0.
      ep(4)=0.
      call pstrans(ep,ey,1)
      pptl(1,nptl)=ep(3)
      pptl(2,nptl)=ep(4)
      pptl(3,nptl)=ep(2)
      pptl(4,nptl)=ep(1)
      do i=1,4
        pptl(i,nptl+1)=ept(i)-pptl(i,nptl)
      enddo
      nptl=nptl+1
        endif

 1000 continue
      if(ish.ge.5)write(ifch,*)'fremny: final nptl',nptl
      call utprix('fremny',ish,ishini,5)
      return
      end

c-----------------------------------------------------------------------
      subroutine psadis(iret)
c-----------------------------------------------------------------------
c psadis - DIS interaction
c-----------------------------------------------------------------------
      double precision ept(4),ept1(4),xx,wpt(2),eprt,pl,plprt,psutz
     *,psuds
      dimension ep3(4),ey(3),ey0(3),bx(6),
     *qmin(2),iqc(2),nqc(2),ncc(2,2),gdv(2),gds(2),dfp(4)
      parameter (mjstr=20000)
      common /psar29/ eqj(4,mjstr),iqj(mjstr),ncj(2,mjstr),nj
      common /psar30/ iorj(mjstr),ityj(mjstr),bxj(6,mjstr)
      double precision pgampr,rgampr
      common/cgampr/pgampr(5),rgampr(4)
      parameter (ntim=1000)
      common/cprt/nprtj,pprt(5,ntim),idprt(ntim),iorprt(ntim)
     &,idaprt(2,ntim)
      common/ciptl/iptl
      include 'epos.inc'
      include 'epos.incsem'

      call utpri('psadis',ish,ishini,3)
      if(ish.ge.3)write (ifch,*)'engy,elepti,iolept:'
      if(ish.ge.3)write (ifch,*)engy,elepti,iolept
      nptl=nptl+1
      idptl(nptl)=1220
      istptl(nptl)=1
      nptlh=nptl
      iptl=nptl
      s00=1.

      pptl(1,nptl)=0.
      pptl(2,nptl)=0.
      pptl(3,nptl)=-engy/2
      pptl(4,nptl)=engy/2
      pptl(5,nptl)=0

1     continue
      if(iolept.eq.1)then
        wtot=engy**2
        engypr=wtot/4./elepti
        gdv01=psdh(ydmax*wtot,qdmin,iclpro,0)
        gdv02=psdh(ydmax*wtot,qdmin,iclpro,1)
        gds01=psdsh(ydmax*wtot,qdmin,iclpro,dqsh,0)
        gds02=psdsh(ydmax*wtot,qdmin,iclpro,dqsh1,1)
        gb0=(1.+(1.-ydmax)**2)*(gdv01+gds01)
     *  +2.*(1.-ydmax)*(gdv02+gds02)

2       continue
        qq=qdmin*(qdmax/qdmin)**rangen()
        yd=ydmin*(ydmax/ydmin)**rangen()
        wd=yd*wtot
        if(ish.ge.4)write (ifch,*)'qq,wd,yd,ydmin,ydmax:'
        if(ish.ge.4)write (ifch,*)qq,wd,yd,ydmin,ydmax
        if(wd.lt.qq)goto 2

        gdv(1)=psdh(wd,qq,iclpro,0)
        gdv(2)=psdh(wd,qq,iclpro,1)
        gds(1)=psdsh(wd,qq,iclpro,dqsh,0)
        gds(2)=psdsh(wd,qq,iclpro,dqsh1,1)
        gbtr=(1.+(1.-yd)**2)*(gdv(1)+gds(1))
        gblong=2.*(1.-yd)*(gdv(2)+gds(2))
c        gblong=0.   !???????
        gb=(gbtr+gblong)/gb0*.7
        if(ish.ge.4)then
          if(gb.gt.1.)write(ifmt,*)'gb,qq,yd,wd',gb,qq,yd,wd
          write (ifch,*)'gb,gdv,gds,gdv0,gds0,yd:'
          write (ifch,*)gb,gdv,gds,gdv01,gds01,
     *    gdv02,gds02,yd
        endif
        if(rangen().gt.gb)goto 2
        
        long=int(rangen()+gblong/(gbtr+gblong))
        elepto=qq/elepti/4.+elepti*(1.-yd)
        costhet=1.-qq/elepti/elepto/2.
        theta=acos(costhet)
        if(theta/pi*180..lt.themin)goto 2
        if(theta/pi*180..gt.themax)goto 2
        if(elepto.lt.elomin)goto 2
        if(ish.ge.3)write (ifch,*)'theta,elepto,elepti,iclpro:'
        if(ish.ge.3)write (ifch,*)theta/pi*180.,elepto,elepti,iclpro
        xbjevt=qq/wd
        qsqevt=qq

        call pscs(bcos,bsin)
        rgampr(1)=-elepto*sin(theta)*bcos
        rgampr(2)=-elepto*sin(theta)*bsin
        rgampr(3)=elepti-elepto*costhet
        rgampr(4)=elepti-elepto

        pgampr(1)=rgampr(1)
        pgampr(2)=rgampr(2)
        pgampr(3)=rgampr(3)-engypr
        pgampr(4)=rgampr(4)+engypr
        sm2=pgampr(4)*pgampr(4)
     *  -pgampr(1)*pgampr(1)-pgampr(2)*pgampr(2)-pgampr(3)*pgampr(3)
        pgampr(5)=sqrt(sm2)
        call utlob2(1,pgampr(1),pgampr(2),pgampr(3),pgampr(4),pgampr(5)
     *  ,rgampr(1),rgampr(2),rgampr(3),rgampr(4),40)
        if(ish.ge.4)write (ifch,*)'rgampr:',rgampr

      elseif(iolept.lt.0)then
21      call lepexp(xbjevt,qsq)
        qq=qsq
        wd=qq/xbjevt
        if(qq.lt.qdmin.or.qq.gt.qdmax)goto21
        
        gdv(1)=psdh(wd,qq,iclpro,0)
        gdv(2)=psdh(wd,qq,iclpro,1)
        gds(1)=psdsh(wd,qq,iclpro,dqsh,0)
        gds(2)=psdsh(wd,qq,iclpro,dqsh1,1)
        yd=wd/engy**2
        gbtr=(1.+(1.-yd)**2)*(gdv(1)+gds(1))
        gblong=2.*(1.-yd)*(gdv(2)+gds(2))
        gblong=0. !????????????
        long=int(rangen()+gblong/(gbtr+gblong))
      else
        stop'wrong iolept'
      endif
      if(ish.ge.3)write (ifch,*)'qq,xbj,wd,gdv,gds,dqsh:'
      if(ish.ge.3)write (ifch,*)qq,xbjevt,wd,gdv,gds,dqsh

      egyevt=sqrt(wd-qq)
      pmxevt=.5*egyevt

      wp0=sqrt(qq)       !breit frame
      wm0=(wd-qq)/wp0
      ey0(1)=egyevt/wp0  !boost to the hadronic c.m.s.
      ey0(2)=1.
      ey0(3)=1.
      do i=1,6
        bx(i)=0.
      enddo

      if(long.eq.0)then
        sdmin=qq/(1.-sqrt(q2ini/qq))
        sqmin=sdmin
      else
        sdmin=4.*max(q2min,qcmass**2)+qq  !minimal mass for born
        xmm=(5.*sdmin-qq)/4.
        sqmin=1.1*(xmm+sqrt(xmm**2-qq*(sdmin-qq-4.*q2ini)))
     *  /2./(1.-4.*q2ini/(sdmin-qq))
      endif
      if(long.eq.1.and.wd.lt.1.001*sdmin)goto 1
      
      proja=210000.
      projb=0.
      call fremnu(ammin,proja,projb,proja,projb,
     *icp1,icp2,icp3,icp4)

      nj=0

      if((rangen().lt.gdv(long+1)/(gdv(long+1)+gds(long+1)).or.
     *egyevt.lt.1.0001*(ammin+sqrt(sdmin-qq))).and.
     *(long.eq.0.or.wd.gt.sqmin))then 
        if(long.eq.0)then
          xd=qq/wd
          tu=psdfh4(xd,q2min,0.,iclpro,1)/2.25
          td=psdfh4(xd,q2min,0.,iclpro,2)/9.
          gdv0=(tu+td)*4.*pi**2*alfe/qq
     *    *sngl(psuds(qq,1)/psuds(q2min,1))
          if(ish.ge.4)write (ifch,*)'gdv0:',gdv0,sdmin
          
          if(rangen().lt.gdv0/gdv(1).or.wd.le.1.0001*sdmin)then    !?????
            if(ish.ge.3)write (ifch,*)'no cascade,gdv0,gdv',gdv0,gdv
            if(rangen().lt.tu/(tu+td))then
              iq1=1
              izh=3
            else
              iq1=2
              izh=6
            endif
            jq=1
            if(ish.ge.8)write (ifch,*)'before call timsh2: ',
     *      'qq,egyevt,iq1',qq,egyevt,iq1
            call timsh2(qq,0.,egyevt,iq1,-iq1)

            nj=nj+1
            iqj(nj)=izh
            nqc(1)=nj
            nqc(2)=0

            ep3(1)=pprt(4,2)
            ep3(2)=pprt(3,2)
            ep3(3)=0.
            ep3(4)=0.
            call pstrans(ep3,ey0,1)
            do i=1,4
              eqj(i,nj)=ep3(i)
            enddo
            s0h=0.
            c0h=1.
            s0xh=0.
            c0xh=1.
            call psreti(nqc,jq,1,ey0,s0xh,c0xh,s0h,c0h)
            goto 17
          endif
        endif

        call psdint(wd,qq,sds0,sdn0,sdb0,sdt0,sdr0,1,long)
        if(ish.ge.3)write (ifch,*)'wd,qq,sds0,sdn0,sdb0,sdt0,sdr0:'
        if(ish.ge.3)write (ifch,*)wd,qq,sds0,sdn0,sdb0,sdt0
        gb10=(sdn0+sdt0)*(1.-qq/wd)
        xdmin=sqmin/wd

3       continue
        xd=(xdmin-qq/wd)/((xdmin-qq/wd)/(1.-qq/wd))
     *  **rangen()+qq/wd
        call psdint(xd*wd,qq,sds,sdn,sdb,sdt,sdr,1,long)
        if(ish.ge.3)write (ifch,*)'wdhard,qq,sds,sdn,sdb,sdt:'
        if(ish.ge.3)write (ifch,*)xd*wd,qq,sds,sdn,sdb,sdt
        tu=psdfh4(xd,q2min,0.,iclpro,1)
        td=psdfh4(xd,q2min,0.,iclpro,2)
        gb1=(sdn*(tu/2.25+td/9.)+sdt*(tu+td)/4.5)
     *  *(1.-qq/wd/xd)/gb10
        if(gb1.gt.1..and.ish.ge.1)write(ifmt,*)'gb1,xd,wd,qq,sdt0,sdt',
     *   gb1,xd,wd,qq,sdt0,sdt
        if(ish.ge.6)write (ifch,*)'gb1,xd,wd,qq,sdt0,sdt:'
        if(ish.ge.6)write (ifch,*)gb1,xd,wd,qq,sdt0,sdt
        if(rangen().gt.gb1)goto 3

        gdres=(sdt-sds)/4.5
        gdrga=sdr/4.5
        gdsin=sds/4.5
        dtu=tu*(sdn/2.25+sdt/4.5)
        dtd=td*(sdn/9.+sdt/4.5)
        if(rangen().lt.dtu/(dtu+dtd))then
          iq1=1
          izh=3
          gdbor=sdb/2.25
          gdnon=sdn/2.25
        else
          iq1=2
          izh=6
          gdbor=sdb/9.
          gdnon=sdn/9.
        endif

        wpi=wp0
        wmi=(xd*wd-qq)/wpi
        iqc(2)=iq1
        nj=nj+1
        iqj(nj)=izh
        eqj(1,nj)=.5*(wm0-wmi)
        eqj(2,nj)=-eqj(1,nj)
        eqj(3,nj)=0.
        eqj(4,nj)=0.
        ncc(1,2)=nj
        ncc(2,2)=0
        if(ish.ge.3)write (ifch,*)'wp0,wm0,wpi,wmi,iqc(2),eqj'
        if(ish.ge.3)write (ifch,*)wp0,wm0,wpi,wmi,iqc(2),eqj(2,nj)

      else
        xdmin=sdmin/wd
        xpmax=((egyevt-ammin)**2+qq)/wd
        iq1=int(3.*rangen()+1.)*(2.*int(.5+rangen())-1.)

        aks=rangen()
        if(long.eq.0.and.aks.lt.dqsh/gds(1).and.
     *  egyevt.gt.ammin+sqrt(s00))then
          if(ish.ge.3)write (ifch,*)'no cascade for q_s',
     *    aks,dqsh/gds(1)
          xd=qq/wd
          xpmin=xd+s00/wd
          jcasc=0
          if(iq1.gt.0)then
            jq=1
          else
            jq=2
          endif
        else
          jcasc=1
          call psdint(xpmax*wd,qq,sds0,sdn0,sdb0,sdt0,sdr0,0,long)
          call psdint(xpmax*wd,qq,sdsq0,sdnq0,sdbq0,sdtq0,sdrq0,1,long)
          if(ish.ge.3)write (ifch,*)
     *    'xpmax*wd,qq,sds0,sdn0,sdb0,sdt0,sdr0:'
          if(ish.ge.3)write (ifch,*)
     *    xpmax*wd,qq,sds0,sdn0,sdb0,sdt0,sdr0
        gb10=sdt0*fzeroGluZZ(0.,iclpro)+(sdnq0+sdtq0)
     *        *fzeroSeaZZ(0.,iclpro)
          gb10=gb10*15.

4         xd=xdmin*(xpmax/xdmin)**rangen()
          xpmin=xd
          call psdint(xd*wd,qq,sds,sdn,sdb,sdt,sdr,0,long)
          call psdint(xd*wd,qq,sdsq,sdnq,sdbq,sdtq,sdrq,1,long)
          if(ish.ge.3)write (ifch,*)'xd*wd,qq,sds,sdn,sdb,sdt,sdr:'
          if(ish.ge.3)write (ifch,*)xd*wd,qq,sds,sdn,sdb,sdt,sdr
          wwg=sdt*fzeroGluZZ(xd,iclpro)
          wwq=(sdnq+sdtq)*fzeroSeaZZ(xd,iclpro)
          gb12=(wwq+wwg)/gb10*(xpmax/xd)**dels
          if(gb12.gt.1..and.ish.ge.1)write(ifmt,*)
     *    'gb12,xpmax*wd,xd*wd,sdt0,sdnq0+sdtq0,sdt,sdnq+sdtq',
     *    gb12,xpmax*wd,xd*wd,sdt0,sdnq0+sdtq0,sdt,sdnq+sdtq,
     *    wwq,wwg,(xpmax/xd)**dels,gb10
          if(ish.ge.5)write (ifch,*)'gb12,xd,xpmax,wwq,wwg:'
          if(ish.ge.5)write (ifch,*)gb12,xd,xpmax,wwq,wwg
          if(rangen().gt.gb12)goto 4
        endif

        if(jcasc.ne.0)then
          gb20=(1.-xd/xpmax)**betpom*sdt*(1.-glusea)+
     *    EsoftQZero(xd/xpmax)*(sdnq+sdtq)*glusea
        else
          gb20=EsoftQZero(xd/xpmax)
        endif
        if(1.+2.*(-alpqua)+dels.ge.0.)then
          xpminl=(1.-xpmax)**(alplea(iclpro)+1.)
          xpmaxl=(1.-xpmin)**(alplea(iclpro)+1.)

5         xp=1.-(xpminl+(xpmaxl-xpminl)*rangen())**
     *    (1./(alplea(iclpro)+1.))
          if(jcasc.ne.0)then
            gb2=((1.-xd/xp)**betpom*sdt*(1.-glusea)+
     *      EsoftQZero(xd/xp)*(sdnq+sdtq)*glusea)*(xp/xpmax)**
     *      (1.+2.*(-alpqua)+dels)/gb20
          else
           gb2=EsoftQZero(xd/xp)*(xp/xpmax)**(1.+2.*(-alpqua)+dels)/gb20
          endif
          if(gb2.gt.1..and.ish.ge.1)then
            write(ifmt,*)'gb2,xp:',gb2,xp
c            read (*,*)
          endif
          if(rangen().gt.gb2)goto 5
        else
          xpmaxl=xpmax**(2.+2.*(-alpqua)+dels)
          xpminl=xpmin**(2.+2.*(-alpqua)+dels)

6         xp=(xpminl+(xpmaxl-xpminl)*rangen())**
     *    (1./(2.+2.*(-alpqua)+dels))
          if(jcasc.ne.0)then
            gb21=((1.-xd/xp)**betpom*sdt*(1.-glusea)+
     *      EsoftQZero(xd/xp)*(sdnq+sdtq)*glusea)*
     *      ((1.-xp)/(1.-xd))**alplea(iclpro)/gb20
          else
          gb21=EsoftQZero(xd/xp)*((1.-xp)/(1.-xd))**alplea(iclpro)/gb20
          endif
          if(gb21.gt.1..and.ish.ge.1)then
            write(ifmt,*)'gb21,xp:',gb21,xp
c            read (*,*)
          endif
          if(rangen().gt.gb21)goto 6
        endif

        wwh=xd*wd-qq
        wwsh=xp*wd-qq
        ammax=(egyevt-sqrt(wwsh))**2
22      call fremnx(ammax,ammin,sm,icp3,icp4,iret)
        if(iret.ne.0.and.ish.ge.1)write(ifmt,*)'iret.ne.0!'
     *                                         ,ammax,ammin**2
        wmn=(1.-xp)*wd/wp0
        wpn=sm/wmn
        pnx=0.
        pny=0.
        wpp=wp0-wpn
        wmp=wm0-wmn
        if(ish.ge.5)write(ifch,*)'wp0,wm0,wpn,wmn,wpp,wmp:'
        if(ish.ge.5)write(ifch,*)wp0,wm0,wpn,wmn,wpp,wmp

        if(jcasc.eq.0.or.rangen().lt.wwq/(wwg+wwq).
     *  and.xd*wd.gt.sqmin.and.wwsh.gt.
     *  (sqrt(wwh)+sqrt(s00))**2)then
          zgmin=xd/xp
          zgmax=1./(1.+wp0/xd/wd/(wpp-wwh/wmp))
          if(zgmin.gt.zgmax)goto 22
23        zg=zgmin-rangen()*(zgmin-zgmax)
          if(rangen().gt.zg**dels*((1.-xd/xp/zg)/ (1.-xd/xp))**betpom)
     *    goto 23
          xg=xd/zg             !w- share for the struck quark
          wmq=wd/wp0*(xg-xd)   !w- for its counterpart
          wpq=s00/wmq          !1. gev^2 / wmq
          wmq=0.
          wpp=wpp-wpq
          wmp=wmp-wmq
          sxx=wpp*wmp
          if(ish.ge.5)write (ifch,*)'wpq,wmq,wpp,wmp,sxx:'
          if(ish.ge.5)write (ifch,*)wpq,wmq,wpp,wmp,sxx

          if(jcasc.eq.0)then
            if(ish.ge.6)write (ifch,*)'before call timsh2: qq,sxx,iq1',
     *      qq,sxx,iq1
            call timsh2(qq,0.,sqrt(sxx),iq1,-iq1)
            ept(1)=.5*(wpp+wmp)
            ept(2)=.5*(wpp-wmp)
            ept(3)=0.
            ept(4)=0.
            call psdeftr(sxx,ept,ey)
            ep3(1)=pprt(4,2)
            ep3(2)=pprt(3,2)
            ep3(3)=0.
            ep3(4)=0.

            call pstrans(ep3,ey,1)
            wmp=ep3(1)-ep3(2)
            goto 24
          endif
        else
          iq1=0
          sxx=wpp*wmp
        endif

        if(ish.ge.3)write (ifch,*)'wwh,wwsh,sxx,wpp,wmp:',
     *  wwh,wwsh,sxx,wpp,wmp

        wpi=wpp
        wmi=wwh/wpp
        wmp=wmp-wmi
24      call fremny(wpn,wmn,pnx,pny,sm,icp1,icp2,icp3,icp4,bx,ey0)

        if((-alpqua).eq.-1.)stop'dis does not work for 1/x'
25      aks=rangen()
        z=.5*aks**(1./(1.+(-alpqua)))
        if(z.lt.1.e-5.or.rangen().gt.(2.*(1.-z))**(-alpqua))goto 25
        if(rangen().gt..5)z=1.-z
        wm2=wmp*z
        wm1=wmp-wm2

        iqc(2)=iq1
        nj=nj+1
        iqj(nj)=-int(2.*rangen()+1.)
        iqj(nj+1)=-iqj(nj)
        eqj(1,nj)=.5*wm1
        eqj(2,nj)=-eqj(1,nj)
        eqj(3,nj)=0.
        eqj(4,nj)=0.
        eqj(1,nj+1)=.5*wm2
        eqj(2,nj+1)=-eqj(1,nj+1)
        eqj(3,nj+1)=0.
        eqj(4,nj+1)=0.
        nj=nj+1

        if(iq1.eq.0)then
          ncc(1,2)=nj-1
          ncc(2,2)=nj
          gdres=sdt-sds
          gdrga=sdr
          gdsin=sds
          gdbor=sdb
          gdnon=sdn
        else
          nj=nj+1
          if(iabs(iq1).eq.3)then
            iqj(nj)=-iq1*4/3
          else
            iqj(nj)=-iq1
          endif
          eqj(1,nj)=.5*(wpq+wmq)
          eqj(2,nj)=.5*(wpq-wmq)
          eqj(3,nj)=0.
          eqj(4,nj)=0.
          if(iq1.gt.0)then
            ncj(1,nj)=nj-1
            ncj(1,nj-1)=nj
            ncj(2,nj)=0
            ncj(2,nj-1)=0
          else
            ncj(1,nj)=nj-2
            ncj(1,nj-2)=nj
            ncj(2,nj)=0
            ncj(2,nj-2)=0
          endif

          if(jcasc.eq.0)then
            if(iq1.gt.0)then
              nqc(1)=nj-2
              nqc(2)=0
            else
              nqc(1)=nj-1
              nqc(2)=0
            endif
            s0h=0.
            c0h=1.
            s0xh=0.
            c0xh=1.
            call psreti(nqc,jq,1,ey,s0xh,c0xh,s0h,c0h)
            goto 17
          else
            gdres=(sdtq-sdsq)/4.5
            gdrga=sdrq/4.5
            gdsin=sdsq/4.5
            gdbor=sdbq/4.5
            gdnon=sdnq/4.5
            if(iq1.gt.0)then
              ncc(1,2)=nj-2
              ncc(2,2)=0
            else
              ncc(1,2)=nj-1
              ncc(2,2)=0
            endif
          endif
        endif

        if(ish.ge.3)write (ifch,*)'wpn,wmn,wpi,wmi,wm1,wm2,nj'
        if(ish.ge.3)write (ifch,*)wpn,wmn,wpi,wmi,wm1,wm2,nj
      endif
      
      si=wpi*wmi+qq
      qmin(2)=q2min                 !effective momentum cutoff below 
      s2min=max(4.*qq,16.*q2min)    !mass cutoff for born scattering

      if(rangen().gt.gdres/(gdres+gdsin+gdnon).or.
     *si.lt.(s2min+qq))goto 12

c---------------------------------------
c hard pomeron (resolved photon)
c---------------------------------------
      if(ish.ge.3)write(ifmt,*)'resolved,gdrga,gdres',gdrga,gdres

      jj=1
      if(rangen().gt.gdrga/gdres.and.si.gt.1.1*s2min+qq)then 
        if(ish.ge.3)write(ifmt,*)'dir-res,si,qq',si,qq     
        pt=0.
        pt2=0.
        iqc(1)=0
        ept(1)=.5*(wpi+wmi)
        ept(2)=.5*(wpi-wmi)
        ept(3)=0.
        ept(4)=0.
        wpt(1)=wpi           !lc+ for the current jet emission
        wpt(2)=si/wpi        !lc- for the current jet emission
      
        qqmin=max(q2min,s2min/(si/qq-1.))
        qqmax=min(si/2.,si-s2min)
        qmin(1)=qqmin
        xmax=1.
        xmin=(s2min+qq)/si
        if(qqmin.ge.qqmax.or.xmin.ge.xmax)stop'min>max'
        gb0=psjti(qmin(1),qq,si-qq,7,iqc(2),1)*psfap(1.d0,0,1)
     
        ncc(1,1)=0
        ncc(2,1)=0
        jgamma=1
        ntry=0
        goto 9
      else
        if(ish.ge.3)write(ifmt,*)'res,si,qq',si,qq     
        qmin(1)=q2min                 !effective momentum cutoff above
        si=si-qq
        zmin=s2min/si
        dft0=psdfh4(zmin,q2min,qq,0,0)*psjti(q2min,qq,si,0,iqc(2),1)
     *  +(psdfh4(zmin,q2min,qq,0,1)+psdfh4(zmin,q2min,qq,0,2)+
     *  psdfh4(zmin,q2min,qq,0,3))*psjti(q2min,qq,si,7,iqc(2),1)
     
7       continue
        z=zmin**rangen()
        do i=1,4
          dfp(i)=psdfh4(z,q2min,qq,0,i-1)
        enddo
        dfp(1)=dfp(1)*psjti(q2min,qq,z*si,0,iqc(2),1)
        dfptot=dfp(1)
        if(iqc(2).eq.0)then
          sjq=psjti(q2min,qq,z*si,1,0,1)
          do i=2,4
            dfp(i)=dfp(i)*sjq
            dfptot=dfptot+dfp(i)
          enddo
        else
          sjqqp=psjti(q2min,qq,z*si,1,2,1)
          do i=2,4
            if(iabs(iqc(2)).eq.i-1)then
              dfp(i)=dfp(i)*(psjti(q2min,qq,z*si,1,1,1)+
     *        psjti(q2min,qq,z*si,1,-1,1))/2.
            else
              dfp(i)=dfp(i)*sjqqp
            endif
            dfptot=dfptot+dfp(i)
          enddo
        endif
        
        if(rangen().gt.dfptot/dft0)goto 7
        
        wpq=wpi*(1.-z)
        wpi=wpi*z
        aks=dfptot*rangen()
        if(aks.lt.dfp(1))then
          iqc(1)=0
          nj=nj+1
          ncc(1,1)=nj
          ncc(2,1)=nj+1

          iqj(nj)=-int(2.*rangen()+1.)
          iqj(nj+1)=-iqj(nj)          
          wpq1=wpq*rangen()
          eqj(1,nj)=.5*wpq1
          eqj(2,nj)=eqj(1,nj)
          eqj(3,nj)=0.
          eqj(4,nj)=0.
          eqj(1,nj+1)=.5*(wpq-wpq1)
          eqj(2,nj+1)=eqj(1,nj+1)
          eqj(3,nj+1)=0.
          eqj(4,nj+1)=0.
          nj=nj+1

        else
          if(aks.lt.dfp(1)+dfp(2))then
            iqc(1)=1
          elseif(aks.lt.dfp(1)+dfp(2)+dfp(3))then
            iqc(1)=2
          else
            iqc(1)=3
          endif
          iqc(1)=iqc(1)*(2*int(2.*rangen())-1)
          nj=nj+1
          ncc(1,1)=nj
          ncc(2,1)=0
          
          iqj(nj)=-iqc(1)
          eqj(1,nj)=.5*wpq
          eqj(2,nj)=eqj(1,nj)
          eqj(3,nj)=0.
          eqj(4,nj)=0.
        endif
                 
        ept(1)=.5*(wpi+wmi)
        ept(2)=.5*(wpi-wmi)
        ept(3)=0.
        ept(4)=0.
        jgamma=0
        ntry=0
      endif

8     continue

c ladder rung
c---------------------------------------
      pt2=ept(3)**2+ept(4)**2
      pt=sqrt(pt2)

      wpt(1)=ept(1)+ept(2)              !lc+ for the current jet emissi
      wpt(2)=ept(1)-ept(2)              !lc- for the current jet emissi

      s2min=max(qmin(1),16.*qmin(2))    !mass cutoff for born 
      s2min=max(s2min,4.*qq)   
      
      if(jj.eq.1)then
        wwmin=2.*s2min-2.*pt*sqrt(q2ini)
        wwmin=(wwmin+sqrt(wwmin**2+4.*pt2*(s2min-q2ini)))
     *  /(1.-q2ini/s2min)/2.  
        sj=psjti(qmin(1),qq,si,iqc(1),iqc(2),1)   !total jet 
        sj2=psjti1(q2min,qmin(1),qq,si,iqc(2),iqc(1),1)
        if(ish.ge.3)write(ifch,*)'resolved - si,wwmin,s2min,sj,sj2:'
        if(ish.ge.3)write(ifch,*)si,wwmin,s2min,sj,sj2
        if(sj.eq.0.)stop'sj=0'
        if(rangen().gt.sj2/sj.and.si.gt.1.1*wwmin)goto 26
        jj=2
      endif
      sj=psjti1(qmin(2),qmin(1),qq,si,iqc(2),iqc(1),1)
      sjb=psbint(qmin(1),qmin(2),qq,si,iqc(1),iqc(2),1) !born parton-parton
      wwmin=17./16*s2min-2.*pt*sqrt(q2ini)
      wwmin=(wwmin+sqrt(wwmin**2+pt2*(s2min/4.-4.*q2ini)))
     */(1.-16.*q2ini/s2min)/2.  
      if(rangen().lt.sjb/sj.or.si.lt.1.1*wwmin)goto 10

26    continue
      wpt(jj)=wpt(jj)-pt2/wpt(3-jj)   
                    
      if(jj.eq.1)then
        discr=(si+2.*pt*sqrt(q2ini))**2-4.*q2ini*(2.*si+pt2)
        if(discr.lt.0..and.ish.ge.1)write(ifmt,*)'discr,si,pt,wwmin',
     *  discr,si,pt,wwmin
        discr=sqrt(discr)
        qqmax=(si+2.*pt*sqrt(q2ini)+discr)/2./(2.+pt2/si)
      else
        discr=(si+2.*pt*sqrt(q2ini))**2-4.*q2ini*(17.*si+pt2)
        if(discr.lt.0..and.ish.ge.1)write(ifmt,*)'discr,si,pt,wwmin',
     *  discr,si,pt,wwmin
        discr=sqrt(discr)
        qqmax=(si+2.*pt*sqrt(q2ini)+discr)/2./(17.+pt2/si)
      endif
      qqmin=2.*q2ini*si/(si+2.*pt*sqrt(q2ini)+discr)
      if(jj.eq.1.and.s2min.gt.qqmin.or.
     *jj.eq.2.and.s2min.gt.16.*qqmin)then
        xmm=.5*(si-s2min+2.*pt*sqrt(q2ini))
        discr=xmm**2-q2ini*(si+pt2)
        if(discr.lt.0..and.ish.ge.1)write(ifmt,*)'discr1,si,pt,wwmin',
     *  discr,si,pt,wwmin
        qqmin=q2ini*si/(xmm+sqrt(discr))
      endif

      xmin=1.-q2ini/qqmin
      xmax=1.-q2ini/qqmax
      if(ish.ge.6)write(ifch,*)'qqmin,qqmax,xmin,xmax',
     *qqmin,qqmax,xmin,xmax
      if(qqmin.lt.qmin(jj))then
        qqmin=qmin(jj)
        xmi=max(1.-((pt*sqrt(qqmin)+sqrt(pt2*qqmin+
     *  si*(si-s2min-qqmin*(1.+pt2/si))))/si)**2,
     *  (s2min+qqmin*(1.+pt2/si)-2.*pt*sqrt(qqmin))/si)
        xmin=max(xmin,xmi)
        if(xmin.le.0.)xmin=(s2min+qqmin*(1.+pt2/si))/si
        if(ish.ge.6)write(ifch,*)'qqmin,qmin(jj),xmin,s2min',
     *  qqmin,qmin(jj),xmin,s2min
      endif

      qm0=qmin(jj)
      xm0=1.-q2ini/qm0
      if(xm0.gt.xmax.or.xm0.lt.xmin)then
        xm0=.5*(xmax+xmin)
      endif
c      s2max=xm0*si
      s2max=xm0*si-qm0*(1.+pt2/si)+2.*pt*sqrt(q2ini)  !new ladder mass squared
      xx=xm0

      if(jj.eq.1)then
        sj0=psjti(qm0,qq,s2max,0,iqc(2),1)*psfap(xx,iqc(1),0)+
     *  psjti(qm0,qq,s2max,7,iqc(2),1)*psfap(xx,iqc(1),1)
        gb0=sj0/log(q2ini/qcdlam)*sngl(psuds(qm0,iqc(1)))*qm0*2.
      else
        sj0=psjti1(qm0,qmin(1),qq,s2max,0,iqc(1),1)*psfap(xx,iqc(2),0)
     *  +psjti1(qm0,qmin(1),qq,s2max,7,iqc(1),1)*psfap(xx,iqc(2),1)
        gb0=sj0/log(q2ini/qcdlam)*sngl(psuds(qm0,iqc(2)))*qm0*2.
      endif
      if(gb0.le.0.)then
         write(ifmt,*)'gb0.le.0.  si,qq,pt2:',si,qq,pt2
         iret=1
         goto 9999
      endif
      if(xm0.le..5)then
        gb0=gb0*xm0**(1.-delh)
      else
        gb0=gb0*(1.-xm0)*2.**delh
      endif

      xmin2=max(.5,xmin)
      xmin1=xmin**delh                 !xmin, xmax are put into powe
      xmax1=min(xmax,.5)**delh       !to simulate x value below
      if(xmin.ge..5)then
        djl=1.
      elseif(xmax.lt..5)then
        djl=0.
      else
        djl=1./(1.+((2.*xmin)**delh-1.)/delh/
     *  log(2.*(1.-xmax)))
      endif
      
      ntry=0
9     continue
      ntry=ntry+1
      if(ntry.ge.10000)then
        print *,"ntry.ge.10000"
        iret=1
        goto 9999
      endif
      if(jgamma.ne.1)then
        if(rangen().gt.djl)then        !lc momentum share in the cur
          x=(xmin1+rangen()*(xmax1-xmin1))**(1./delh)
        else
          x=1.-(1.-xmin2)*((1.-xmax)/(1.-xmin2))**rangen()
        endif
        q2=qqmin/(1.+rangen()*(qqmin/qqmax-1.))
        qt2=q2*(1.-x)
        if(ish.ge.6)write(ifch,*)'jj,q2,x,qt2',jj,q2,x,qt2
        if(qt2.lt.q2ini)goto 9
      else
        x=xmin+rangen()*(xmax-xmin)
        q2=qqmin*(qqmax/qqmin)**rangen()
        qt2=(q2-x*qq)*(1.-x)
        if(ish.ge.6)write(ifch,*)'jj,q2,x,qt2',jj,q2,x,qt2
        if(qt2.lt.0.)goto 9
      endif

      qt=sqrt(qt2)
      call pscs(bcos,bsin)
c ep3 is now 4-vector for s-channel gluon produced in current ladder run
      ep3(3)=qt*bcos
      ep3(4)=qt*bsin
      ptnew=(ept(3)-ep3(3))**2+(ept(4)-ep3(4))**2
      if(jj.eq.1)then
        s2min2=max(q2,s2min)
      else
        s2min2=max(s2min,16.*q2)
      endif

      if(jgamma.ne.1)then
        s2=x*si-q2*(1.+pt2/si)-ptnew+pt2+qt2  !new ladder mass squared
        if(s2.lt.s2min2)goto 9      !rejection in case of too low mass
        xx=x
        
        if(jj.eq.1)then
          sj1=psjti(q2,qq,s2,0,iqc(2),1)
          if(iqc(1).ne.0)then
            sj2=psjti(q2,qq,s2,iqc(1),iqc(2),1)
          elseif(iqc(2).eq.0)then
            sj2=psjti(q2,qq,s2,1,0,1)
          else
            sj2=psjti(q2,qq,s2,1,1,1)/6.+
     *      psjti(q2,qq,s2,-1,1,1)/6.+
     *      psjti(q2,qq,s2,2,1,1)/1.5
          endif
        else
          sj1=psjti1(q2,qmin(1),qq,s2,0,iqc(1),1)
          if(iqc(2).ne.0)then
            sj2=psjti1(q2,qmin(1),qq,s2,iqc(2),iqc(1),1)
          elseif(iqc(1).eq.0)then
            sj2=psjti1(q2,qmin(1),qq,s2,1,0,1)
          else
            sj2=psjti1(q2,qmin(1),qq,s2,1,1,1)/6.+
     *      psjti1(q2,qmin(1),qq,s2,-1,1,1)/6.+
     *      psjti1(q2,qmin(1),qq,s2,2,1,1)/1.5
          endif
        endif
c gb7 is the rejection function for x and q**2 simulation
        gb7=(sj1*psfap(xx,iqc(jj),0)+sj2*psfap(xx,iqc(jj),1))
     *  /log(qt2/qcdlam)*sngl(psuds(q2,iqc(jj)))*q2/gb0

        if(x.le..5)then
          gb7=gb7*x**(1.-delh)
        else
          gb7=gb7*(1.-x)*2.**delh
        endif
      else
        s2=x*si-q2               !new ladder mass squared
        if(s2.lt.s2min2)goto 9   !rejection in case of too low mass

        sj1=0.
        xx=x
        if(iqc(2).eq.0)then
          sj2=psjti(q2,qq,s2,1,0,1)
        else
          sj2=psjti(q2,qq,s2,1,1,1)/naflav/2.+
     *    psjti(q2,qq,s2,-1,1,1)/naflav/2.+
     *    psjti(q2,qq,s2,2,1,1)*(1.-1./naflav)
        endif
        gb7=sj2*psfap(xx,0,1)/gb0  !????*(1.-x*qq/q2) 
      endif
      if(gb7.gt.1..or.gb7.lt.0..and.ish.ge.1)write(ifmt,*)'gb7,q2,x,gb0'
     *,gb7,q2,x,gb0
      if(rangen().gt.gb7)goto 9

      if(ish.ge.6)write(ifch,*)'res: jj,iqc,ncc:',
     *jj,iqc(jj),ncc(1,jj),ncc(2,jj)

      nqc(2)=0
      iqnew=iqc(jj)
      if(jgamma.ne.1)then
        if(rangen().lt.sj1/(sj1+sj2))then
          if(iqc(jj).eq.0)then
            jt=1
            jq=int(1.5+rangen())
            nqc(1)=ncc(jq,jj)
          else
            jt=2
            if(iqc(jj).gt.0)then
              jq=1
            else
              jq=2
            endif
            nqc(1)=0
            iqnew=0
          endif
          iq1=iqc(jj)
        else
          if(iqc(jj).ne.0)then
            iq1=0
            jt=3
            if(iqc(jj).gt.0)then
              jq=1
            else
              jq=2
            endif
            nqc(1)=ncc(1,jj)
          else
            jt=4
            jq=int(1.5+rangen())
            iq1=int(naflav*rangen()+1.)*(3-2*jq)
            nqc(1)=ncc(jq,jj)
            iqnew=-iq1
          endif
        endif
      else
        jt=5
        jq=int(1.5+rangen())
        iq1=int(naflav*rangen()+1.)*(3-2*jq)
        iqnew=-iq1
        nqc(1)=0
      endif
      eprt=max(1.d0*qt,
     *.5d0*((1.d0-x)*wpt(jj)+qt2/(1.d0-x)/wpt(jj)))
      pl=((1.d0-x)*wpt(jj)-eprt)*(3-2*jj)
      if(iq1.eq.0)then
        iq2ini=9
      else
        iq2ini=iq1
      endif
27    call timsh1(q2,sngl(eprt),iq2ini)
      amprt=pprt(5,1)**2
      plprt=eprt**2-amprt-qt2
      if(plprt.lt.-1d-6)goto 27
      ep3(1)=eprt
      ep3(2)=dsqrt(max(0.d0,plprt))
      if(pl.lt.0.d0)ep3(2)=-ep3(2)
      ey(1)=1.
      ey(2)=1.
      ey(3)=1.
      do i=1,4
        ept1(i)=ept(i)-ep3(i)
      enddo
      call psdefrot(ep3,s0xh,c0xh,s0h,c0h)
      if(ish.ge.6)then
      write(ifch,*)'q2,amprt,qt2',q2,amprt,qt2
      write(ifch,*)'eprt,plprt',eprt,plprt
      write(ifch,*)'ep3',ep3
      write(ifch,*)'ept',ept
      write(ifch,*)'ept1',ept1
      endif
      s2new=psnorm(ept1)

      if(s2new.gt.s2min2)then
        if(jj.eq.1)then
          gb=psjti(q2,qq,s2new,iqnew,iqc(2),1)
        else
          gb=psjti1(q2,qmin(1),qq,s2new,iqnew,iqc(1),1)
        endif
        if(iqnew.eq.0)then
          gb=gb/sj1
        else
          gb=gb/sj2
        endif
        if(ish.ge.1)then
          if(gb.gt.1.)write (ifch,*)'gb,s2new,s2,q2,iqnew',
     *    gb,s2new,s2,q2,iqnew
        endif
        if(rangen().gt.gb)goto 9
      else
        goto 9
      endif
      jgamma=0

      call psreti(nqc,jq,1,ey,s0xh,c0xh,s0h,c0h)

      if(jt.eq.1)then
        ncc(jq,jj)=nqc(2)
      elseif(jt.eq.2)then
        ncc(jq,jj)=ncc(1,jj)
        ncc(3-jq,jj)=nqc(1)
      elseif(jt.eq.3)then
        ncc(1,jj)=nqc(2)
      elseif(jt.eq.4)then
        ncc(1,jj)=ncc(3-jq,jj)
      elseif(jt.eq.5)then
        ncc(1,jj)=nqc(1)
        ncc(2,jj)=0
      endif
      iqc(jj)=iqnew
      if(ish.ge.6)write(ifch,*)'qt2,amprt,ncc:',
     *qt2,amprt,ncc(1,jj),ncc(2,jj)

      do i=1,4       
        ept(i)=ept1(i)
      enddo        
c c.m. energy squared, minimal  4-momentum transfer square and gluon 4-v
c for the next ladder run
      qmin(jj)=q2
      si=s2new
      if(ish.ge.3)write (ifch,*)'res: new jet - iqj,ncj,ep3,ept',
     *iqj(nj),ncj(1,nj),ncj(2,nj),ep3,ept

      goto 8            !next simulation step will be considered

10    continue
      if(ish.ge.3)write(ifch,*)'res: iqc,si,ept:',iqc,si,ept

c highest virtuality subprocess in the ladder
c---------------------------------------
      qqs=max(qmin(1)/4.,4.*qmin(2))
      qqs=max(qqs,qq)
      call psabor(si,qqs,iqc,ncc,ept,1,nptlh,bx)
      goto 17
      
12    continue
c---------------------------------------
c hard pomeron (direct photon)
c---------------------------------------
      ept(1)=.5*(wpi+wmi)
      ept(2)=.5*(wpi-wmi)
      ept(3)=0.
      ept(4)=0.
      if(ish.ge.3)write (ifch,*)'direct photon - ept,si,qq:',ept,si,qq

13    continue

c ladder rung
c---------------------------------------
      pt2=ept(3)**2+ept(4)**2
      pt=sqrt(pt2)
      wpt(1)=ept(1)+ept(2)
      wpt(2)=si/wpt(1)

      gdbor=psdbin(qmin(2),qq,si,iqc(2),long)
      gdtot=psdsin(qmin(2),qq,si,iqc(2),long)
      if(iqc(2).ne.0)then
        if(ish.ge.8)write (ifch,*)'qmin(2),qq,si',qmin(2),qq,si
        gdnon=psdnsi(qmin(2),qq,si,long)
        if(iabs(iqc(2)).eq.1.or.iabs(iqc(2)).eq.4)then
          gdbor=gdbor/2.25
          gdtot=gdnon/2.25+gdtot/4.5
        else
          gdbor=gdbor/9.
          gdtot=gdnon/9.+gdtot/4.5
        endif
      else
        gdnon=0.
      endif

      if(long.ne.0.or.qmin(2).ge.qq)then
        s2min=qq+4.*max(qmin(2),qcmass**2)
        wwmin=(5.*s2min-qq)/4.-2.*pt*sqrt(q2ini)
        wwmin=(wwmin+sqrt(wwmin**2-(qq-pt2)*(s2min-qq-4.*q2ini)))
     *  /2./(1.-4.*q2ini/(s2min-qq))
      else
        s2min=qq/(1.-sqrt(q2ini/qq))
        wwmin=s2min+qq-2.*pt*sqrt(q2ini)
        wwmin=(wwmin+sqrt(wwmin**2-4.*(qq-pt2)*(qq-q2ini)))
     *  /2./(1.-q2ini/qq)
      endif
      
      if(ish.ge.3)write(ifch,*)'si,s2min,wwmin,qmin(2),gdtot,gdbor:'
      if(ish.ge.3)write(ifch,*)si,s2min,wwmin,qmin(2),gdtot,gdbor

      if((rangen().lt.gdbor/gdtot.or.si.lt.1.1*wwmin).and.
     *(long.eq.0.and.qmin(2).lt.qq.or.iqc(2).eq.0))goto 15
      if(si.lt.1.1*wwmin)stop'si<1.1*wwmin'

      qqmax=0.
      qqmin=0.

      xmm=si+2.*sqrt(q2ini)*pt-qq
      discr=xmm**2-4.*q2ini*(5.*si-qq+pt2)
      if(discr.lt.0.)goto 29
      discr=sqrt(discr)
      qqmax=(xmm+discr)/2./(5.-(qq-pt2)/si)
      qqmin=2.*q2ini*si/(xmm+discr)
      
29    continue
      if(4.*qqmin.lt.s2min-qq.or.long.eq.0.and.
     *qmin(2).lt.qq)then
        xmm=si-s2min+2.*sqrt(q2ini)*pt
        qqmin=2.*q2ini*si/(xmm+sqrt(xmm**2-4.*q2ini*(si-qq+pt2)))
      endif
      xmin=1.-q2ini/qqmin

      if(qqmin.lt.qmin(2))then
        qqmin=qmin(2)
        xmi=max(1.-((pt*sqrt(qqmin)+sqrt(pt2*qqmin+
     *  si*(si-s2min-qqmin*(1.-(qq-pt2)/si))))/si)**2,
     *  (s2min+qqmin*(1.-(qq-pt2)/si)-2.*pt*sqrt(qqmin))/si)
        xmin=max(xmin,xmi)
      endif
      if(xmin.le.qq/si)xmin=1.001*qq/si

      if(long.eq.0.and.qmin(2).lt.qq)qqmax=max(qqmax,qq)
      xmax=1.-q2ini/qqmax
      
      if(ish.ge.6)write(ifch,*)'qqmax,qqmin,xmax,xmin:',
     *qqmax,qqmin,xmax,xmin
      if(qqmax.lt.qqmin)stop'qqmax<qqmin'

      qm0=qqmin
      xm0=1.-q2ini/qm0
      s2max=si*xm0-qm0*(1.-qq/si)

      sds=psdsin(qm0,qq,s2max,0,long)/4.5
      sdv=psdsin(qm0,qq,s2max,1,long)/4.5

      sdn=psdnsi(qm0,qq,s2max,long)
      if(iqc(2).eq.0)then
        sdn=sdn/4.5
      elseif(iabs(iqc(2)).eq.1.or.iabs(iqc(2)).eq.4)then
        sdn=sdn/2.25
      else
        sdn=sdn/9.
      endif
      sdv=sdv+sdn
      xx=xm0

      sj0=sds*psfap(xx,iqc(2),0)+sdv*psfap(xx,iqc(2),1)
      gb0=sj0/log(q2ini/qcdlam)*sngl(psuds(qm0,iqc(2)))*qm0*5.
      if(gb0.le.0.)then
         write(ifmt,*)'gb0.le.0.  si,qq,pt2:',si,qq,pt2
         iret=1
         goto 9999
      endif

      if(xm0.le..5)then
        gb0=gb0*(xm0-qq/si)/(1.-2.*qq/si)
      else
        gb0=gb0*(1.-xm0)
      endif

      xmin2=max(.5,xmin)
      xmax1=min(xmax,.5)
      if(xmin.ge..5)then
        djl=1.
      elseif(xmax.lt..5)then
        djl=0.
      else
        djl=1./(1.-(1.-2.*qq/si)*log((.5-qq/si)/(xmin-qq/si))/
     *  log(2.*(1.-xmax)))
      endif

14    continue
      if(rangen().gt.djl)then        !lc momentum share in the cur
        x=(xmin-qq/si)*((xmax1-qq/si)/(xmin-qq/si))**rangen()+qq/si
      else
        x=1.-(1.-xmin2)*((1.-xmax)/(1.-xmin2))**rangen()
      endif

      q2=qqmin/(1.+rangen()*(qqmin/qqmax-1.))

      qt2=q2*(1.-x)
      if(ish.ge.9)write(ifch,*)'q2,x,qt2,qq,qqmin,qqmax:',
     *q2,x,qt2,qq,qqmin,qqmax
      if(qt2.lt.q2ini)goto 14   !p_t check

      if(long.ne.0.or.q2.ge.qq)then
        s2min2=max(4.*q2+qq,s2min)
      else
        s2min2=s2min
      endif
      qt=sqrt(qt2)
      call pscs(bcos,bsin)
c ep3 is now 4-vector for s-channel gluon produced in current ladder run
      ep3(3)=qt*bcos
      ep3(4)=qt*bsin
      ptnew=(ept(3)-ep3(3))**2+(ept(4)-ep3(4))**2

      s2=x*si-ptnew+pt2-q2*(x-(qq-pt2)/si)
      if(s2.lt.s2min2)goto 14   !check of the kinematics
      sds=psdsin(q2,qq,s2,0,long)/4.5
      sdv0=psdsin(q2,qq,s2,1,long)/4.5
      if(ish.ge.8)write (ifch,*)'q2,qq,s2',q2,qq,s2
      sdn0=psdnsi(q2,qq,s2,long)

      if(iqc(2).eq.0)then
        sdn=sdn0/4.5
      else
        if(iabs(iqc(2)).eq.1.or.iabs(iqc(2)).eq.4)then
          sdn=sdn0/2.25
        else
          sdn=sdn0/9.
        endif
      endif
      sdv=sdv0+sdn

      xx=x
      sj1=sds*psfap(xx,iqc(2),0)
      sj2=sdv*psfap(xx,iqc(2),1)

c gb7 is the rejection function for x and q**2 simulation.
      gb7=(sj1+sj2)/log(qt2/qcdlam)*sngl(psuds(q2,iqc(2)))/gb0*q2
      if(x.le..5)then
        gb7=gb7*(x-qq/si)/(1.-2.*qq/si)
      else
        gb7=gb7*(1.-x)
      endif
 
      if(gb7.gt.1..and.ish.ge.1)write(ifmt,*)'gb7,q2,x,qt2,iqc(2),'
     * ,'gb0,sj1,sj2',gb7,q2,x,qt2,iqc(2),gb0,sj1,sj2
      if(ish.ge.3)write (ifch,*)'gb7,q2,x,qt2,iqc(2),gb0,sj1,sj2,long',
     * gb7,q2,x,qt2,iqc(2),gb0,sj1,sj2,long
      if(rangen().gt.gb7)goto 14


      if(ish.ge.6)write(ifch,*)'iqc,ncc:',iqc(2),ncc(1,2),ncc(2,2)
      iqcnew=iqc(2)
      nqc(2)=0         !emitted parton color connections
      if(rangen().lt.sj1/(sj1+sj2).or.(long.ne.0.or.q2.ge.qq).and.
     *s2.lt.1.5*s2min2)then
        if(iqc(2).eq.0)then
          jt=1
          jq=int(1.5+rangen())
          nqc(1)=ncc(jq,2)
        else
          jt=2
          if(iqc(2).gt.0)then
            jq=1
          else
            jq=2
          endif
          nqc(1)=0
        endif
        iq1=iqc(2)
        iqcnew=0

      else
        if(iqc(2).ne.0)then
          jt=3
          iq1=0
          if(iqc(2).gt.0)then
            jq=1
          else
            jq=2
          endif
          nqc(1)=ncc(1,2)

        else
          tu=sdn0/2.25+sdv0
          if(naflav.eq.4)tu=tu*2.
          td=sdn0/9.+sdv0
          if(rangen().lt.tu/(tu+2.*td))then
            if(naflav.eq.3)then
              iq1=1
            else
              iq1=1+3*int(.5+rangen())
            endif
          else
            iq1=int(2.5+rangen())
          endif
          jq=int(1.5+rangen())
          iq1=iq1*(3-2*jq)
          iqcnew=-iq1
          jt=4
          nqc(1)=ncc(jq,2)
        endif
      endif

      eprt=max(1.d0*qt,
     *.5d0*((1.d0-x)*wpt(2)+qt2/(1.d0-x)/wpt(2)))
      pl=eprt-(1.d0-x)*wpt(2)
      if(iq1.eq.0)then
        iq2ini=9
      else
        iq2ini=iq1
      endif
28    call timsh1(q2,sngl(eprt),iq2ini)
      amprt=pprt(5,1)**2
      plprt=eprt**2-amprt-qt2
      if(plprt.lt.-1d-6)goto 28
      ep3(1)=eprt
      ep3(2)=dsqrt(max(0.d0,plprt))
      if(pl.lt.0.d0)ep3(2)=-ep3(2)
      ey(1)=1.
      ey(2)=1.
      ey(3)=1.
      do i=1,4
        ept1(i)=ept(i)-ep3(i)
      enddo
      call psdefrot(ep3,s0xh,c0xh,s0h,c0h)
      call psrotat(ep3,s0xh,c0xh,s0h,c0h)
      s2new=psnorm(ept1)+qq

      if((long.ne.0.or.q2.ge.qq).and.iqcnew.ne.0)then
        xmm=(5.*s2min2-qq)/4.-2.*sqrt(ptnew*q2ini)
        s2min2=1.1*(xmm+sqrt(xmm**2-(qq-ptnew)*
     *  (s2min2-qq-4.*q2ini)))/2./(1.-4.*q2ini/(s2min2-qq))
      endif
      if(s2new.gt.s2min2)then
        sds1=psdsin(q2,qq,s2new,iqcnew,long)/4.5
        if(iqcnew.eq.0)then
          gb=sds1/sds
        else
          if(ish.ge.8)write (ifch,*)'q2,qq,s2new',q2,qq,s2new
          sdn1=psdnsi(q2,qq,s2new,long)
          if(iabs(iqcnew).eq.1.or.iabs(iqcnew).eq.4)then
            sdn1=sdn1/2.25
            sdv=sdv0+sdn0/2.25
          else
            sdn1=sdn1/9.
            sdv=sdv0+sdn0/9.
          endif
          gb=.9999*(sds1+sdn1)/sdv
        endif
        if(ish.ge.3.and.gb.gt.1..and.ish.ge.1)write(ifmt,*)'gbs2',gb
        if(rangen().gt.gb)goto 14
      else
        goto 14
      endif

      call psreti(nqc,jq,1,ey,s0xh,c0xh,s0h,c0h)

      iqc(2)=iqcnew
      if(jt.eq.1)then      !current parton color connections
        ncc(jq,2)=nqc(2)
      elseif(jt.eq.2)then
        ncc(jq,2)=ncc(1,2)
        ncc(3-jq,2)=nqc(1)
      elseif(jt.eq.3)then
        ncc(1,2)=nqc(2)
      elseif(jt.eq.4)then
        ncc(1,2)=ncc(3-jq,2)
        ncc(2,2)=0
      endif

      do i=1,4
        ept(i)=ept1(i)
      enddo
      if(ish.ge.3)write (ifch,*)'new jet - iqj,ncj,ep3,ept',
     *iqj(nj),ncj(1,nj),ncj(2,nj),ep3,ept
c c.m. energy squared, minimal  4-momentum transfer square and gluon 4-v
c for the next ladder run
      qmin(2)=q2
      si=s2new
      goto 13            !next simulation step will be considered

15    continue
      if(ish.ge.3)write(ifch,*)'iqc,si,qmin(2),nj:',
     *iqc(2),si,qmin(2),nj
c highest virtuality subprocess in the ladder
c---------------------------------------
      gb01=0.
      tmax=0.
      tmin=si
      if(iqc(2).eq.0.and.si.gt.qq+4.*max(qcmass**2,qmin(2)))then
        qminn=max(qcmass**2,qmin(2))
        tmin1=2.*qminn/(1.-qq/si)/(1.+sqrt(1.-4.*qminn/(si-qq)))
        tmin=tmin1
        tmax=si/2.
        fb01=psdbom(si,si/2.,si/2.,qq,long)
        if(long.eq.0)fb01=fb01*si/2.
        gb01=fb01/log(qminn/qcdlam)*sngl(psuds(qminn,iqc(2)))/si**2
        gb0=gb01
      endif
      
      if(long.eq.0.and.qmin(2).lt.qq)then
        tmax=max(tmax,qq) 
        tmin=max(qmin(2),
     *  2.*q2ini/(1.-qq/si)/(1.+sqrt(1.-4.*q2ini/(si-qq))))
        ze=qq/si+tmin/si*(1.-qq/si)
        xx=ze
        qt2=tmin*(1.-ze)
        if(qt2.lt..999*q2ini.and.ish.ge.1)write(ifmt,*)'bor-dir:qt20'
     *                                                 ,qt2
        gb0=gb01+psfap(xx,iqc(2),1)/log(qt2/qcdlam)
     *  *sngl(psuds(tmin,iqc(2))/psuds(tmin,1)*psuds(qq,1))
     *  /si*(1.-tmin*qq/si**2/ze)
      endif
      gb0=gb0*2.

      call psdeftr(si-qq,ept,ey)

      if(ish.ge.6)write(ifch,*)'tmin,tmax,qq,si-qq,gb0:'
      if(ish.ge.6)write(ifch,*)tmin,tmax,qq,si-qq,psnorm(ept),gb0

c------------------------------------------------
16    continue
      if(long.eq.0)then
        t=tmin*(tmax/tmin)**rangen()
      else
        t=tmin+(tmax-tmin)*rangen()
      endif
      
      u=si-t
      ze=qq/si+t/si*(1.-qq/si)
      qt2=t*(1.-ze)
      if(t.le.qq.and.long.eq.0)then
        xx=ze
        gb=psfap(xx,iqc(2),1)/log(qt2/qcdlam)*sngl(psuds(t,iqc(2))
     *  /psuds(t,1)*psuds(qq,1))/si*(1.-t*qq/si**2/ze)/gb0
      else
        gb=0.
      endif
      
      gb1=0.
      if(iqc(2).eq.0..and.si.gt.qq+4.*max(qcmass**2,qmin(2)).
     *and.qt2.gt.qcmass**2.and.t.le.si/2..and.t.ge.tmin1)then
        fb1=psdbom(si,t,u,qq,long)
        if(long.eq.0)fb1=fb1*t
        gb1=fb1/log(qt2/qcdlam)*sngl(psuds(qt2,iqc(2)))/si**2/gb0
c        gb1=0.  !???????????????????????
        gb=gb+gb1
      endif
      

      if(ish.ge.6)write(ifch,*)'gb,t,iqc(2),si,qq,qmin(2),long:',
     *gb,t,iqc(2),si,qq,qmin(2),long
      if (ish.ge.1) then  
        if(gb.gt.1.)write(*,*)'gb,gb1,gb0,gb01',
     *  ',t,iqc(2),si,qq,qmin(2),long:',
     *  gb,gb1,gb0,gb01,fb1,fb01,t,iqc(2),si,qq,qmin(2),long
      endif 
     
      if(rangen().gt.gb)goto 16
      if(ish.ge.3)write(ifch,*)'born:t,qt2:',t,qt2

      nqc(2)=0
      if(iqc(2).eq.0)then
        jq=int(1.5+rangen())
        jq2=3-jq
        if(rangen().gt.gb1/gb)then
          iq1=(1+int(3.*rangen()))*(3-2*jq)
        else
          iq1=4*(3-2*jq)
        endif
        iq2=-iq1                             !quark flavors
        nqc(1)=ncc(jq,2)
      else
        if(iqc(2).gt.0)then
          jq=1
        else
          jq=2
        endif
        jq2=jq
        iq1=0
        iq2=iqc(2)
        nqc(1)=ncc(1,2)
      endif

      call pscs(bcos,bsin)
      z=sngl(psutz(dble(si-qq),dble(qt2),dble(qt2)))
      if(t.lt..5*si)z=1.-z
      wp3=z*sqrt(si-qq)
      wm3=qt2/wp3
      if(iabs(iq1).eq.4)qt2=qt2-qcmass**2
      qt=sqrt(qt2)
      ep3(1)=.5*(wp3+wm3)
      ep3(2)=.5*(wp3-wm3)
      ep3(3)=qt*bcos
      ep3(4)=qt*bsin
      call psdefrot(ep3,s0xh,c0xh,s0h,c0h)
      if(iq1.eq.0)then
        iq2ini1=9
      else
        iq2ini1=iq1
      endif
      if(iq2.eq.0)then
        iq2ini2=9
      else
        iq2ini2=iq2
      endif
      if(ish.ge.5)write (ifch,*)'jq,jt,iq2ini1,iq2ini2',
     *jq,jt,iq2ini1,iq2ini2
     
      if(t.lt.qq.and.iabs(iq1).ne.4)then
        qq1=t*(1.-ze)
        qq2=qq
      else
        qq1=qt2
        qq2=qt2
      endif
      call timsh2(qq1,qq2,sqrt(si-qq),iq2ini1,iq2ini2)
      nfprt=1
      call psreti(nqc,jq,nfprt,ey,s0xh,c0xh,s0h,c0h)

      if(iqc(2).eq.0)then
        nqc(1)=ncc(3-jq,2)
        nqc(2)=0
      else
        nqc(1)=nqc(2)
        nqc(2)=0
      endif

      nfprt=2
      call psreti(nqc,jq2,nfprt,ey,s0xh,c0xh,s0h,c0h)

17    continue
      if(ish.ge.3)write (ifch,*)'nj',nj
      if(nj.gt.0)then
                   
          ityj(i)=30   
          iorj(i)=nptlh   
        do n=1,nj
          do i=1,4
            ep3(i)=eqj(i,n)
          enddo
          call pstrans(ep3,ey0,-1)         !boost to the c.m. system
          do i=1,4
            eqj(i,n)=ep3(i)
          enddo
          do l=1,6
            bxj(l,n)=bx(l)
          enddo
          ityj(n)=0
          iorj(n)=1
        enddo
      endif
      call psjarr(jfl)       !kinky strings formation
      if(ish.ge.3)write (ifch,*)'jfl',jfl
      if(jfl.eq.0)then
        iret=1
      else
        iret=0
        ep3(4)=egyevt
        ep3(2)=0.
        ep3(3)=0.
        ep3(1)=0.
        do i=2,nptl
          do l=1,4
            ep3(l)=ep3(l)-pptl(l,i)
          enddo
        enddo
        if(ish.ge.3)write(ifch,*)'energy-momentum balance:'
        if(ish.ge.3)write(ifch,*)ep3
        if(abs(ep3(4)).gt.3.e-2)write(*,*)'energy-momentum balance:',ep3
      endif
 9999 call utprix('psadis',ish,ishini,3)
      return
      end
      
c-----------------------------------------------------------------------
      subroutine psaevc
c-----------------------------------------------------------------------
      include 'epos.inc'
c structure functions calculation
      logical lcalc
      double precision xx,xmax
      dimension evs(21,21,135)
      common /psar2/  edmax,epmax
      common /psar31/ evk0(21,21,54)
      common /psar32/ evk(21,21,135)
      include 'epos.incsem'

      inquire(file=fnie(1:nfnie),exist=lcalc)
      if(lcalc)then
       if(inicnt.eq.1)then
        write(ifmt,'(3a)')'read from ',fnie(1:nfnie),' ...'
        open(1,file=fnie(1:nfnie),status='old')
        read (1,*)qcdlam0,q2min0,q2ini0,naflav0,epmax0
        if(qcdlam0.ne.qcdlam)write(ifmt,'(a)')'iniev: wrong qcdlam'
        if(q2min0 .ne.q2min )write(ifmt,'(a)')'iniev: wrong q2min'
        if(q2ini0 .ne.q2ini )write(ifmt,'(a)')'iniev: wrong q2ini'
        if(naflav0.ne.naflav)write(ifmt,'(a)')'iniev: wrong naflav'
        if(epmax0 .ne.epmax )write(ifmt,'(a)')'iniev: wrong epmax'
        if(qcdlam0.ne.qcdlam.or.q2min0.ne.q2min.or.q2ini0.ne.q2ini
     *  .or.naflav0.ne.naflav.or.epmax0.ne.epmax)then
           write(6,'(//a//)')'   iniev has to be reinitialized!!!'
           stop
        endif   
        read (1,*)evk0,evk
        close(1)
       endif
       goto 101
      endif

      write(ifmt,'(a)')'iniev does not exist -> calculate tables  ...'
      xmax=1.d0-2.d0*q2ini/epmax
      do l=1,27
        if(l.le.12)then
          xx=.1d0*exp(l-13.d0)
        elseif(l.le.21)then
          xx=.1d0*(l-12.d0)
        else
          xx=1.d0-.1d0*(10.d0*(1.d0-xmax))**((l-21)/6.)
        endif

        qmin=max(1.d0*q2min,q2ini/(1.-xx))
      do i=1,21
        qq=qmin*(.5*epmax/qmin)**((i-1)/20.)
      do j=1,21
        qj=qmin*(qq/qmin)**((j-1)/20.)
        if(l.eq.27.or.i.eq.1.or.j.eq.21)then
          evk0(i,j,l)=0.
          evk0(i,j,l+27)=0.
          do k=1,5
            evk(i,j,l+27*(k-1))=0.
          enddo
        else
          do k=1,2
            evk0(i,j,l+27*(k-1))=log(psev0(qj,qq,xx,k))
          enddo
        endif
      enddo
      enddo
      enddo

      n=1

1     n=n+1
      write(ifmt,2)n
2     format(5x,i2,'-th order contribution')

      do l=1,26
        write(ifmt,*)'l',l
        if(l.le.12)then
          xx=.1d0*exp(l-13.d0)
        elseif(l.le.21)then
          xx=.1d0*(l-12.d0)
        else
          xx=1.d0-.1d0*(10.d0*(1.d0-xmax))**((l-21)/6.)
        endif

        qmin=max(1.d0*q2min,q2ini/(1.d0-xx))
      do i=2,21
        qq=qmin*(.5*epmax/qmin)**((i-1)/20.)
      do j=1,20
        qj=qmin*(qq/qmin)**((j-1)/20.)
        do m=1,3
        do k=1,2
          if(m.ne.3)then
            ev=psev(qj,qq,xx,m,k,n)
            ev0=psevi0(qj,qq,xx,m,k)
            evs(i,j,l+27*(m-1)+54*(k-1))=log((ev+ev0)/psfap(xx,m-1,k-1)
     *      /log(log(qq*(1.d0-xx)/qcdlam)/log(qj*(1.d0-xx)/qcdlam))*4.5)
          elseif(k.ne.1)then
            evs(i,j,l+108)=log((psev(qj,qq,xx,m,k,n)+
     *      psevi0(qj,qq,xx,2,2))/psfap(xx,2,2)
     *      /log(log(qq*(1.d0-xx)/qcdlam)/log(qj*(1.d0-xx)/qcdlam))*4.5)
          endif
        enddo
        enddo
      enddo
      enddo
      enddo

      jec=0
      do i=2,21
      do j=1,20
      do l=1,26
      do k=1,5
        if(n.eq.2.or.evs(i,j,l+27*(k-1)).ne.0..and.
     *  abs(1.-evk(i,j,l+27*(k-1))/evs(i,j,l+27*(k-1))).gt.1.e-2)then
          jec=1
          evk(i,j,l+27*(k-1))=evs(i,j,l+27*(k-1))
        endif
      enddo
      enddo
      enddo
      enddo

      if(jec.ne.0)goto 1

      write(ifmt,'(a)')'write to iniev ...'
      open(1,file=fnie(1:nfnie),status='unknown')
      write (1,*)qcdlam,q2min,q2ini,naflav,epmax
      write (1,*)evk0,evk
      close(1)

101   continue
      return
      end

c------------------------------------------------------------------------
      function psdbom(s,t,u,qq,long)
c-----------------------------------------------------------------------
c psdbom - integrand for DIS c-quark cross-sections (matrix element squared)
c s  - total c.m. energy squared for the scattering (for n=2: s+qq),
c t  - invariant variable for the scattering |(p1-p3)**2|
c u  - invariant variable for the scattering |(p1-p4)**2|
c qq - photon virtuality
c long: 0 - contr. to (F2-F_L), 1 - contr. to F_L
c-----------------------------------------------------------------------
      include 'epos.incsem'
      if(long.eq.0)then       !F2-F_L
        psdbom=(2.*(t/u+u/t)*(qq**2+(s-qq)**2)/s**2+
     *  4.*(qcmass*s/t/u)**2*(qq-2.*qcmass**2)+
     *  8.*qcmass**2/t/u*(s-2.*qq))   *2.    !=4.5/2.25
      else                    !F_L_C
        psdbom=16.*qq*((s-qq)/s**2-qcmass**2/t/u) *2.  !=4.5/2.25
      endif
      return
      end

c------------------------------------------------------------------------
      function psdbin(q1,qq,s,m1,long)
c-----------------------------------------------------------------------
c psdbin - DIS born cross-section
c q1      - virtuality cutoff for current end of the ladder
c qq      - photon virtuality
c s=2(pq) - s_true + qq
c s2min   - mass cutoff for born scattering
c m1       - incoming parton type (0 - g, 1,2 - q)
c-----------------------------------------------------------------------
      double precision xx
      include 'epos.incsem'
      include 'epos.inc'

      psdbin=0.
      q2mass=qcmass**2
      s2min=4.*max(q1,q2mass)+qq
      if(m1.eq.0.and.s.gt.s2min.and.(idisco.eq.0.or.idisco.eq.2))then
        tmax=s/2.
        qtq=4.*max(q2mass,q1)/(s-qq)
        if(qtq.lt.1.)then
          tmin=.5*s*qtq/(1.+sqrt(1.-qtq))
        else
          tmin=.5*s
        endif
        psdbin=psdbin+psdbor(q1,qq,s,long)*(1./tmin-1./tmax)
      endif
     
      if(long.eq.0.and.q1.lt.qq.and.s.gt.qq/(1.-q2ini/qq)
     *.and.(idisco.eq.0.or.idisco.eq.1))then
        m=min(1,iabs(m1))+1
        xx=qq/s
        psdbin=psdbin+psevi0(q1,qq,xx,m,2)*4.*pi**2*alfe/s
      endif
      return
      end

c------------------------------------------------------------------------
      function psdbor(q1,qq,s,long)
c-----------------------------------------------------------------------
c psdbor - DIS born cross-section
c q1      - virtuality cutoff for current end of the ladder
c qq      - photon virtuality
c s=2(pq) - s_true + qq
c s2min   - mass cutoff for born scattering
c-----------------------------------------------------------------------
      common /ar3/    x1(7),a1(7)
      include 'epos.inc'
      include 'epos.incsem'
      double precision psuds

      psdbor=0.
      q2mass=qcmass**2
      qtq=4.*max(q2mass,q1)/(s-qq)
      j=0   !Gluon

      tmax=s/2.
      if(qtq.lt.1.)then
        tmin=.5*s*qtq/(1.+sqrt(1.-qtq))
      else
        tmin=.5*s
      endif
      if(tmax.lt.tmin.and.ish.ge.1)write(ifmt,*)'s,q1,qq,tmin,tmax',
     *s,q1,qq,tmin,tmax

      ft=0.
      do i=1,7
      do m=1,2
        t=2.*tmin/(1.+tmin/tmax+(2*m-3)*x1(i)*(1.-tmin/tmax))
        u=s-t

        qt=t*u/s*(1.-qq/s)
        if(qt.lt..999*max(q2mass,q1).and.ish.ge.1)
     &  write(ifmt,*)'psdbor:qt,q1',qt,q1
        fb=psdbom(s,t,u,qq,long)*t**2
        ft=ft+a1(i)*fb*pssalf(qt/qcdlam)*sngl(psuds(qt,j))
      enddo
      enddo
      psdbor=ft/s**2*pi**2*alfe/sngl(psuds(q1,j))
      return
      end

c------------------------------------------------------------------------
      subroutine psdint(s,qq,sds,sdn,sdb,sdt,sdr,m1,long)
c-----------------------------------------------------------------------
c psdint - dis cross-sections interpolation - for minimal
c effective momentum cutoff in the ladder
c s   - total c.m. energy squared for the ladder,
c qq  - photon virtuality,
c sds - dis singlet cross-section,
c sdn - dis nonsinglet cross-section,
c sdb - dis born cross-section,
c sdt - dis singlet+resolved cross-section,
c m1  - parton type at current end of the ladder (0 - g, 1,2 - q)
c-----------------------------------------------------------------------
      double precision xx
      dimension wk(3),wj(3)
      common /psar2/  edmax,epmax
      common /psar27/ csds(21,26,4),csdt(21,26,2),csdr(21,26,2)
      include 'epos.incsem'
      include 'epos.inc'
      
      sds=0.
      sdn=0.
      sdt=0.
      sdr=0.
      sdb=psdbin(q2min,qq,s,m1,long)

      m=min(1,iabs(m1))+1
      qlj=log(qq/q2min)*2.+1.
      j=int(qlj)
      if(j.lt.1)j=1
      if(j.gt.19)j=19
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.)*.5
      wj(1)=1.-wj(2)+wj(3)
      wj(2)=wj(2)-2.*wj(3)

      s2min=4.*max(q2min,qcmass**2)+qq
      if(m1.ne.0)s2min=s2min/(1.-4.*q2ini/(s2min-qq))
      if(s.le.s2min.or.idisco.ne.0.and.idisco.ne.2)goto 1
      
      qtq=4.*max(q2min,qcmass**2)/(s-qq)
      if(qtq.lt.1.)then
        tmin=.5*s*qtq/(1.+sqrt(1.-qtq))
      else
        tmin=.5*s
      endif
      tmax=s/2.

      sl=log(s/s2min)/log(edmax/s2min)*25.+1.
      k=int(sl)
      if(k.lt.1)k=1
      if(k.gt.24)k=24
      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      do k1=1,3
        k2=k+k1-1
      do j1=1,3
        sds=sds+csds(j+j1-1,k2,m+2*long)*wj(j1)*wk(k1)
      enddo
      enddo
      if(m.eq.1)then
        sds=exp(sds)*(1./tmin-1./tmax)
      else
        sds=max(sds,0.)     
      endif
      
1     continue
      s2min=max(4.*qq,16.*q2min)+qq
      if(s.le.s2min.or.long.ne.0.or.idisco.ne.0.and.idisco.ne.3)then
        sdt=sds
        goto 2
      endif

      sl=log(s/s2min)/log(edmax/s2min)*25.+1.
      k=int(sl)
      if(k.lt.1)k=1
      if(k.gt.24)k=24
      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)

      do k1=1,3
        k2=k+k1-1
      do j1=1,3
        sdr=sdr+csdr(j+j1-1,k2,m)*wj(j1)*wk(k1)
        sdt=sdt+csdt(j+j1-1,k2,m)*wj(j1)*wk(k1)
      enddo
      enddo
      
      sdr=max(sdr,0.)
      sdt=max(sds,sds+sdt) 
      sdt=sdt+sdr
            
2     continue
      if(long.eq.0.and.q2min.lt.qq.and.s.gt.qq/(1.-q2ini/qq)
     *.and.(idisco.eq.0.or.idisco.eq.1))then
        xx=qq/s
        dsi=psevi(q2min,qq,xx,m,2)*4.*pi**2*alfe/s
        if(m1.eq.0)then
          sds=sds+dsi
          sdt=sdt+dsi
        else
          dnsi=psevi(q2min,qq,xx,3,2)*4.*pi**2*alfe/s
          sdn=sdn+dnsi
          sds=sds+max(dsi-dnsi,0.)
          sdt=sdt+max(dsi-dnsi,0.)
        endif
      endif
      
      if(m1.eq.0)then
        sds=max(sds,sdb)
        sdt=max(sdt,sdb)
      else
        sdn=max(sdn,sdb)
      endif
      return
      end

c-----------------------------------------------------------------------
      function psdnsi(q1,qq,s,long)
c-----------------------------------------------------------------------
c psdnsi - DIS nonsinglet cross-section interpolation
c q1 - effective momentum cutoff for current end of the ladder,
c qq - photon virtuality,
c s - total c.m. energy squared for the ladder,
c-----------------------------------------------------------------------
      double precision xx
      include 'epos.incsem'
      include 'epos.inc'

      psdnsi=0.
      if(long.eq.0.and.q1.lt.qq.and.s.gt.qq/(1.-q2ini/qq))then
        xx=qq/s
        psdnsi=psdnsi+max(0.,psevi(q1,qq,xx,3,2)*4.*pi**2*alfe/s)
      endif
      return
      end

c-----------------------------------------------------------------------
      function psdrga(qq,s,s2min,j)
c-----------------------------------------------------------------------
c psdrga - DIS resolved cross-section (photon sf)
c qq    - photon virtuality
c s     - total c.m. energy squared for the process
c s2min - mass cutoff for born scattering
c j     - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------
      common /ar3/   x1(7),a1(7)
      include 'epos.incsem'

      psdrga=0.
      if(s.le.s2min)return

      xmin=s2min/s
      do i=1,7          
      do m=1,2
        z=xmin**(.5+(m-1.5)*x1(i))
        tu=psdfh4(z,q2min,qq,0,1)
        td=psdfh4(z,q2min,qq,0,2)
        ts=psdfh4(z,q2min,qq,0,3)
        tg=psdfh4(z,q2min,qq,0,0)
        if(j.eq.0)then
          sj=tg*psjti(q2min,qq,z*s,0,j,1)+
     *    (tu+td+ts)*psjti(q2min,qq,z*s,1,j,1)
        else
          sj=tg*psjti(q2min,qq,z*s,0,j,1)+
     *    (tu+td)*(psjti(q2min,qq,z*s,1,1,1)/4.+
     *    psjti(q2min,qq,z*s,-1,1,1)/4.+
     *    psjti(q2min,qq,z*s,2,1,1)/2.)+
     *    ts*psjti(q2min,qq,z*s,2,1,1)
        endif
        psdrga=psdrga+a1(i)*sj
      enddo
      enddo
      psdrga=-psdrga*log(xmin)*alfe/2.  *4.5 !mean e^2 is taken out
      return
      end

c-----------------------------------------------------------------------
      function psdres(qq,s,s2min,j)
c-----------------------------------------------------------------------
c psdres - DIS resolved photon cross-section
c qq    - photon virtuality
c s     - total w squared for the ladder (s+qq)
c s2min - mass cutoff for born scattering
c j     - parton type at current end of the ladder (0 - g, 1,2 etc. - q)
c-----------------------------------------------------------------------
      double precision xx
      common /ar3/   x1(7),a1(7)
      include 'epos.inc'
      include 'epos.incsem'

      psdres=0.
      if(s.le.s2min+qq)return

      qmin=max(q2min,s2min/(s/qq-1.))
      qmax=min(s-s2min,s/2.)

c numerical integration over transverse momentum squared;
c gaussian integration is used
      do i=1,7
      do m=1,2
        qi=2.*qmin/(1.+qmin/qmax+(2*m-3)*x1(i)*(1.-qmin/qmax))

        zmax=min(1.,qi/qq)
        zmin=(max(qi,s2min)+qi)/s

        fsj=0.
        if(zmax.gt.zmin)then
          do i1=1,7
          do m1=1,2
            z=.5*(zmax+zmin+(2*m1-3)*x1(i1)*(zmax-zmin))
            s2=z*s-qi
            xx=z
            if(j.eq.0)then
              sj=psfap(xx,0,1)*psjti(qi,qq,s2,1,j,1)
            else
              sj=psfap(xx,0,1)*(psjti(qi,qq,s2,1,1,1)/6.+
     *        psjti(qi,qq,s2,-1,1,1)/6.+
     *        psjti(qi,qq,s2,2,1,1)/1.5)
            endif
            fsj=fsj+a1(i1)*sj*qi    !????????(qi-z*qq)
          enddo
          enddo
          fsj=fsj*(zmax-zmin)
        elseif(ish.ge.1)then
          write(ifmt,*)'psdres:zmax,zmin',zmax,zmin
        endif
        psdres=psdres+a1(i)*fsj
      enddo
      enddo
      psdres=psdres*(1./qmin-1./qmax)*alfe*.75/pi  !alpha_s -> 6 alpha_e
      return
      end

c------------------------------------------------------------------------
      function psds(q1,qq,s,j,long)
c-----------------------------------------------------------------------
c psds - DIS singlet cross-section
c q1      - virtuality cutoff for current end of the ladder
c qq      - photon virtuality
c s=2(pq) - s_true + qq
c s2min   - mass cutoff for born scattering
c-----------------------------------------------------------------------
      double precision xxe,xmax,xmin,xmax1,xmin1
      common /ar3/    x1(7),a1(7)
      include 'epos.inc'
      include 'epos.incsem'

      psds=0.
      q2mass=qcmass**2
      s2min=4.*max(q1,q2mass)
      smin=(s2min+qq)/(1.-4.*q2ini/s2min)
      if(s.le.1.001*smin)return

      xmax=.5d0*(1.d0+qq/s+dsqrt((1.d0-qq/s)**2-16.d0*q2ini/s))
      xmin=max(1.d0+qq/s-xmax,1.d0*(s2min+qq)/s)
      if(xmin.gt.xmax.and.ish.ge.1)write(ifmt,*)'xmin,xmax,q1,qq,s,smin'
     *,xmin,xmax,q1,qq,s,smin

      fx1=0.
      fx2=0.
      if(xmax.gt..9d0)then
        xmin1=max(xmin,.9d0)
        do i=1,7
        do m=1,2
          xxe=1.d0-(1.d0-xmax)*((1.d0-xmin1)/(1.d0-xmax))**
     *    (.5d0-x1(i)*(m-1.5))
          xx=xxe

          sh=xx*s
          qtmin=max(1.d0*max(q2mass,q1),q2ini/(1.d0-xxe))
          qtq=4.*qtmin/(sh-qq)

          tmin=.5*sh*qtq/(1.+sqrt(1.-qtq))
          tmax=.5*sh
          if(tmin.gt.tmax.and.ish.ge.1)write(ifmt,*)'psds:tmin,tmax'
     &                                              ,tmin,tmax

          ft=0.
          do i1=1,7
          do m1=1,2
            t=.5*(tmin+tmax+(2*m1-3)*x1(i1)*(tmin-tmax))
            u=sh-t
            qt=t*u/sh*(1.-qq/sh)
            if(qt.lt.qtmin.and.ish.ge.1)write(ifmt,*)'psds:qt,qtmin'
     &                                               ,qt,qtmin
            
            fb=psdsj(q1,xxe,sh,qt,t,u,qq,j,long)
            ft=ft+a1(i1)*fb*pssalf(qt/qcdlam)
          enddo
          enddo
          ft=ft*(tmax-tmin)
          fx1=fx1+a1(i)*ft*(1.-xx)/sh**2
        enddo
        enddo
        fx1=fx1*log((1.d0-xmin1)/(1.d0-xmax))
      endif

      if(xmin.lt..9d0)then
        xmax1=min(xmax,.9d0)
        do i=1,7
        do m=1,2
          xxe=xmin*(xmax1/xmin)**(.5-x1(i)*(m-1.5))
          xx=xxe

          sh=xx*s
          qtmin=max(1.d0*max(q2mass,q1),q2ini/(1.d0-xxe))
          qtq=4.*qtmin/(sh-qq)

          tmin=.5*sh*qtq/(1.+sqrt(1.-qtq))
          tmax=.5*sh
          if(tmin.gt.tmax.and.ish.ge.1)write(ifmt,*)'psds:tmin,tmax'
     *                                              ,tmin,tmax

          ft=0.
          do i1=1,7
          do m1=1,2
            t=(.5*(tmin+tmax+(2*m1-3)*x1(i1)*
     *      (tmin-tmax)))
            u=sh-t
            qt=t*u/sh*(1.-qq/sh)
            if(qt.lt.qtmin.and.ish.ge.1)write(ifmt,*)'psds:qt,qtmin'
     *                                               ,qt,qtmin
            
            fb=psdsj(q1,xxe,sh,qt,t,u,qq,j,long)
            ft=ft+a1(i1)*fb*pssalf(qt/qcdlam)
          enddo
          enddo
          ft=ft*(tmax-tmin)
          fx2=fx2+a1(i)*ft*xx/sh**2
        enddo
        enddo
        fx2=fx2*log(xmax1/xmin)
      endif
      psds=(fx1+fx2)*pi**2*alfe
      return
      end

c-----------------------------------------------------------------------
      function psdsj(q1,xx,s,qt,t,u,qq,j,long)
c-----------------------------------------------------------------------
c psdsj - integrand for dis singlet cross-section
c q1 - virtuality cutoff for current end of the ladder
c xx - lc momentum ratio between initial (j) and final (l) partons
c s  - c.m. energy squared for the born scattering,
c t  - invariant variable for the born scattering |(p1-p3)**2|
c u  - invariant variable for the born scattering |(p1-p4)**2|
c qq - photon virtuality
c j  - initial parton at the end of the ladder (0 - g, 1,2 - q)
c-----------------------------------------------------------------------
      double precision xx
      include 'epos.incsem'

      fb=psdbom(s,t,u,qq,long)
      psdsj=psevi(q1,qt,xx,min(1,iabs(j))+1,1)*fb
      return
      end

c------------------------------------------------------------------------
      function psdsin(q1,qq,s,m1,long)
c-----------------------------------------------------------------------
c psdsin - DIS singlet cross-section interpolation
c q1 - effective momentum cutoff for current end of the ladder,
c qq - photon virtuality,
c s -  total c.m. energy squared for the ladder,
c m1 - parton type at current end of the ladder (0 - g, 1,2 - q)
c-----------------------------------------------------------------------
      double precision xx
      dimension wi(3),wj(3),wk(3)
      common /psar2/  edmax,epmax
      common /psar25/ csdsi(21,21,104)
      include 'epos.incsem'
      include 'epos.inc'

      psdsin=0.
      m=min(1,iabs(m1))+1

      q2mass=qcmass**2
      s2min=4.*max(q2min,q2mass)+qq
      sdmin=4.*max(q1,q2mass)+qq
      if(m1.ne.0)then
        s2min=s2min/(1.-4.*q2ini/(s2min-qq))
        sdmin=sdmin/(1.-4.*q2ini/(sdmin-qq))
      endif
c      if(s.le.1.e8*sdmin)goto 2  !????????????????
      if(s.le.sdmin)goto 2
      
      qmin=q2min
      qmax=(s-qq)/4.
      if(m1.ne.0)qmax=(s-qq+sqrt((s-qq)**2-16.*s*q2ini))/8.
      qtq=4.*max(q2mass,q1)/(s-qq)
      if(qtq.lt.1.)then
        tmin=.5*s*qtq/(1.+sqrt(1.-qtq))
      else
        tmin=.5*s
      endif
      tmax=s/2.

      qlj=log(qq/q2min)*2.+1.
      j=int(qlj)
      if(j.lt.1)j=1
      if(j.gt.19)j=19
      wj(2)=qlj-j
      wj(3)=wj(2)*(wj(2)-1.)*.5
      wj(1)=1.-wj(2)+wj(3)
      wj(2)=wj(2)-2.*wj(3)

      qli=log(q1/qmin)/log(qmax/qmin)*20.+1.
      i=int(qli)
      if(i.lt.1)i=1
      if(i.gt.19)i=19
      wi(2)=qli-i
      wi(3)=wi(2)*(wi(2)-1.)*.5
      wi(1)=1.-wi(2)+wi(3)
      wi(2)=wi(2)-2.*wi(3)

      sl=log(s/s2min)/log(edmax/s2min)*25.+1.
      k=int(sl)
      if(k.lt.1)k=1
      if(k.gt.24)k=24
      wk(2)=sl-k
      wk(3)=wk(2)*(wk(2)-1.)*.5
      wk(1)=1.-wk(2)+wk(3)
      wk(2)=wk(2)-2.*wk(3)
      
      dsin1=0.
      do k1=1,3
        k2=k+k1-1+26*(m-1)+52*long
      do i1=1,3
      do j1=1,3
        dsin1=dsin1+csdsi(i+i1-1,j+j1-1,k2)*wi(i1)*wj(j1)*wk(k1)
      enddo
      enddo
      enddo
      if(m1.eq.0)then
        psdsin=psdsin+exp(dsin1)*(1./tmin-1./tmax)
      else
        psdsin=psdsin+max(0.,dsin1)
      endif

2     continue
      if(long.eq.0.and.q1.lt.qq.and.s.gt.qq/(1.-q2ini/qq))then
        xx=qq/s
        dsi=psevi(q1,qq,xx,m,2)*4.*pi**2*alfe/s
        if(m1.eq.0)then
          psdsin=psdsin+max(dsi,0.)
        else
          dnsi=psevi(q1,qq,xx,3,2)*4.*pi**2*alfe/s
          psdsin=psdsin+max(dsi-dnsi,0.)
        endif
      endif
      return
      end














c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################
c
c                        unused 3P
c
c###########################################################################
c###########################################################################
c###########################################################################
c###########################################################################


c------------------------------------------------------------------------
      function psvy(xpp0,xpr0,xpm0,xmr0,b,iqq)
c-----------------------------------------------------------------------
c psvy - 3p-contributions to the interaction eikonal
c xpp - lc+ for the pomeron,
c xpr - lc+ for the remnant,
c xpm - lc- for the pomeron,
c xpr - lc- for the remnant,
c b   - impact parameter,
c iqq=1  - Y-proj-uncut
c iqq=2  - Y-proj-1-cut
c iqq=3  - Y-proj-2-cut
c iqq=4  - Y-proj-soft-cut
c iqq=5  - Y-proj-gss-cut
c iqq=6  - Y-proj-qss-cut
c iqq=7  - Y-proj-ssg-cut
c iqq=8  - Y-proj-ssq-cut
c iqq=9  - Y-proj-difr
c iqq=-1 - Y-targ-uncut
c iqq=-2 - Y-targ-1-cut
c iqq=-3 - Y-targ-2-cut
c iqq=-4 - Y-targ-soft-cut
c iqq=-5 - Y-targ-gss-cut
c iqq=-6 - Y-targ-qss-cut
c iqq=-7 - Y-targ-ssg-cut
c iqq=-8 - Y-targ-ssq-cut
c iqq=-9 - Y-targ-difr
c------------------------------------------------------------------------
      
      psvy=0. 
      return 
      end
             
c------------------------------------------------------------------------
      function psvx(xpp,xpr,xpm,xmr,b,iqq)
c-----------------------------------------------------------------------
c psvx - 4p-contributions to the interaction eikonal
c xpp - lc+ for the pomeron,
c xpr - lc+ for the remnant,
c xpm - lc- for the pomeron,
c xpr - lc- for the remnant,
c b   - impact parameter,
c iqq=0  - X-uncut
c iqq=1  - X-1-cut
c iqq=2  - X-Y+cut
c iqq=-2 - X-Y-cut
c iqq=3  - X-1-cut-soft
c iqq=4  - X-1-cut-gss
c iqq=-4 - X-1-cut-gss
c iqq=5  - X-1-cut-qss
c iqq=-5 - X-1-cut-qss
c iqq=6  - X-difr+
c iqq=-6 - X-difr-
c------------------------------------------------------------------------
      
      psvx=0. 
      return
      end
           




c
c
c
c
c
cc------------------------------------------------------------------------
c      function psftig(x,xpomr,jj)
cc-----------------------------------------------------------------------
c
c        psftig=0.
c      end
c      
c         
cc------------------------------------------------------------------------
c      function psftih(zz,ddd)
cc-----------------------------------------------------------------------
c
c      psftih=0.
c      return
c      end          
c     
cc------------------------------------------------------------------------
c      function psftij(zz,ddd)
cc------------------------------------------------------------------------
c
c      psftij=0.
c      return
c      end
c         
cc------------------------------------------------------------------------
c      function psftik(xp,del1,rh1,rh2,rh3,rp,z)
cc-----------------------------------------------------------------------
c
c      psftik=0.
c      return
c      end          
c     
cc------------------------------------------------------------------------
c      function psftim(xp1,xp,del1,rh1,rh2,rh3,rp,z)
cc-----------------------------------------------------------------------
c
c      psftim=0.
c      return
c      end          
c     
cc------------------------------------------------------------------------
c      function psftil(xp,del1,rh1,rh2,rh3,rp,z)
cc------------------------------------------------------------------------
c
c      psftil=0.
c      return
c      end
c
cc------------------------------------------------------------------------
c      function psftigt(x)
cc-----------------------------------------------------------------------
c
c      psftigt=0.
c       return
c       end
c                  
cc------------------------------------------------------------------------
c      function psftist(x)
cc-----------------------------------------------------------------------
c      psftist=0.
c       return
c       end
c                  
cc------------------------------------------------------------------------
c      function psftig1(x,xpomr,jj)
cc-----------------------------------------------------------------------
c      psftig1=0.
c      return
c      end
c           
cc------------------------------------------------------------------------
c      function psftis1(x,xpomr,jj)
cc-----------------------------------------------------------------------
c      psftis1=0.
c      return
c      end
c
cc------------------------------------------------------------------------
c      function psd3p1(xd,xpomr,qq,jj)
cc-----------------------------------------------------------------------
cc psd3p1 - df2difr/dx_pomr
cc xd - bjorken x,
cc xpomr - pomeron x,
cc qq - photon virtuality
cc jj=1 - 1st order
cc-----------------------------------------------------------------------
c
c      psd3p1=0.
c
c      return
c      end
c      
cc------------------------------------------------------------------------
c      function psv3p(sy,xpp,xpm,zb)
cc-----------------------------------------------------------------------
cc psv3p - 3p-contributions to the interaction eikonal
cc sy  - energy squared for the hard interaction,
cc xpp - lc+ for the sh pomeron,
cc xpm - lc- for the sh pomeron,
cc z   - impact parameter factor, z=exp(-b**2/4*rp),
cc------------------------------------------------------------------------
c      
c      psv3p=0.      
c      return 
c      end
c                  
cc------------------------------------------------------------------------
c      function psfro(xpomr,zb,iclp,icdpro)
cc-----------------------------------------------------------------------
cc psfro - generalized froissaron between proj. and 3p-vertex
cc xpp   - lc+ for the proj. side,
cc xpomr - lc+ for the vertex,
cc zb    - impact parameter factor, z=exp(-b**2/4*rp),
cc------------------------------------------------------------------------
c      psfro=0.
c      return
c      end
c                       
cc------------------------------------------------------------------------
c      function psv2(xpomr,zb,iclp,iqq)
cc-----------------------------------------------------------------------
cc psv2 - 2-pom contribution to the froissaron 
cc xpomr - lc+ for the vertex,
cc zb    - impact parameter factor, z=exp(-b**2/4*rp),
cc------------------------------------------------------------------------
c      
c      psv2=0.
c      return
c      end
c           
cc------------------------------------------------------------------------
c      function psvfro(xpp,xpr,xpomr,zb,iclp,icdpro,iqq)
cc-----------------------------------------------------------------------
cc psvfro - effective froissaron contributions
cc xpomr - lc+ for the vertex,
cc zb    - impact parameter factor, z=exp(-b**2/4*rp),
cc iqq=0 - total uncut
cc iqq=1 - total 1-cut
cc iqq=2 - total 2-cut
cc iqq=3 - soft 1-cut
cc iqq=4 - gg 1-cut
cc iqq=5 - qg 1-cut
cc iqq=6 - difr
cc------------------------------------------------------------------------
c      
c      psvfro=0.
c      return
c      end
c          
cc------------------------------------------------------------------------
c      function psfroin(xpomr,z,iclp,icdpro)
cc-----------------------------------------------------------------------
cc psfroin - interpolation of effective froissaron corrections
cc xpomr - lc+ for the 3p-vertex,
cc z   - impact parameter factor, z=exp(-b**2/4*rp),
cc-----------------------------------------------------------------------
c      
c      psfroin=0.
c      return 
c      end
c          
cc------------------------------------------------------------------------
c      function psvnorm(b)
cc-----------------------------------------------------------------------
cc psvnorm - X-contribution normalization
cc b   - impact parameter
cc-----------------------------------------------------------------------
c      
c      psvnorm=0.
c      return 
c      end
c                    
cc------------------------------------------------------------------------
c      function psvxb(dxpp,dxpr,dxpm,dxmr,dxpomr,bb1,bb2
c     *,iclp,iclt,icdpro,icdtar,iqq)
cc-----------------------------------------------------------------------
cc psvxb   - integrand for X-contributions
cc dxpomr  - lc+ for the vertex,
cc bb1,bb2 - impact parameters to proj(targ),
cc iqq=0   - uncut
cc iqq=1   - 1-cut
cc iqq=2   - Y-cut
cc iqq=3   - soft 1-cut
cc iqq=4   - gg 1-cut
cc iqq=5   - qg 1-cut
cc iqq=6   - difr
cc------------------------------------------------------------------------
c      double precision dxpp,dxpr,dxpm,dxmr,dxpomr
c      
c      psvxb=0.
c      return
c      end
c      
cc------------------------------------------------------------------------
c      function pscoef(zz,alp1,alp2,iclp,iqq)
cc-----------------------------------------------------------------------
cc pscoef - integrated vertexes
cc zz=xpomr/xpr(xpp)
cc iqq=0 - 2-cut
cc iqq=1 - 1-uncut
cc iqq=2 - 2-uncut
cc------------------------------------------------------------------------
c      
c      pscoef=0.
c      return
c      end
c      
cc------------------------------------------------------------------------
c      function pscoefi(z,i1,i2,iclp,iqq)
cc-----------------------------------------------------------------------
cc pscoefi - interpolation of integrated vertexes
cc z=xpomr/xpr(xpp)
cc iqq=0 - 2-cut
cc iqq=1 - 1-uncut
cc iqq=2 - 2-uncut
cc------------------------------------------------------------------------
c      pscoefi=0.
c      return 
c      end
c       
c                    
