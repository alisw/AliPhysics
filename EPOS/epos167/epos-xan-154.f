c---------------------------------------------------------------------
      subroutine xiniall
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr)
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1
      parameter(mxxhis=50)
      common/varhis/icorrtrig(0:mxxhis),ihardevent(0:mxxhis)
     &,ijetfind1(0:mxxhis),ijetfind2(0:mxxhis)

      nhis=0
      nfra=0
      imulty1=0
      do n=0,mxxhis
        icorrtrig(n)=0
        ihardevent(n)=0
        ijetfind1(n)=0
        ijetfind2(n)=0
      enddo
      do n=1,mxhis
        xpara(1,n)=0
      enddo
      ncontrall=0
      noerrall=0
      
      end

c---------------------------------------------------------------------
      subroutine xini
c---------------------------------------------------------------------
c  called after beginhisto
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1

      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      character line*160,cvar*6
      logical go
      common/nl/noplin
      character*160 cline
      common/cjjj/jjj,cline
      
      call utpri('xini  ',ish,ishini,5)

      i=1
                                !      iapl=0
                                !      nhis=0
      j=jjj     !-1
      line=cline
                                !      nfra=1
                                !      ifra(1)=iframe
      iapl=0
      if(nfra.eq.0)then     
        nfra=1
        ifra(1)=iframe
      endif
      nhis=nhis+1
      if(nhis.gt.mxhis)stop'xini: mxhis too small.       '
      noweak(nhis)=0
      ionoerr=0
c      newfra=0
      indfra=1
c      nepfra=0
      inpfra=1
 1    call utword(line,i,j,0)
      if(line(i:j).eq.'application')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'analysis')then
                                !iapl=0
                                !nhis=nhis+1
                                !newfra=0
                                !indfra=1
                                !nepfra=0
                                !inpfra=1
        else
          iapl=1
        endif
      elseif(line(i:j).eq.'input')then !-----------
        call utword(line,i,j,0)
        if(nopen.ge.0)then
         nopen=nopen+1
         if(nopen.gt.9)stop'too many nested input commands'
         open(unit=20+nopen,file=line(i:j),status='old')
         if(iprmpt.eq.1)iprmpt=-1
        endif 
      elseif(line(i:j).eq.'runprogram')then !-----------
        if(iapl.eq.0)then
        else
          goto 9999
        endif
      elseif(line(i:j).eq.'frame'.or.line(i:j).eq.'frame+')then !------
        ifp=1
        if(line(i:j).eq.'frame+')ifp=2
        call utword(line,i,j,1)
        if(line(i:j).eq.'total')then
          nfp=iframe
        elseif(line(i:j).eq.'nucleon-nucleon')then
          nfp=11
        elseif(line(i:j).eq.'target')then
          nfp=12
        elseif(line(i:j).eq.'gamma-nucleon')then
          nfp=21
        elseif(line(i:j).eq.'lab')then
          nfp=22
        elseif(line(i:j).eq.'breit')then
          nfp=23
        elseif(line(i:j).eq.'thrust')then
          nfp=33
        elseif(line(i:j).eq.'sphericity')then
          nfp=32
        endif
        go=.true.
        do l=1,nfra
          if(ifra(l).eq.nfp)then
            inl=l
            go=.false.
          endif
        enddo
        if (go) then
          nfra=nfra+1
          inl=nfra
          ifra(nfra)=nfp
        endif
        if(ifp.eq.1)then
          indfra=inl
c          newfra=nfp
          ivfra(1,nhis)=indfra
          ivfra(2,nhis)=indfra
        else
          inpfra=inl
c          nepfra=nfp
        endif
      elseif(line(i:j).eq.'binning')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'lin')then
          iologb=0
          iocnxb=0
        elseif(line(i:j).eq.'log')then
          iologb=1
          iocnxb=0
        elseif(line(i:j).eq.'clin')then
          iologb=0
          iocnxb=1
        elseif(line(i:j).eq.'clog')then
          iologb=1
          iocnxb=1
        else
          print *, 'what the heck is ',line(i:j),' binning?'
          print *, 'I will use the linear (lin) one'
        endif
      elseif(line(i:j).eq.'setm')then !-----------
        if(iapl.eq.0) then
          print *,"You should use histogram instead of setm, please"
          stop
        endif
      elseif(line(i:j).eq.'set')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'iologb')then
          call utword(line,i,j,1)
          read(line(i:j),*) iologb
        elseif(line(i:j).eq.'iocnxb')then
          call utword(line,i,j,1)
          read(line(i:j),*) iocnxb
        elseif(line(i:j).eq.'etacut')then
          call utword(line,i,j,1)
          read(line(i:j),*) etacut
        elseif(line(i:j).eq.'nemsi')then
          call utword(line,i,j,1)
          read(line(i:j),*)nemsi
        endif
      elseif(line(i:j).eq.'xpara')then !-----------
        call utword(line,i,j,1) 
        read(line(i:j),*)ipara
        if(ipara.gt.mxpara)stop'mxpara too small.         '
        call utword(line,i,j,1) 
        read(line(i:j),*)val
        xpara(ipara,nhis)=val
      elseif(line(i:j).eq.'echo')then !-----------
        call utword(line,i,j,1)
        if(line(i:j).eq.'on')iecho=1
        if(line(i:j).eq.'off')iecho=0
        if(line(i:j).ne.'on'.and.line(i:j).ne.'off')
     *  stop'invalid option'
      elseif(line(i:j).eq.'noweak')then !-----------
        noweak(nhis)=1
      elseif(line(i:j).eq.'histogram')then !-----------
        nac(nhis)=1
        call utword(line,i,j,1) !xvaria
        cvar='      '
        cvar=line(i:j)
        call xtrans(cvar,inom,ifrnew,nhis)
        if(inom.eq.-1)then
          if(line(i:i).ge.'0'.and.line(i:i).le.'9')then
            inom=298
            read(line(i:j),*) sval(1,nhis)
          endif
        endif
        ivar(1,nhis)=inom
        if(ifrnew.ne.0)then     !check frame for e+e- event 
          go=.true.             !shape variables
          do l=1,nfra
            if(ifra(l).eq.ifrnew)then
              indfra=l
              go=.false.        !have it already
            endif
          enddo
          if (go) then
            nfra=nfra+1
            ifra(nfra)=ifrnew
            indfra=nfra
          endif
        endif
        call utword(line,i,j,1) !yvaria
        cvar='      '
        cvar=line(i:j)
        call xtrans(cvar,inom,ifrnew,nhis)
        ivar(2,nhis)=inom
        if(inom.eq.-1)then
          if(line(i:i).ge.'0'.and.line(i:i).le.'9')then
            inom=299
            read(line(i:j),*) sval(2,nhis)
          endif
        endif
        if(inom.eq.-1)ivar(1,nhis)=inom

        ivfra(1,nhis)=indfra
        ivfra(2,nhis)=indfra

        call utword(line,i,j,1) !normation
        read(line(i:j),*) inorm(nhis)

        call utword(line,i,j,1) !xmin
        if(line(i:j).eq.'egy')then
         if(engy.gt.0)then
          egy=engy
         elseif(ecms.gt.0.)then
          egy=ecms   
         elseif(elab.gt.0)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*elab*atg+atg**2+apj**2 )
         elseif(ekin.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*(ekin+apj)*atg+atg**2+apj**2 )
         elseif(pnll.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*sqrt(pnll**2+apj**2)*atg+atg**2+apj**2 )
         else
          stop'pb in xini (1).   '
         endif
         xmin(nhis)=egy-0.001
        else
         read(line(i:j),*) xmin(nhis)
        endif
        
        call utword(line,i,j,1) !xmax
        if(line(i:j).eq.'egy')then
         if(engy.gt.0)then
          egy=engy
         elseif(ecms.gt.0.)then
          egy=ecms   
         elseif(elab.gt.0)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*elab*atg+atg**2+apj**2 )
         elseif(ekin.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*(ekin+apj)*atg+atg**2+apj**2 )
         elseif(pnll.gt.0.)then
          call idmass(idproj,apj)
          call idmass(idtarg,atg)
          egy=sqrt( 2*sqrt(pnll**2+apj**2)*atg+atg**2+apj**2 )
         else
          stop'pb in xini (2).   '
         endif
         xmax(nhis)=egy+0.001
        else
         read(line(i:j),*) xmax(nhis)
        endif

        call utword(line,i,j,1) !nbin
        read(line(i:j),*) nbin(nhis)
        do l=1,nbin(nhis)
          bin(l,nac(nhis),nhis)=0.
          zcbin(l,nac(nhis),nhis)=0
        enddo
        lookcontr(nhis)=0
        lookcontrx(nhis)=0
        inoerr(nhis)=0
      elseif(line(i:j).eq.'idcode')then !-----------
        call utword(line,i,j,1) !idcode
        if(line(i:i+2).eq.'995')stop'xini: idcode 995 not supported'
        if(line(i:i+2).eq.'994')stop'xini: idcode 994 not supported'
        nidcod(nhis)=nidcod(nhis)+1
        read(line(i:j),*) idcod(nidcod(nhis),nhis)
        idmod(nidcod(nhis),nhis)=.false.
      elseif(line(i:j).eq.'idcode+')then !-----------
        stop'xini: idcode+ not supported'
        call utword(line,i,j,1) !idcode
        if(line(i:i+2).eq.'995')stop'xini: idcode 995 not supported'
        if(line(i:i+2).eq.'994')stop'xini: idcode 994 not supported'
        nidcod(nhis)=nidcod(nhis)+1
        read(line(i:j),*) idcod(nidcod(nhis),nhis)
        idmod(nidcod(nhis),nhis)=.true.
      elseif(line(i:j).eq.'trigger')then !-----------
        call utword(line,i,j,1) 
        ntc=1
        imo=1
        ncontr=0
        icontrtyp(nhis)=0
        if(line(i:j).eq.'or'.or.line(i:j).eq.'contr')then
          imo=2
          if(line(i:j).eq.'contr')imo=3
          call utword(line,i,j,1) 
          read(line(i:j),*)ztc
          ntc=nint(ztc)
          call utword(line,i,j,1) 
          if(imo.eq.3)then
            ncontr=ntc
            ncontrall=ncontrall+ncontr
            if(ncontrall.gt.mxcontr)stop'xini: mxcontr too small.     '
            if(ncontr.gt.mxcnt)stop'xini: mxcnt too small.     '
            lookcontr(nhis)=ncontrall-ncontr+1
            lookcontrx(nhis)=ncontrall
            do nb=1,nbin(nhis)
              do nn=1,ncontr
                bbin(nb,nac(nhis),lookcontr(nhis)-1+nn)=0.d0
                zbbin(nb,nac(nhis),lookcontr(nhis)-1+nn)=0.d0
              enddo
            enddo
            do nn=1,ncontr
                    nccevt(lookcontr(nhis)-1+nn)=0
            enddo
          endif  
        endif  
        do n=1,ntc
          if(n.ne.1)call utword(line,i,j,1) !trigger-name
          cvar='      '
          ifp=1
          if(line(j:j).eq.'+')then
            cvar=line(i:j-1)
            ifp=2
          else
            cvar=line(i:j)
            ifp=1
          endif
          call xtrans(cvar,inom,ifrnew,nhis)
          if(inom.gt.0)then
            ntri(nhis)=ntri(nhis)+1
            if(ntc.eq.1)then
              ntrc(ntri(nhis),nhis)=1
            elseif(n.eq.1)then
              ntrc(ntri(nhis),nhis)=2
            elseif(n.eq.ntc)then
              ntrc(ntri(nhis),nhis)=3
            else  
              ntrc(ntri(nhis),nhis)=0
            endif
            if(imo.eq.3)then
              ntrc(ntri(nhis),nhis)=-1
              if(n.eq.1)then
                icontrtyp(nhis)=1+inom/100
              else        
                if(1+inom/100.ne.icontrtyp(nhis))
     *               stop'xini: type mismatch'                
              endif
            endif  
            itri(ntri(nhis),nhis)=inom
            if(ifp.eq.1)then
              itfra(ntri(nhis),nhis)=indfra
            else
              itfra(ntri(nhis),nhis)=inpfra
            endif
            xmitrp(ntri(nhis),nhis)=100.
            xmatrp(ntri(nhis),nhis)=100.
            call utword(line,i,j,1) !-----------xmin----------
            if(line(i:j).eq.'inf')then
              xmitri(ntri(nhis),nhis)=1e30
            elseif(line(i:j).eq.'-inf')then
              xmitri(ntri(nhis),nhis)=-1e30
            elseif(line(i:j).eq.'A')then
              xmitri(ntri(nhis),nhis)=maproj
            elseif(line(i:j).eq.'A+1')then
              xmitri(ntri(nhis),nhis)=maproj+1
            elseif(line(i:j).eq.'A+B')then
              xmitri(ntri(nhis),nhis)=maproj+matarg
            elseif(line(i:j).eq.'A+B+1')then
              xmitri(ntri(nhis),nhis)=maproj+matarg+1
            elseif(line(i:j).eq.'lead')then    !leading particle (neads Standard Variable)
              xmitri(ntri(nhis),nhis)=-123456
              imulty1=1 
            else
              kk=0
              do k=i+1,j-1
                if(line(k:k).eq.'%')kk=k
              enddo
              if(kk.eq.0)then
                read(line(i:j),*)xmitri(ntri(nhis),nhis)
              else
                read(line(i:kk-1),*)xmitrp(ntri(nhis),nhis)
                read(line(kk+1:j),*)xmitri(ntri(nhis),nhis) 
              endif        
            endif
            call utword(line,i,j,1) !-----------xmax------------
            if(line(i:j).eq.'inf')then
              xmatri(ntri(nhis),nhis)=1e30
            elseif(line(i:j).eq.'-inf')then
              xmatri(ntri(nhis),nhis)=-1e30
            elseif(line(i:j).eq.'A')then
              xmatri(ntri(nhis),nhis)=maproj
            elseif(line(i:j).eq.'A+1')then
              xmatri(ntri(nhis),nhis)=maproj+1
            elseif(line(i:j).eq.'A+B')then
              xmatri(ntri(nhis),nhis)=maproj+matarg
            elseif(line(i:j).eq.'A+B+1')then
              xmatri(ntri(nhis),nhis)=maproj+matarg+1
            elseif(line(i:j).eq.'lead')then    !leading particle (neads Standard Variable)
              xmatri(ntri(nhis),nhis)=-123456
              imulty1=1 
            else
              kk=0
              do k=i+1,j-1
                if(line(k:k).eq.'%')kk=k
              enddo
              if(kk.eq.0)then
                read(line(i:j),*)xmatri(ntri(nhis),nhis)
                xmatrp(ntri(nhis),nhis)=100.
              else
                read(line(i:kk-1),*)xmatrp(ntri(nhis),nhis)
                read(line(kk+1:j),*)xmatri(ntri(nhis),nhis)
              endif        
            endif
            !---exchange min-max------------------
            if(xmitri(ntri(nhis),nhis).gt.xmatri(ntri(nhis),nhis))then
              xmatri_save=xmatri(ntri(nhis),nhis)
              xmatrp_save=xmatrp(ntri(nhis),nhis)
              xmatri(ntri(nhis),nhis)=xmitri(ntri(nhis),nhis)
              xmatrp(ntri(nhis),nhis)=xmitrp(ntri(nhis),nhis)
              xmitri(ntri(nhis),nhis)=xmatri_save
              xmitrp(ntri(nhis),nhis)=xmatrp_save
            endif
            !-------------------------------------
          else
            ivar(1,nhis)=-1
            call utword(line,i,j,1) !xmin
            call utword(line,i,j,1) !xmax
          endif
        enddo
      elseif(line(i:j).eq.'noerrorbut')then !-----------
        ionoerr=ionoerr+1
        if(ionoerr.gt.2)stop'xini: not more than 2 noerrorbut !   '
        noerrall=noerrall+1
        if(noerrall.gt.mxhis/2)stop'xini: to many noerrorbut     '
        
        call utword(line,i,j,1) !variable-name
        cvar=line(i:j)
        call xtrans(cvar,inom,ifrnew,nhis)
        if(inom.gt.0)then
          if(inom.gt.100)then
            write(*,*)'xini: noerrorbut can not be used with :',cvar
            stop'xini: error with noerrorbut!'
          endif
          noerrhis(nhis)=noerrall-ionoerr+1
          noerr(noerrhis(nhis),ionoerr)=inom
          do nb=1,nbin(nhis)
             ebin(nb,nac(nhis),ionoerr-1+noerrhis(nhis))=0.d0
                  zebin(nb,nac(nhis),ionoerr-1+noerrhis(nhis))=0.d0
          enddo
        else
          ionoerr=ionoerr-1
          noerrall=noerrall-1
        endif
        inoerr(nhis)=ionoerr
      elseif(line(i:j).eq.'write')then !-----------
        call utword(line,i,j,1)
      elseif(line(i:j).eq.'writearray')then !-----------
        call utword(line,i,j,1)
        iologb=0
        iocnxb=0
      elseif(line(i:j).eq.'writehisto')then !-----------
        call utword(line,i,j,1)
        iologb=0
        iocnxb=0
      elseif(line(i:j).eq.'endhisto')then   !-----------
        ilog(nhis)=.false.
        icnx(nhis)=.false.
        if(iologb.eq.1)ilog(nhis)=.true.
        if(iocnxb.eq.1)icnx(nhis)=.true.
        if(ilog(nhis))then
          xinc(nhis)=1./log(xmax(nhis)/xmin(nhis))*nbin(nhis)
        else
          xinc(nhis)=float(nbin(nhis))/(xmax(nhis)-xmin(nhis))
        endif
        iologb=0
        iocnxb=0
        jjj=j
        cline=line
        goto 9999
      endif
      goto 1
      
 9999 continue 
      if(ish.ge.5)then
        do n=1,nhis
          write (ifch,*) ivar(1,n),ivar(2,n),'(',ivfra(1,n),ivfra(2,n)
     $         ,')',inorm(n)
     $         ,xmin(n),xmax(n),ilog(n),icnx(n)
     $         ,nbin(n),(idcod(j,n),j=1,nidcod(n))
     $         ,' tri:',ntri(n),(itri(j,n),j=1,ntri(n)),'('
     $         ,(itfra(j,n),j=1,ntri(n)),')'
     $         ,(xmitri(j,n),j=1,ntri(n)) ,(xmatri(j,n),j=1,ntri(n))
        enddo
        write (ifch,*) (ifra(j),j=1,nfra)
      endif
      call utprix('xini  ',ish,ishini,5)
      return
      end


c---------------------------------------------------------------------
      subroutine xana
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr)
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
     
      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
      double precision bofra
      common/dfra/bofra(5,mxfra)
      parameter (ntim=1000)
      common/cprt/nprtj,pprt(5,ntim),idprt(ntim),iorprt(ntim)
     &     ,idaprt(2,ntim)

      double precision pgampr,rgampr
      common/cgampr/pgampr(5),rgampr(4)

      dimension ten(4,3)
      logical go,goo(mxcnt)

      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1
      parameter(mxxhis=50)
      common/varhis/icorrtrig(0:mxxhis),ihardevent(0:mxxhis)
     &,ijetfind1(0:mxxhis),ijetfind2(0:mxxhis)
     
      call utpri('xana  ',ish,ishini,4)  
      
      do n=1,nhis
      do i=1,mypara
        ypara(i,n)=0
      enddo        
      enddo
          

      if(ish.ge.5)write(ifch,*)'frames ...'      

      if(iappl.eq.6)then
        if(mod(iolept/10,10).eq.1) call gakjet(1)
        if(mod(iolept/100,10).eq.1) call gakjet(2)
      endif

      do l=1,nfra
        emax(l)=egyevt/2
        if(ifra(l).eq.12)emax(l)=sqrt(pnll**2+prom**2)
        if(ifra(l).eq.iframe)then
          imofra(1,l)=0
          imofra(2,l)=0
          imofra(3,l)=0
        elseif(ifra(l).eq.11.or.ifra(l).eq.12)then
          imofra(1,l)=0
          imofra(2,l)=0
          bofra(1,l)=0d0
          bofra(2,l)=0d0
          bofra(3,l)=dsinh(dble(yhaha))
          bofra(4,l)=dcosh(dble(yhaha))
          bofra(5,l)=1d0
          if(ifra(l).eq.11.and.iframe.eq.12)then
            imofra(3,l)=1       ! target -> NN 
          elseif(ifra(l).eq.12.and.iframe.eq.11)then
            imofra(3,l)=-1      ! NN -> target
          else
            imofra(3,l)=0       ! not known
          endif
        elseif(ifra(l).eq.21)then
          if(iframe.ne.21)then
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon for event run'
            stop'bye bye'
          endif
        elseif(ifra(l).eq.22)then
          if(iframe.eq.21)then                                
            imofra(1,l)=-1      !'  trafo gN -> lab'
            r1fra(1,l)=rgampr(1)
            r1fra(2,l)=rgampr(2)
            r1fra(3,l)=rgampr(3)
            imofra(2,l)=0
            imofra(3,l)=-1    
            bofra(1,l)=pgampr(1)
            bofra(2,l)=pgampr(2)
            bofra(3,l)=pgampr(3)
            bofra(4,l)=pgampr(4)
            bofra(5,l)=pgampr(5)
          elseif(iframe.eq.22)then
                                ! nothing to do already gN-frame
          else
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon or lab for event run'
            stop'bye bye'
          endif  
        elseif(ifra(l).eq.23)then
          if(iframe.eq.21)then                                
            imofra(1,l)=0       ! gN -> breit-frame
            r1fra(1,l)=rgampr(1)
            r1fra(2,l)=rgampr(2)
            r1fra(3,l)=rgampr(3)
            imofra(2,l)=0
            imofra(3,l)=1    
            bofra(1,l)=0d0
            bofra(2,l)=0d0
            bofra(3,l)=rgampr(4)
            bofra(4,l)=sqrt(rgampr(1)**2+rgampr(2)**2+rgampr(3)**2)
            bofra(5,l)=sqrt( bofra(4,l)**2-rgampr(4)**2)
          elseif(iframe.eq.23)then
                                ! nothing to do already breit-frame
          else
            print *, 'invalid frame request'
            print *, 'choose frame gamma-nucleon or lab for event run'
            stop'bye bye'
          endif
        elseif(ifra(l).eq.33.or.ifra(l).eq.36)then
          if(ifra(l).eq.33)then
            call gakthru(ten,2)
          else
            call gakthru(ten,3)
          endif
          if(ten(4,1).lt.0.)then
            imofra(1,l)=0
            imofra(2,l)=0
            imofra(3,l)=0
          else
            arox=ten(1,1)
            aroy=ten(2,1)
            aroz=ten(3,1)
            brox=ten(1,2)
            broy=ten(2,2)
            broz=ten(3,2)
            call utrota(1,arox,aroy,aroz,brox,broy,broz)
            imofra(1,l)=1
            r1fra(1,l)=arox
            r1fra(2,l)=aroy
            r1fra(3,l)=aroz
            imofra(2,l)=1
            r2fra(1,l)=brox
            r2fra(2,l)=broy
            r2fra(3,l)=broz
            imofra(3,l)=0       !no boost
          endif
          bofra(1,l)=dble(ten(4,1)) !usually this is for boosting
          bofra(2,l)=dble(ten(4,2)) !I abuse it to store the eigenvalues
          bofra(3,l)=dble(ten(4,3)) !
        elseif(ifra(l).eq.32.or.ifra(l).eq.34.or.ifra(l).eq.35)then
          if(ifra(l).eq.32)then
            call gaksphe(ten,2.,2)
          elseif(ifra(l).eq.34)then
            call gaksphe(ten,1.,2)
          else
            call gaksphe(ten,2.,3)
          endif
          if(ten(4,1).lt.0.)then
            imofra(1,l)=0
            imofra(2,l)=0
            imofra(3,l)=0
          else
            arox=ten(1,1)
            aroy=ten(2,1)
            aroz=ten(3,1)
            brox=ten(1,2)
            broy=ten(2,2)
            broz=ten(3,2)
            call utrota(1,arox,aroy,aroz,brox,broy,broz)
            imofra(1,l)=1
            r1fra(1,l)=arox
            r1fra(2,l)=aroy
            r1fra(3,l)=aroz
            imofra(2,l)=1
            r2fra(1,l)=brox
            r2fra(2,l)=broy
            r2fra(3,l)=broz
            imofra(3,l)=0         
          endif
          bofra(1,l)=dble(ten(4,1))
          bofra(2,l)=dble(ten(4,2))
          bofra(3,l)=dble(ten(4,3))
        endif
      enddo
      
      do n=1,nhis
        itrevt(n)=.false.
        if(ivar(1,n).ge.100.and.ivar(1,n).le.199) sval(1,n)=0.
        if(ivar(2,n).ge.100.and.ivar(2,n).le.199) sval(2,n)=0.
        if(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
          call xval(n,ivar(1,n),ivfra(1,n),0,x) !initializing of  variables
        endif
        if(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
          call xval(n,ivar(2,n),ivfra(2,n),0,y) !
        endif
        do j=1,ntri(n)
          valtri(j,n)=0.
        enddo
        do j=1,nbin(n)  !copy bins
          bin(j,3-nac(n),n)=bin(j,nac(n),n)
          zcbin(j,3-nac(n),n)=zcbin(j,nac(n),n)
        enddo
        if(lookcontr(n).gt.0)then
          do j=1,nbin(n)
            do loo=lookcontr(n),lookcontrx(n)
                    bbin(j,3-nac(n),loo)=bbin(j,nac(n),loo)
                    zbbin(j,3-nac(n),loo)=zbbin(j,nac(n),loo)
            enddo
          enddo
        endif
        if(inoerr(n).gt.0)then
          do j=1,nbin(n)
            do nn=1,inoerr(n)
              ebin(j,3-nac(n),nn-1+noerrhis(n))=ebin(j,nac(n),
     &                                          nn-1+noerrhis(n))
              zebin(j,3-nac(n),nn-1+noerrhis(n))=zebin(j,nac(n),
     &                                          nn-1+noerrhis(n))
            enddo
          enddo
        endif
      enddo

      if(imulty1.eq.1)then
        if(ish.ge.5)write(ifch,*)'Calculate standard variables ...'
        call StandardVariables
      endif
      if(ish.ge.5)write(ifch,*)'Call corrtrig ...'
      do n=1,icorrtrig(0)
        call corrtrig(icorrtrig(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call hardevent ...'
      do n=1,ihardevent(0)
        call hardevent(ihardevent(n))
      enddo
      if(ish.ge.5)write(ifch,*)'Call jetfind ...'
      do n=1,ijetfind1(0)        
        call jetfind(1,ijetfind1(n))
      enddo
      do n=1,ijetfind2(0)        
        call jetfind(2,ijetfind2(n))
      enddo
      
c...........................loop nptl...................................      
      if(ish.ge.5)write(ifch,*)'Loop nptl ...'
      do j=1,nptl
        if(iorptl(j).lt.0.or.istptl(j).gt.istmax)goto8
        if(ish.ge.5)write(ifch,*)'ptl :',j
        call idchrg(idptl(j),ch)
        do i=1,nfra
          iffra(i)=0            !flag if frame calculated or not
        enddo
        do n=1,nhis
          if(ivar(1,n).eq.-1.or.ivar(2,n).eq.-1)goto 9

c...........check ids
          go=nidcod(n).eq.0
          do i=1,nidcod(n)
            if(istptl(j).eq.0.and.idcod(i,n).eq.9990)then
              if(abs(idptl(j)).ge.100
     $         .and.abs(idptl(j)).lt.10000) go=.true.
            elseif(istptl(j).eq.0.and.idcod(i,n).eq.9970)then
              if(abs(ch).gt.0.1.and.abs(idptl(j)).ge.100
     $         .and.abs(idptl(j)).lt.10000) go=.true.
            elseif(istptl(j).eq.0.and.idcod(i,n).eq.-9960)then
              if(ch.lt.-0.1.and.abs(idptl(j)).ge.100
     $             .and.abs(idptl(j)).lt.10000)go=.true.
            elseif(istptl(j).eq.0.and.idcod(i,n).eq.9960)then
              if(ch.gt.0.1.and.abs(idptl(j)).ge.100
     $         .and.abs(idptl(j)).lt.10000)go=.true.
            elseif((istptl(j).le.1.or.istptl(j).ge.10)
     $            .and.idcod(i,n).eq.idptl(j))then
              go=.true.        
            endif
          enddo
          if(ish.ge.10)write(ifch,*)j,' id,ist',idptl(j),istptl(j),go

c...........check weak decay 
          if(go)then
            if(noweak(n).eq.1)then  !do not consider weak decay products
             if(iorptl(j).ne.0)then
              idora=abs( idptl(iorptl(j)) )
              if(  idora.eq.20   .or.idora.eq.2130
     &               .or.idora.eq.2230 .or.idora.eq.1130  
     &         .or.idora.eq.2330 .or.idora.eq.1330  
     &         .or.idora.eq.3331 )go=.false.
             ! print *, j,n, '   ', idptl(j),idora,go
             endif
            endif
          endif
          
c...........check triggers
          if(go)then
            if(ish.ge.7)write(ifch,*)'  check triggers in histogram ',n
            ncontr=0
            do i=1,ntri(n)
              if(ish.ge.7)write(ifch,*)'  trigger variable: ',itri(i,n)
              if(itri(i,n).lt.100)then
                call xval(n,itri(i,n),itfra(i,n),j,x)
                if(ntrc(i,n).ne.-1)then
                  call triggercondition(i,n,x,go)
                else
                  ncontr=ncontr+1 
                  goo(ncontr)=.true.
                  call triggercondition(i,n,x,goo(ncontr))
                  if((ivar(1,n).gt.100.and.ivar(1,n).lt.200)
     .                 .or.(ivar(2,n).gt.100.and.ivar(2,n).lt.200))then
                    print*,'!-----------------------------------------'
                    print*,'!  100-199 event variables can not be used'
                    print*,'! in connection with "trigger contr ..."  '
                    print*,'!-----------------------------------------'
                    stop'in xana (1).                      '
                  endif
                endif  
              elseif(itri(i,n).lt.200)then
                if(ntrc(i,n).eq.-1)then
                    print*,'!-----------------------------------------'
                    print*,'!  100-199 event variables can not be used'
                    print*,'! in connection with "trigger contr ..."  '
                    print*,'!-----------------------------------------'
                    stop'in xana (2).                      '
                endif
                call xval(n,itri(i,n),itfra(i,n),j,x)
                valtri(i,n)=valtri(i,n)+x
              endif
            enddo
          endif
                 
c............fill histogram 
          if(go)then
            if(ish.ge.7)write(ifch,*)'  fill histogram '
     &            ,n,ivar(1,n),ivar(2,n),ivfra(2,n)
            if(ivar(1,n).lt.100.or.ivar(2,n).lt.100)then
              if(ivar(2,n).lt.100)then
                call xval(n,ivar(2,n),ivfra(2,n),j,y)
                sval(2,n)=y
              endif
              if(ivar(1,n).lt.100)then
                call xval(n,ivar(1,n),ivfra(1,n),j,x)
                if(x.ge.xmin(n).and.x.le.xmax(n))then
                  norm3=mod(inorm(n)/100,10)
                  if(norm3.eq.1)then
                    y=y*x
                  elseif(norm3.eq.2.and.ivar(1,n).eq.63.and.x.ne.0.)then
                    y=y/(x+pptl(5,j))/2/pi
                  elseif(norm3.eq.2.and.ivar(1,n).ne.63.and.x.ne.0.)then
                    y=y/x/2/pi
                  elseif(norm3.eq.4.and.x.ne.0.)then
                    y=y/x**1.5
                  elseif(norm3.eq.5.and.x.ne.0.)then
                    y=y/x
                  elseif(norm3.eq.7.and.x.ne.0.)then
                    y=y/x/sqrt(x-pptl(5,j))
                  endif
                  if(icnx(n))then
                    call fillhistoconex(n,x,y,ivfra(2,n),j)   !for conex
                  else
                    if(ilog(n))then
                      nb=1+int(log(x/xmin(n))*xinc(n)) 
                    else
                      nb=1+int((x-xmin(n))*xinc(n))
                    endif
                    bin(nb,nac(n),n)=bin(nb,nac(n),n)+y
                    if(ncontr.gt.0)then  !ptl trigger contr
                      do nn=1,ncontr
                        if(goo(nn))then
                             bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &                  bbin(nb,nac(n),lookcontr(n)-1+nn)+y
                             zbbin(nb,nac(n),lookcontr(n)-1+nn)=
     &                  zbbin(nb,nac(n),lookcontr(n)-1+nn)+1
                        endif
                      enddo
                    endif
                    if(inoerr(n).gt.0)then
                      do nn=1,inoerr(n)
                       call xval(n,noerr(noerrhis(n),nn),ivfra(2,n),j,y)
                        ebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       ebin(nb,nac(n),nn-1+noerrhis(n))+y
                        zebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       zebin(nb,nac(n),nn-1+noerrhis(n))+1
                      enddo
                  endif
                    zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
                  endif
                  itrevt(n)=.true.
                endif
              endif
            endif
            if(ivar(1,n).gt.100.and.ivar(1,n).lt.200)then
              call xval(n,ivar(1,n),ivfra(1,n),j,x)
              sval(1,n)=sval(1,n)+x
            endif
            if(ivar(2,n).gt.100.and.ivar(2,n).lt.200)then
              call xval(n,ivar(2,n),ivfra(2,n),j,y)
              sval(2,n)=sval(2,n)+y
            endif
            if(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
              call xval(n,ivar(1,n),ivfra(1,n),j,x)
            endif
            if(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
              call xval(n,ivar(2,n),ivfra(2,n),j,y)
            endif
            if(ish.ge.8)write (ifch,*) '   ---> histo n,x,y:',n,x,y
          endif
   9      continue
        enddo
  8     continue        
      enddo
c...........................end loop nptl...........................      


      do n=1,nhis
      if(ivar(1,n).eq.-1.or.ivar(2,n).eq.-1)goto 99

c........check event triggers

       go=.true.
        ncontr=0
        do i=1,ntri(n)
          if(itri(i,n).gt.100)then
            if(itri(i,n).lt.200)then
              x=valtri(i,n)
            else               
              call xval(n,itri(i,n),itfra(i,n),0,x) 
            endif
            if(ntrc(i,n).ne.-1)then
              call triggercondition(i,n,x,go)
            else
              ncontr=ncontr+1 
              goo(ncontr)=.true.
              call triggercondition(i,n,x,goo(ncontr))
            endif  
          endif
        enddo
        
c........event variables > 200

        if(go)then
          if(ivar(1,n).gt.100)then
            if(ivar(1,n).gt.200.and.ivar(1,n).lt.300)then
              call xval(n,ivar(1,n),ivfra(1,n),0,x)
            elseif(ivar(1,n).gt.300.and.ivar(1,n).lt.400)then
              call xval(n,ivar(1,n),ivfra(1,n),nptl+1,x)
            elseif(ivar(1,n).gt.100.and.ivar(1,n).lt.200)then
              x=sval(1,n)
            else
              call xval(n,ivar(1,n),ivfra(1,n),0,x)
            endif
            if(ivar(2,n).gt.200.and.ivar(2,n).lt.300)then
              call xval(n,ivar(2,n),ivfra(2,n),0,y)
            elseif(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
              call xval(n,ivar(2,n),ivfra(2,n),nptl+1,y)
            elseif(ivar(2,n).gt.0.and.ivar(2,n).lt.200)then
              y=sval(2,n)
            else             !inom>500
              call xval(n,ivar(2,n),ivfra(2,n),0,y) 
            endif
c The following doesn't work for ivar(2,n)<100, since particle number is not defined !
c            if(ivar(2,n).gt.200.and.ivar(2,n).lt.300)then
c              call xval(n,ivar(2,n),ivfra(2,n),0,y)
c            elseif(ivar(2,n).gt.300.and.ivar(2,n).lt.400)then
c              call xval(n,ivar(2,n),ivfra(2,n),nptl+1,y)
c            elseif(ivar(2,n).gt.100.and.ivar(2,n).lt.200)then
c              y=sval(2,n)
c            else
c              call xval(n,ivar(2,n),ivfra(2,n),0,y) 
c            endif
            if(mod(inorm(n)/100,10).eq.1)y=y*x
            if(mod(inorm(n)/100,10).eq.2.and.x.ne.0.)y=y/x/2/pi
            if(mod(inorm(n)/100,10).eq.4.and.x.ne.0.)y=y/x**1.5
            if(mod(inorm(n)/100,10).eq.5.and.x.ne.0.)y=y/x                  
            sval(1,n)=x
            sval(2,n)=y
            if(ish.ge.6) then
              write (ifch,*) 'histo n,x,y:',n,x,y
            endif
            if(x.ge.xmin(n).and.x.le.xmax(n))then
              if(ilog(n))then
                nb=1+int(log(x/xmin(n))*xinc(n)) 
              else
                nb=1+int((x-xmin(n))*xinc(n))
              endif
              bin(nb,nac(n),n)=bin(nb,nac(n),n)+y
              if(ncontr.gt.0)then
                do nn=1,ncontr
                  if(goo(nn))
     &                    bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &             bbin(nb,nac(n),lookcontr(n)-1+nn)+y
                enddo
              endif
              zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
              itrevt(n)=.true.
            endif
          endif
        endif
        
c........particle variables 

        if(go)then
          if(ivar(1,n).le.100)then
            if(ncontr.gt.0)then  !event trigger contr
              do nb=1,nbin(n)
                do nn=1,ncontr
                  if(goo(nn))
     &                     bbin(nb,nac(n),lookcontr(n)-1+nn)=
     &              bbin(nb,nac(n),lookcontr(n)-1+nn)
     &              +bin(nb,nac(n),n)-bin(nb,3-nac(n),n)
                enddo
              enddo        
            endif
          endif
        endif
          
c............event ok (increase ncevt) or not (take copy)
       
        if(go)then
          ncevt(n)=ncevt(n)+1 
          if(ncontr.gt.0)then
            do nn=1,ncontr
              loo=lookcontr(n)-1+nn
              if(goo(nn))
     &        nccevt(loo)=nccevt(loo)+1 
            enddo
          endif
        else
          nac(n)=3-nac(n)       
          itrevt(n)=.false.
        endif
        
 99   continue
      enddo

      call utprix('xana  ',ish,ishini,4)      
      end
      
c--------------------------------------------------------------------             
      subroutine triggercondition(i,n,x,go)             
c--------------------------------------------------------------------             
c ntrc is used to distinguish the different usage of trigger:
c
c    trigger var xmin xmax
c             ntrc=1
c    trigger or n var1 xmin1 xmax1 var2 xmin2 xmax2 ... varn xminn xmaxn
c          1  ntrc=2  
c          2  ntrc=0
c              ...
c         n-1 ntrc=0
c          n  ntrc=3
c    trigger contr n var1 xmin1 xmax1 var2 xmin2 xmax2 ... varn xminn xmaxn
c             ntrc=-1
c--------------------------------------------------------------------             
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/crvar/idlead
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      logical go,gox,ok,goz
                xmn=xmitri(i,n)
                xmx=xmatri(i,n)
                if(xmn.eq.-123456.and.xmx.eq.-123456)then  !for leading part
                  xmn=float(idlead)
                  xmx=float(idlead)
                endif
                pmn=xmitrp(i,n)
                pmx=xmatrp(i,n)
                if(abs(ntrc(i,n)).eq.1)then
                  goz=.true.
                  if(pmn.gt.99.999.and.pmx.gt.99.999)then
                    if(x.lt.xmn.or.x.gt.xmx)goz=.false.
                  else  
                    if(x.lt.xmn-0.5.or.x.gt.xmx+0.5)goz=.false.
                    ok=rangen().le.xmitrp(i,n)/100.
                    if(.not.ok.and.x.lt.xmn+0.5)goz=.false.
                    ok=rangen().le.xmatrp(i,n)/100.
                    if(.not.ok.and.x.gt.xmx-0.5)goz=.false.
                  endif
                  if(.not.goz)go=.false.
                else
                  if(ntrc(i,n).eq.2)gox=.false.
                  goz=.true.
                  if(pmn.gt.99.999.and.pmx.gt.99.999)then
                    if(x.lt.xmn.or.x.gt.xmx)goz=.false.
                  else  
                    if(x.lt.xmn-0.5.or.x.gt.xmx+0.5)goz=.false.
                    ok=rangen().le.xmitrp(i,n)/100.
                    if(.not.ok.and.x.lt.xmn+0.5)goz=.false.
                    ok=rangen().le.xmatrp(i,n)/100.
                    if(.not.ok.and.x.gt.xmx-0.5)goz=.false.
                  endif
                  if(goz)gox=.true.
                  if(ntrc(i,n).eq.3.and..not.gox)go=.false.
                endif
                end

c-----------------------------------------------------------------------
      subroutine fillhistoconex(n,x,y,lf,j)   !for conex
c-----------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall

           if(.not.(mod(inorm(n),10).ne.4     
     &      .and.mod(inorm(n),10).ne.6
     &      .and.mod(inorm(n)/100,10).ne.3))return
                    if(ilog(n))then
                      c=(xmax(n)/xmin(n))**(1./real(nbin(n)))
                      nde=nint(1./log10(c))
                      nb=max(1,1+int(log10(x/xmin(n))*nde))
                      xmb=xmin(n)*c**(nb-0.5)
                      if(x.gt.xmb.and.nb.lt.nbin(n))then
                        if(x.gt.xmax(n))
     &                    write(ifmt,*)'xana max ?',x,xmax(n),nb
                        nbx=1
                        xmx=c*xmb
                      elseif(x.lt.xmb.and.nb.gt.1)then
                        if(x.lt.xmin(n))write(ifmt,*)'xana min ?',x,nb
                        nbx=-1
                        xmx=xmb/c
                      else
                        nbx=0
                        xmx=0.
                      endif
                    else
                      c=(xmax(n)-xmin(n))/real(nbin(n))
                      nb=max(1,1+int((x-xmin(n))/c))
                      xmb=xmin(n)+c*(nb-0.5)
                      if(x.gt.xmb)then
                        nbx=1
                        xmx=c+xmb
                      elseif(x.lt.xmb)then
                        nbx=-1
                        xmx=xmb-c
                      else
                        nbx=0
                        xmx=0.
                      endif
                    endif
                    xc=(x-xmx)/(xmb-xmx)
                    xc=max(0.,min(1.,xc))
                    bin(nb,nac(n),n)=bin(nb,nac(n),n)+xc*y
                    if(nbx.ne.0)bin(nb+nbx,nac(n),n)
     &                  =bin(nb+nbx,nac(n),n)+(1.-xc)*y
                    zcbin(nb,nac(n),n)=zcbin(nb,nac(n),n)+1
                    if(inoerr(n).gt.0)then
                      do nn=1,inoerr(n)
                        call xval(n,noerr(noerrhis(n),nn),lf,j,y2)
                        ebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       ebin(nb,nac(n),nn-1+noerrhis(n))+y2
                        zebin(nb,nac(n),nn-1+noerrhis(n))=
     &                       zebin(nb,nac(n),nn-1+noerrhis(n))+1
                      enddo
                    endif
      end                    

c---------------------------------------------------------------------
      subroutine xhis(n)
c---------------------------------------------------------------------
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      dimension xx(mxbin)

      double precision histoweight
      common/chiswei/histoweight
      common/cyield/yield
      common/csigma/sigma
      double precision dcel
      common/ems3/dcel,ad
      common/geom/rmproj,rmtarg,bmax,bkmx
      save cnormx

      if(ivar(1,n).eq.-1)then
        nrbins=0
        goto9999
      endif

c.......here normalization.......................................
c           see also   "..........fill histogram"
c.................................................................
c     the norm ( inorm(n) ) is a number hijk which normalizes to:
c     
c  k  0:  * 1
c     1:  / number of events
c     2:  / number of triggered events
c     4:  / bin-counts
c     6:  / number of summed bin-counts (yield=1.)
c     7:  uses same normalization as one histo before
c     
c  j  0:  * 1
c     1:  / bin-width
c     2:  * sigma_total / bin-width
c     3:  * sigma_diff / bin-width
c     
c  i  0:  * 1
c     1:  y => y*x
c     2:  y => y/x/2/pi (modified for mt0)
c     3:  kno-scaling
c     4:  y => y/x**1.5
c     5:  y => y/x
c     6:  y => y*xi (for conex, xi=x of the bin)
c     7:  y => y/x/(x-m)
c
c  h  0: normal
c     1: accumulated
c     
c.................................................................
      norm1=mod(inorm(n),10)
      norm2=mod(inorm(n)/10,10)
      norm3=mod(inorm(n)/100,10)
      norm4=mod(inorm(n)/1000,10)
      nctbin=0
      do l=1,nbin(n)
        nctbin=nctbin+zcbin(l,nac(n),n)
        if(norm1.eq.4.and.zcbin(l,nac(n),n).ne.0d0)then
          bin(l,nac(n),n)=bin(l,nac(n),n)/zcbin(l,nac(n),n)
          if(lookcontr(n).gt.0)then
            do loo=lookcontr(n),lookcontrx(n)
              if(zbbin(l,nac(n),loo).ne.0.)
     &           bbin(l,nac(n),loo)=bbin(l,nac(n),loo)
     &            /zbbin(l,nac(n),loo)
            enddo
          endif
        endif
        if(ilog(n))then
          xx(l)=xmin(n)*(xmax(n)/xmin(n))**((float(l)-.5)/nbin(n))
        else
          xx(l)=(float(l)-0.5)*(xmax(n)-xmin(n))/nbin(n)+xmin(n)
        endif
      enddo
      cnorm=1.
      if(norm1.eq.1)cnorm=1./float(nevent)
      if(norm1.eq.2)then
        if(ncevt(n).ne.0)then
          cnorm=1./float(ncevt(n))
        else
          cnorm=0.
        endif
      endif
      if(norm1.eq.6.and.nctbin.ne.0)cnorm=1./float(nctbin)
      if(norm1.eq.7)cnorm=cnormx
      cnormx=cnorm
      if(ntevt.ne.0)
     &   sigma=10.*pi*bmax**2.*nevent/ntevt !total (untriggered) sigma
      if(norm2.eq.3)then      !differential (triggered) sigma
        if(ntevt.ne.0)
     &     sigma=10.*pi*bmax**2.*ncevt(n)/ntevt
      endif
      if(norm3.eq.3)then      !kno
        first=0.
        secnd=0.
        do l=1,nbin(n)
          if(nctbin.ne.0)first=first+xx(l)*zcbin(l,nac(n),n)/nctbin
          if(nctbin.ne.0)secnd=secnd
     $           +xx(l)**2*zcbin(l,nac(n),n)/nctbin
        enddo
      else
        first=1.
      endif
      if(ilog(n))then
        if(norm2.eq.2.or.norm2.eq.3) cnorm=cnorm*sigma
      else
        if(norm2.ge.1.and.norm2.le.3) cnorm=cnorm*xinc(n)
        if(norm2.eq.2.or.norm2.eq.3) cnorm=cnorm*sigma
      endif
      do l=1,nbin(n)
        if(ilog(n).and.norm2.ge.1.and.norm2.le.3)
     $      bin(l,nac(n),n) =  bin(l,nac(n),n)
     $      /(xmin(n)*exp(float(l)/xinc(n))*(1.-exp(-1./xinc(n))))
        bin(l,nac(n),n) =  bin(l,nac(n),n) * cnorm
      enddo
      f=first
      nrbins=nbin(n)
      nctbin=0
      yield=0.
      shft=0
       if(nint(xpara(1,n)).eq.999963)shft=xpara(2,n)
      do ii=1,nbin(n)
        g=1
        if(norm3.eq.1.and.xx(ii).ne.0.)g=1./xx(ii)
        if(norm3.eq.2)g=2*pi*(xx(ii)+shft)
        if(norm3.eq.4)g=xx(ii)**1.5
        if(norm3.eq.5)g=xx(ii)
        if(norm3.eq.7)g=0 
        yield=yield+bin(ii,nac(n),n)/xinc(n)*hisfac*f*g
      enddo
      do l=1,nbin(n)
        x=(xx(l)+xshift)      !*xhfact
        ar(l,1)=x/f
        sigbin=0
        if(zcbin(l,nac(n),n).ne.0d0)
     *   sigbin=bin(l,nac(n),n)*hisfac*f/sqrt(zcbin(l,nac(n),n))
        if(norm4.eq.0.or.l.eq.1)then
          ar(l,3)=bin(l,nac(n),n)*hisfac*f
          if(lookcontr(n).gt.0)then
           do loo=lookcontr(n),lookcontrx(n)
             r=1
             if(norm1.eq.2.and.nccevt(loo).ne.0.)
     *          r=float(ncevt(n))/nccevt(loo)
             lo=loo-lookcontr(n)+1
             ary(l,lo)=bbin(l,nac(n),loo)*hisfac*f*cnorm*r
             if(zbbin(l,nac(n),loo).gt.0.)then
               ardy(l,lo)=ary(l,lo)/sqrt(zbbin(l,nac(n),loo))
             else
               ardy(l,lo)=0   
             endif
             if(norm1.eq.4)ardy(l,lo)=zbbin(l,nac(n),loo)   
            enddo
          endif
          if(norm3.eq.6)then   !conex
           ar(l,3)=ar(l,3)*xx(l)
          endif
          ar(l,4)=sigbin
        else
          ar(l,3)=ar(l-1,3)+bin(l,nac(n),n)*hisfac*f
          ar(l,4)=sqrt(ar(l-1,4)**2+sigbin**2)
        endif
        if(inoerr(n).ge.1)then
          if(zebin(l,nac(n),noerrhis(n)).gt.0.d0)then
         ar(l,4)=ebin(l,nac(n),noerrhis(n))/zebin(l,nac(n),noerrhis(n))
          else
            ar(l,4)=0.
          endif
        endif
        if(inoerr(n).eq.2)then
          if(zebin(l,nac(n),noerrhis(n)+1).gt.0.d0)then
      ar(l,5)=ebin(l,nac(n),noerrhis(n)+1)/zebin(l,nac(n),noerrhis(n)+1)
          else
            ar(l,5)=0.
          endif
        endif
        if(norm1.eq.4)ar(l,4)=zcbin(l,nac(n),n)           
      enddo
      ionoerr=inoerr(n)
      histoweight=dble(ncevt(n))
      if(norm1.eq.4)histoweight=0d0

 9999 hisfac=1.
      xshift=0
      end
      
c-----------------------------------------------------------------------
      integer function nsdiff(insdif,now)
c-----------------------------------------------------------------------
c  returns  1 if trigger condition for NSD fulfilled and 0 otherwise
c  for  UA1 (insdif=1) or CDF (insdif=2) or STAR (insdif=3,4)  
C  or BRAHMS (insdif=5) 
c  now ... noweak(histogram number)
c-----------------------------------------------------------------------
      include 'epos.inc'
      nsdiff=0
           if(insdif.ge.1)then 
      iii1=0
      iii2=0
      ipos=0
      ineg=0
         do npts=1,nptl
        if(istptl(npts).ne.0)goto 60
        if(   idptl(npts).ne.120 .and.idptl(npts).ne.-120
     *   .and.idptl(npts).ne.130 .and.idptl(npts).ne.-130
     *   .and.idptl(npts).ne.1120.and.idptl(npts).ne.-1120
     *   .and.idptl(npts).ne.1130.and.idptl(npts).ne.-1130
     *   .and.idptl(npts).ne.2230.and.idptl(npts).ne.-2230
     *   .and.idptl(npts).ne.2330.and.idptl(npts).ne.-2330
     *   .and.idptl(npts).ne.3331.and.idptl(npts).ne.-3331)goto 60
c        if(now.eq.1)then  !do not consider weak decay products
c         if(iorptl(npts).ne.0)then
c          idora=abs( idptl(iorptl(npts)) )
c          if(  idora.eq.20   .or.idora.eq.2130
c     &     .or.idora.eq.2230 .or.idora.eq.1130  
c     &     .or.idora.eq.2330 .or.idora.eq.1330  
c     &     .or.idora.eq.3331 )goto 60
c          endif
c         endif
        ppp=sqrt(pptl(1,npts)**2+pptl(2,npts)**2+pptl(3,npts)**2)
        if(ppp.gt.abs(pptl(3,npts)))then
          yyy=.5*log((ppp+pptl(3,npts))/(ppp-pptl(3,npts)))
        else
          yyy=sign(100.,pptl(3,npts))
        endif
        if(insdif.eq.1)then
          if(yyy.gt.2.   .and. yyy.lt.5.6)iii1=1
          if(yyy.gt.-5.6 .and. yyy.lt.-2.)iii2=1
        elseif(insdif.eq.2)then
          if(yyy.gt.3.2  .and. yyy.lt.5.9)iii1=1
          if(yyy.gt.-5.9 .and. yyy.lt.-3.2)iii2=1
          if(yyy.gt.0.   .and. yyy.lt.3.0)ipos=ipos+1
          if(yyy.gt.-3.0 .and. yyy.lt.0. )ineg=ineg+1
        elseif(insdif.eq.3)then
          if(yyy.gt.-5.0 .and. yyy.lt.-3.3 )iii1=1   
          if(yyy.gt. 3.3 .and. yyy.lt. 5.0 )iii2=1   
        elseif(insdif.eq.4)then
          if(yyy.gt.-5.0 .and. yyy.lt.-3.1 )iii1=1   
          if(yyy.gt. 3.1 .and. yyy.lt. 5.0 )iii2=1   
        elseif(insdif.eq.5)then
          if(yyy.gt.-5.25 .and. yyy.lt.-3.26 )iii1=1   
          if(yyy.gt. 3.26 .and. yyy.lt. 5.25 )iii2=1   
        endif
60      continue        
         enddo
        if(insdif.eq.1)then
          if(iii1.eq.1 .and. iii2.eq.1) nsdiff=1
        elseif(insdif.eq.2)then
          if((iii1.eq.1 .and. iii2.eq.1) .and.
     *    ((ipos.ne.0 .and. ineg.ne.0) .and. ipos+ineg.ge.4)) nsdiff=1
        elseif(insdif.eq.3.or.insdif.eq.4.or.insdif.eq.5)then
          if(iii1.eq.1 .and. iii2.eq.1) nsdiff=1
        endif
           else
      stop'in nsdiff. argument of nsdiff not authorized.        '  
           endif
      end

c----------------------------------------------------------------------
      subroutine xtrans(cvar,inom,ifr,n)
c----------------------------------------------------------------------
      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1
      parameter(mxxhis=50)
      common/varhis/icorrtrig(0:mxxhis),ihardevent(0:mxxhis)
     &,ijetfind1(0:mxxhis),ijetfind2(0:mxxhis)

      character*6 cvar
      ifr=0
      if(cvar.eq.'numptl')then
        inom=1
      elseif(cvar.eq.'npaptl')then
        inom=2      
      elseif(cvar.eq.'npmptl')then
        inom=3      
      elseif(cvar.eq.'ispptl')then
        inom=4   
      elseif(cvar.eq.'rapx')then
        inom=5   
      elseif(cvar.eq.'iptlfr')then
        inom=6   
      elseif(cvar.eq.'rinp')then
        inom=7   
      elseif(cvar.eq.'eco')then
        inom=8
      elseif(cvar.eq.'absrap')then
        inom=12
      elseif(cvar.eq.'rap')then
        inom=13
      elseif(cvar.eq.'xp')then
        inom=14
      elseif(cvar.eq.'xe')then
        inom=15
      elseif(cvar.eq.'pt')then
        inom=16
      elseif(cvar.eq.'p1a')then
        inom=17
      elseif(cvar.eq.'p2a')then
        inom=18
      elseif(cvar.eq.'xi')then
        inom=19
      elseif(cvar.eq.'xf')then
        inom=20
      elseif(cvar.eq.'t')then
        inom=21
      elseif(cvar.eq.'rapmi')then
        inom=22
      elseif(cvar.eq.'eta')then
        inom=23
      elseif(cvar.eq.'theta')then
        inom=24
      elseif(cvar.eq.'pt2')then
        inom=25
      elseif(cvar.eq.'et')then
        inom=26
      elseif(cvar.eq.'idptl')then
        inom=27
      elseif(cvar.eq.'istptl')then
        inom=28
      elseif(cvar.eq.'mass')then
        inom=29
      elseif(cvar.eq.'idaptl')then
        inom=30
      elseif(cvar.eq.'egy')then
        inom=31
      elseif(cvar.eq.'rapwro')then
        inom=32
      elseif(cvar.eq.'mt')then
        inom=33
      elseif(cvar.eq.'pplus')then
        inom=34
      elseif(cvar.eq.'pminus')then
        inom=35
      elseif(cvar.eq.'p5')then
        inom=36
      elseif(cvar.eq.'pa')then
        inom=37
      elseif(cvar.eq.'sob')then
        inom=38
      elseif(cvar.eq.'idpom')then
        inom=39
      elseif(cvar.eq.'p3a')then
        inom=40
      elseif(cvar.eq.'cmass')then
        inom=41
      elseif(cvar.eq.'arappi')then
        inom=42
      elseif(cvar.eq.'itsptl')then
        inom=50
      elseif(cvar.eq.'ityptl')then
        inom=51
      elseif(cvar.eq.'idoptl')then
        inom=52
      elseif(cvar.eq.'iptl')then
        inom=53
      elseif(cvar.eq.'index')then
        inom=54
      elseif(cvar.eq.'p1')then
        inom=55
      elseif(cvar.eq.'p2')then
        inom=56
      elseif(cvar.eq.'p3')then
        inom=57
      elseif(cvar.eq.'p4')then
        inom=58
      elseif(cvar.eq.'xg')then
        inom=59
      elseif(cvar.eq.'ek')then
        inom=60
      elseif(cvar.eq.'beta')then
        inom=61
      elseif(cvar.eq.'mt0')then
        inom=63
      elseif(cvar.eq.'qsqptl')then
        inom=64
      elseif(cvar.eq.'xelab')then
        inom=65
      elseif(cvar.eq.'hgtc05')then
        inom=66
        imulty1=1             !to switch on the calculation of "Standard variable"
      elseif(cvar.eq.'hadtyp')then
        inom=67
        imulty1=1 
      elseif(cvar.eq.'hgtc1')then
        inom=68
      elseif(cvar.eq.'x4')then
        inom=69
      elseif(cvar.eq.'npn')then
        inom=70
      elseif(cvar.eq.'routp')then
        inom=71
      elseif(cvar.eq.'hgtc3')then
        inom=72
        imulty1=1 
      elseif(cvar.eq.'mu14')then
        inom=73
        imulty1=1 
      elseif(cvar.eq.'delphi')then
        inom=74
        iok=0
        !------------------------------------------------------------
        !icorrtrig sores the histogram numbers of those histograms which 
        !use the delphi variable (and therfore require a call corrtrig
        !------------------------------------------------------------
        do i=1,icorrtrig(0)
         if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1 
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
        icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'v2')then
        inom=75
      elseif(cvar.eq.'pt4')then
        inom=76
      elseif(cvar.eq.'rin')then
        inom=77
      elseif(cvar.eq.'mulevt')then
        inom=101
      elseif(cvar.eq.'etevt')then
        inom=102
      elseif(cvar.eq.'enevt')then
        inom=103
      elseif(cvar.eq.'ev6evt')then
        inom=104
      elseif(cvar.eq.'xenevt')then
        inom=105
      elseif(cvar.eq.'netevt')then
        inom=106
      elseif(cvar.eq.'ptevt')then
        inom=107
      elseif(cvar.eq.'pmxevt')then
        inom=108
      elseif(cvar.eq.'numevt')then
        inom=201
      elseif(cvar.eq.'egyevt')then
        inom=202
      elseif(cvar.eq.'bimevt')then
        inom=203
      elseif(cvar.eq.'xbjevt')then
        inom=204
      elseif(cvar.eq.'qsqevt')then
        inom=205
      elseif(cvar.eq.'yevt')then
        inom=206
      elseif(cvar.eq.'eloevt')then
        inom=207
      elseif(cvar.eq.'nd1evt')then
        inom=208
      elseif(cvar.eq.'nd2evt')then
        inom=209
      elseif(cvar.eq.'theevt')then
        inom=210
      elseif(cvar.eq.'nspevt')then
        inom=211
      elseif(cvar.eq.'nhpevt')then
        inom=212
      elseif(cvar.eq.'sigtot')then
        inom=213
      elseif(cvar.eq.'sigela')then
        inom=214
      elseif(cvar.eq.'sloela')then
        inom=215
      elseif(cvar.eq.'nrgevt')then
        inom=216
      elseif(cvar.eq.'qevt')then
        inom=217
      elseif(cvar.eq.'qtlevt')then
        inom=218
      elseif(cvar.eq.'threvt')then
        inom=220
        ifr=33                  !set thrust-frame 
      elseif(cvar.eq.'omtevt')then
        inom=221
        ifr=33                  !set thrust-frame 
      elseif(cvar.eq.'tmaevt')then
        inom=222
        ifr=33                  !set thrust-frame 
      elseif(cvar.eq.'tmievt')then
        inom=223
        ifr=33                  !set thrust-frame 
      elseif(cvar.eq.'oblevt')then
        inom=224
        ifr=33                  !set thrust-frame 
      elseif(cvar.eq.'sphevt')then
        inom=230
        ifr=32                  !set sph-frame 
      elseif(cvar.eq.'aplevt')then
        inom=231
        ifr=32                  !set sph-frame 
      elseif(cvar.eq.'cpaevt')then
        inom=232
        ifr=34                  !set sph2-frame 
      elseif(cvar.eq.'dpaevt')then
        inom=233
        ifr=34                  !set sph2-frame 
      elseif(cvar.eq.'npoevt')then
        inom=234
      elseif(cvar.eq.'npnevt')then
        inom=235
  
  !....unused 236, 237

      elseif(cvar.eq.'npxevt')then
        inom=238
      elseif(cvar.eq.'mu1evt')then
        inom=240
        imulty1=1 
      elseif(cvar.eq.'muievt')then
        inom=241
        imulty1=1 
      elseif(cvar.eq.'hgtevt')then
        inom=242
        imulty1=1 
      elseif(cvar.eq.'difevt')then
        inom=243
      elseif(cvar.eq.'dixevt')then
        inom=244
      elseif(cvar.eq.'qinevt')then
        inom=250
      elseif(cvar.eq.'qfievt')then
        inom=251
      elseif(cvar.eq.'einevt')then
        inom=252
      elseif(cvar.eq.'efievt')then
        inom=253
      elseif(cvar.eq.'pinevt')then
        inom=254
      elseif(cvar.eq.'pfievt')then
        inom=255
      elseif(cvar.eq.'pxfevt')then    ! leading proton xf in cms
        inom=256
      elseif(cvar.eq.'pi+xf')then     ! pi+xf: pi+ yield at cms xf>0.01 
        inom=257
      elseif(cvar.eq.'pi-xf')then     ! pi-xf: pi- yield at cms xf>0.01 
        inom=258
      elseif(cvar.eq.'sigcut')then
        inom=260
      elseif(cvar.eq.'keu')then
        inom=261
      elseif(cvar.eq.'ked')then
        inom=262
      elseif(cvar.eq.'kes')then
        inom=263
      elseif(cvar.eq.'kolevt')then
        inom=265
      elseif(cvar.eq.'sigsd')then
        inom=266
      elseif(cvar.eq.'nglevt')then
        inom=267        
      elseif(cvar.eq.'kppevt')then   ! collision numbers per participant
        inom=268        
      elseif(cvar.eq.'npievt')then   ! pion + multiplicity per event
        inom=269        
      elseif(cvar.eq.'np2evt')then   ! pion + multiplicity per participant
        inom=270        
      elseif(cvar.eq.'sigdif'.or.cvar.eq.'sigdifr')then
        inom=271
      elseif(cvar.eq.'koievt')then
        inom=272
      elseif(cvar.eq.'ineevt')then
        inom=273
      elseif(cvar.eq.'elaevt')then
        inom=274
      elseif(cvar.eq.'itgevt')then
        inom=275
        iok=0
        do i=1,icorrtrig(0)
          if(icorrtrig(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          icorrtrig(0)=icorrtrig(0)+1 
          if(icorrtrig(0).gt.mxxhis)stop'mxxhis too small'
          icorrtrig(icorrtrig(0))=n
        endif
      elseif(cvar.eq.'hrdevt')then
        inom=276
        iok=0
        do i=1,ihardevent(0)
          if(ihardevent(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          ihardevent(0)=ihardevent(0)+1 
          if(ihardevent(0).gt.mxxhis)stop'mxxhis too small'
          ihardevent(ihardevent(0))=n
        endif
      elseif(cvar(2:6).eq.'j1evt'.or.cvar(2:6).eq.'j2evt')then
        iok=0
        do i=1,ijetfind1(0)
          if(ijetfind1(i).eq.n)iok=1
        enddo
        if(iok.eq.0)then
          ijetfind1(0)=ijetfind1(0)+1 
          if(ijetfind1(0).gt.mxxhis)stop'mxxhis too small'
          ijetfind1(ijetfind1(0))=n
        endif
        if(cvar.eq.'ej1evt')inom=277
        if(cvar.eq.'pj1evt')inom=278
        if(cvar(2:6).eq.'j2evt')then
          iok=0
          do i=1,ijetfind2(0)
            if(ijetfind2(i).eq.n)iok=1
          enddo
          if(iok.eq.0)then
            ijetfind2(0)=ijetfind2(0)+1 
            if(ijetfind2(0).gt.mxxhis)stop'mxxhis too small'
            ijetfind2(ijetfind2(0))=n
          endif
          if(cvar.eq.'ej2evt')inom=279
          if(cvar.eq.'pj2evt')inom=280
        endif
      elseif(cvar.eq.'zppevt')then
        inom=281
      elseif(cvar.eq.'zptevt')then
        inom=282
      elseif(cvar.eq.'***not used***')then
        inom=283        
      elseif(cvar.eq.'nd3evt')then
        inom=284
      elseif(cvar.eq.'nd4evt')then
        inom=285
      elseif(cvar.eq.'nd5evt')then
        inom=287
      elseif(cvar.eq.'mubevt')then
        inom=286
        imulty1=1 
      elseif(cvar.eq.'aimevt')then
        inom=301
      elseif(cvar.eq.'wjbevt')then
        inom=302
      elseif(cvar.eq.'njbevt')then
        inom=303
      elseif(cvar.eq.'djbevt')then
        inom=304
      elseif(cvar.eq.'tjbevt')then
        inom=305
      elseif(cvar.eq.'hjmevt')then
        inom=306
      elseif(cvar.eq.'ljmevt')then
        inom=307
      elseif(cvar.eq.'djmevt')then
        inom=308
      elseif(cvar.eq.'ybal')then
        inom=310
      elseif(cvar.eq.'yabal')then
        inom=310
      elseif(cvar.eq.'sigine')then
        inom=312
      elseif(cvar.eq.'sigiaa')then
        inom=313
      elseif(cvar.eq.'alpdsf')then
        inom=314
      elseif(cvar.eq.'alpdsh')then
        inom=315
      elseif(cvar.eq.'betdsf')then
        inom=316
      elseif(cvar.eq.'betdsh')then
        inom=317
      elseif(cvar.eq.'rexdip')then
        inom=318
      elseif(cvar.eq.'rexdit')then
        inom=319
      elseif(cvar.eq.'m14evt')then
        inom=320
      elseif(cvar.eq.'ht3evt')then
        inom=321
      elseif(cvar.eq.'sigiex')then
        inom=322
      elseif(cvar.eq.'sigdex')then
        inom=323
      elseif(cvar.eq.'sigsex')then
        inom=324
      elseif(cvar.eq.'ox1evt')then
        inom=501
      elseif(cvar.eq.'ox2evt')then
        inom=502
      elseif(cvar.eq.'ox3evt')then
        inom=503
      elseif(cvar.eq.'ox4evt')then
        inom=504
      else
        print *,' '
        print *,'              xtrans: unknown variable ',cvar
        print *,' '
c       inom=-1
        stop
      endif
      end

c----------------------------------------------------------------------
      subroutine xval(n,inom,lf,j,x)
c----------------------------------------------------------------------
c   n ...... histogram index
c   inom ... variable index 
c              1-100 particle variables
c              101-200 accumulative event variables
c              > 200 other event variables
c   lf ..... frame index
c   particle index (used for particle variables)
c----------------------------------------------------------------------
      include 'epos.inc'
      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1
      parameter(mxxhis=50)
      common/varhis/icorrtrig(0:mxxhis),ihardevent(0:mxxhis)
     &,ijetfind1(0:mxxhis),ijetfind2(0:mxxhis)
      common/zeus2/qtl

      parameter (ntim=1000)
      common/cprt/nprtj,pprt(5,ntim),idprt(ntim),iorprt(ntim)
     &     ,idaprt(2,ntim)

      common/cxyzt/xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
     *,optl(mxptl),uptl(mxptl),sptl(mxptl),rptl(mxptl,3)

      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 
      double precision ebin,zebin
      common/errbins/ebin(mxbin,2,mxhis/2),zebin(mxbin,2,mxhis/2),
     $inoerr(mxhis),noerr(mxhis/2,2),noerrhis(mxhis/2),noerrall
      parameter (mxfra=5)
      common/pfra/nfra,ifra(mxfra),ivfra(2,mxhis),itfra(mxtri,mxhis)
     $     ,imofra(3,mxfra),iffra(mxfra),r1fra(3,mxfra),r2fra(3,mxfra)
     $     ,emax(mxfra)
 
      double precision bofra
      common/dfra/bofra(5,mxfra)
      dimension p(5,mxfra),aimuni(10,mxhis),xor(5,mxfra)
      save p,aimuni,xor
      if(iffra(lf).eq.0.and.j.ne.0)then
        do l=1,5
          p(l,lf)=pptl(l,j)
        enddo         
        do l=1,4
          xor(l,lf)=xorptl(l,j)
        enddo         
        if(imofra(1,lf).ne.0)then
          call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $         ,p(1,lf),p(2,lf),p(3,lf))
          call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $         ,xor(1,lf),xor(2,lf),xor(3,lf))
        endif
        if(imofra(2,lf).ne.0)then !the x-z exchanged is ok !!
          call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $         ,p(3,lf),p(2,lf),p(1,lf))
          call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $         ,xor(3,lf),xor(2,lf),xor(1,lf))
        endif
        if(imofra(3,lf).ne.0)then
          call utlob4(imofra(3,lf),bofra(1,lf),bofra(2,lf) 
     $         ,bofra(3,lf) ,bofra(4,lf),bofra(5,lf)     
     $         ,p(1,lf),p(2,lf),p(3,lf),p(4,lf))
          call utlob4(imofra(3,lf),bofra(1,lf),bofra(2,lf) 
     $         ,bofra(3,lf) ,bofra(4,lf),bofra(5,lf)     
     $         ,xor(1,lf),xor(2,lf),xor(3,lf),xor(4,lf))
        endif
        iffra(lf)=1
      endif

c--------------------------------- 1 - 100 ----------------------------
      if(inom.eq.1)then
        x=1.
      elseif(inom.eq.2)then
        x=isign(1,idptl(j))
      elseif(inom.eq.3)then
        chrg=0
        if(iabs(idptl(j)).le.9999
     $       .and.mod(iabs(idptl(j)),10).le.1)
     $       call idchrg(idptl(j),chrg)
        if(chrg.eq.0.)then
          x=0
        else
          x=int(sign(1.,chrg))
        endif
      elseif(inom.eq.4)then
        iad=abs(idptl(j))
        jspin=mod(iad,10)
        x=0.
        if (iad.ge.100.and.iad.lt.1000) x=1./(1.+2.*jspin)
        if (iad.ge.1000.and.iad.lt.9999) x=1./(2.+2*jspin)
      elseif(inom.eq.5)then    !'rapx'  !st-rap for string segments only !!!!!!!!!!
        x=dezptl(j)
      elseif(inom.eq.6)then                                      !'iptlfr'
        x=0
        if(j.ge.minfra.and.j.le.maxfra)x=1  
      elseif(inom.eq.7)then                                        !'rinp'
        aa=cos(phievt)
        bb=sin(phievt)
        x=xptl(j)*aa+yptl(j)*bb
      elseif(inom.eq.8)then                        !'eco' !engy in comoving frame
        x=0
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then
          amt=sqrt(amt)
          rap=sign(1.,p(3,lf))*alog((p(4,lf)+abs(p(3,lf)))/amt)
          rapx=dezptl(j)
          x=amt*cosh(rap-rapx)
        endif
      elseif(inom.eq.12)then                                      !'absrap'
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then
          amt=sqrt(amt)
          x=alog((p(4,lf)+abs(p(3,lf)))/amt)
        else
          x=0.                  !
        endif
      elseif(inom.eq.13)then    !'rap'
        amt=p(5,lf)**2+p(1,lf)**2+p(2,lf)**2
        if(amt.gt.0..and.p(4,lf)+abs(p(3,lf)).gt.0.d0)then
          amt=sqrt(amt)
          x=sign(1.,p(3,lf))*alog((p(4,lf)+abs(p(3,lf)))/amt)
        else
          x=0.                  !
        endif
      elseif(inom.eq.14)then                                         !'xp'
        x=sqrt(p(3,lf)**2+p(2,lf)**2+p(1,lf)**2)/emax(lf)
      elseif(inom.eq.15)then    !'xe'
        x=p(4,lf)/emax(lf)
      elseif(inom.eq.16)then                                         !'pt'
        x=sqrt(p(2,lf)**2+p(1,lf)**2)
      elseif(inom.eq.17)then
        x=abs(p(1,lf))
      elseif(inom.eq.18)then
        x=abs(p(2,lf))
      elseif(inom.eq.19)then
        x=-log(sqrt(p(3,lf)**2+p(2,lf)**2+p(1,lf)**2)/emax(lf))
      elseif(inom.eq.20)then                                     !'xf'
        m=mod(ifra(lf)/10,10)
        if(m.eq.1)then
c          pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
          if(ifra(lf).eq.12)pmax=pnll
          x=p(3,lf)/pmax
        else
          x=p(3,lf)/emax(lf)
        endif
      elseif(inom.eq.21)then
        pmax=pmxevt
c        pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        if(ifra(lf).eq.12)pmax=pnll
        x=-(amproj**2-2.*sqrt(amproj**2+pmax**2)*p(4,lf)
     *      +2.*abs(pmax*p(3,lf))+p(5,lf)**2)
      elseif(inom.eq.22)then
        amt=sqrt(p(5,lf)**2+p(1,lf)**2+p(2,lf)**2)
        if(amt.ne.0.)then
          x=-sign(1.,p(3,lf))*alog((p(4,lf)+abs(p(3,lf)))/amt)
        else
          x=0.                  !
        endif
      elseif(inom.eq.23)then                                     !'eta'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        x=0
        if(p(3,lf).ne.0..and.pt.ne.0.)x=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
      elseif(inom.eq.24)then                                     !'theta'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        x=90 
        if(p(3,lf).ne.0.)x=atan(pt/p(3,lf))/pi*180.
        if(x.lt.0.)x=180.+x
      elseif(inom.eq.25)then                                     !'pt2'
        x=p(2,lf)**2+p(1,lf)**2
      elseif(inom.eq.26)then                                     !'et'
        pt=sqrt(p(2,lf)**2+p(1,lf)**2)
        x=0
        eef=p(4,lf)
c        if(idptl(j).ge.1000)eef=eef-prom
c        if(idptl(j).le.-1000)eef=eef+prom
        p2=p(3,lf)**2+p(2,lf)**2+p(1,lf)**2
        if(p2.ne.0.)x=eef*pt/sqrt(p2)
      elseif(inom.eq.27)then                                     !'idptl'
        x=idptl(j)
      elseif(inom.eq.28)then    !istptl
        x=istptl(j)
      elseif(inom.eq.29)then    !mass
        x=p(5,lf)
        if(istptl(j).le.1)call idmass(idptl(j),x)
      elseif(inom.eq.30)then    !idaptl
        x=abs(idptl(j))
      elseif(inom.eq.31)then    !egy
        x=egyevt
      elseif(inom.eq.32)then    !arapwro
        x=0
        pt2=p(2,lf)**2+p(1,lf)**2
        if(p(3,lf).ne.0.)x=sign(1.,p(3,lf))*
     *       alog((sqrt(p(3,lf)**2+pt2+.13957**2)+abs(p(3,lf)))
     *       /sqrt(pt2+.13957**2)) 
      elseif(inom.eq.33)then                                  !'mt' 
        x=sqrt(p(2,lf)**2+p(1,lf)**2+p(5,lf)**2)
      elseif(inom.eq.34)then                                  !'pplus' 
        x=sign(1.,p(3,lf)) * (p(4,lf)+abs(p(3,lf)))
      elseif(inom.eq.35)then                                  !'pminus' 
        x=sign(1.,p(3,lf)) * (p(4,lf)-abs(p(3,lf)))
      elseif(inom.eq.36)then                                  !'p5' (mass)
        x=p(5,lf)
      elseif(inom.eq.37)then                                  !pa
        x=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
      elseif(inom.eq.38)then                                  !'pa'
        if(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2.ne.0)
     *       x=egyevt**2/sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)*p(4,lf)
      elseif(inom.eq.39)then                                  !idpom
        x=idptl(j)/1000000
      elseif(inom.eq.40)then                                  !p3a
        x=abs(p(3,lf))
      elseif(inom.eq.41)then                                  
        cm2=p(4,lf)**2-p(3,lf)**2-p(2,lf)**2-p(1,lf)**2         !cmass
        x=sign(sqrt(abs(cm2)),cm2)
      elseif(inom.eq.42)then    !arappi
        x=0
        pt2=p(2,lf)**2+p(1,lf)**2
        if(p(3,lf).ne.0.)
     *       x=alog((sqrt(p(3,lf)**2+pt2+.13957**2)+abs(p(3,lf)))
     *       /sqrt(pt2+.13957**2)) 
      elseif(inom.eq.50)then
        x=itsptl(j)
      elseif(inom.eq.51)then
        x=ityptl(j)
      elseif(inom.eq.52)then
        x=0.
        if(iorptl(j).ne.0) x=idptl(iorptl(j))
      elseif(inom.eq.53)then
        x=j
      elseif(inom.eq.54)then                       !sloela
        call idflav(idptl(j),ifl1,ifl2,ifl3,jspin,indx)
        x=indx
      elseif(inom.eq.55)then                       !p_x
        x=p(1,lf)
      elseif(inom.eq.56)then                       !p_y
        x=p(2,lf)
      elseif(inom.eq.57)then                       !p_z
        x=p(3,lf)
      elseif(inom.eq.58)then                       !'e'
        x=p(4,lf)
      elseif(inom.eq.59)then                       !E/p_max
c        pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        if(ifra(lf).eq.12)pmax=pnll
        x=p(4,lf)/pmax
      elseif(inom.eq.60)then                       !e_kin
        x=p(4,lf)-p(5,lf)
      elseif(inom.eq.61)then                       !'beta'
        x=p(3,lf)/p(4,lf)
      elseif(inom.eq.63)then                       !'mt0'
        x=sqrt(p(2,lf)**2+p(1,lf)**2+p(5,lf)**2)-p(5,lf)
      elseif(inom.eq.64)then                       !qsqptl
        x=qsqptl(j)
      elseif(inom.eq.65)then                       !xelab=Elab/Eolab 
        x=p(4,lf)/(engy**2/2/prom-prom)          
        if(x.gt.0.9999999) x=.9999999
      elseif(inom.eq.66)then    !hgtc05 ... charged ptl mult |[c]|<0.5
        x=multc05
      elseif(inom.eq.67)then    !hadtyp ... primary (1) or secondary (2) hadron
        if(j.le.nbdky)then
          x=1
        else
          x=2
        endif
      elseif(inom.eq.68)then    !hgtc1
        x=multc1
      elseif(inom.eq.69)then                       !'x4'
        x=xor(4,lf)
      elseif(inom.eq.70)then                       !'npn'
        x=npjevt+ntgevt    
      elseif(inom.eq.71)then                       !'routp'
        cc=-sin(phievt)
        dd=cos(phievt)
        x=xptl(j)*cc+yptl(j)*dd
      elseif(inom.eq.72)then    !hgtc3 ... charged ptl mult |eta|<3.15  /6.3
        x=multc3/6.3
      elseif(inom.eq.73)then    !mu14 ... charged ptl mult |eta|<1  pt>.4  
        x=multc14
      elseif(inom.eq.74)then    !delphi ... azimuthhal correlation
        x=10000.
        pt=sqrt(p(1,lf)**2+p(2,lf)**2)

        if(nint(ypara(1,n)).ne.0.and.j.ne.nint(ypara(1,n)).and.
     $           pt.gt.0)then 
           phi=sign(1.,p(2,lf))*acos(p(1,lf)/pt)          
               x=phi-ypara(2,n)
           if(x.lt.-2.5*pi)then
            x=x+4*pi
           elseif(x.lt.-0.5*pi)then
            x=x+2*pi
          elseif(x.gt.3.5*pi)then
            x=x-4*pi
          elseif(x.gt.1.5*pi)then
            x=x-2*pi
          endif
        endif  
      elseif(inom.eq.75)then    !'v2'
        aa=cos(phievt)
        bb=sin(phievt)
        cc=-sin(phievt)
        dd=cos(phievt)
        px=p(1,lf)*aa+p(2,lf)*bb
        py=p(1,lf)*cc+p(2,lf)*dd
        pt2=p(2,lf)**2+p(1,lf)**2
        x=0
        if(pt2.gt.0.)x=(px**2-py**2)/pt2
      elseif(inom.eq.76)then                                     !'pt4'
        x=(p(2,lf)**2+p(1,lf)**2)**2
      elseif(inom.eq.77)then                                     !'rin'
        x=rinptl(j)


c--------------------------------- 101 - 200 ----------------------------

      elseif(inom.eq.101)then           !mulevt
        x=1.
      elseif(inom.eq.102)then                      !'etevt'
        x=0
        if(istptl(j).eq.0)then
         eef=p(4,lf) 
         if(idptl(j).ge.1000)eef=eef-prom
         if(idptl(j).le.-1000)eef=eef+prom
         pp=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
         if(pp.ne.0.)x=eef*sqrt(p(1,lf)**2+p(2,lf)**2)/pp
         if(x.ne.x)then
           write(ifch,*)x,eef,p(1,lf),p(2,lf),p(3,lf),pp,prom,idptl(j),j
           call alist('xan&',1,nptl)
           stop 'probleme dans xan'
         endif
         endif
      elseif(inom.eq.103)then
        x=p(4,lf)/1000.
      elseif(inom.eq.104)then                       !'ev6evt'
        x=0
        if(istptl(j).eq.0)then
         pt=sqrt(p(2,lf)**2+p(1,lf)**2)
         eta=0
         if(p(3,lf).ne.0..and.pt.ne.0.)eta=sign(1.,p(3,lf))*
     *   alog((sqrt(p(3,lf)**2+pt**2)+abs(p(3,lf)))/pt)
         if(pt.eq.0.)eta=sign(1e5,p(3,lf))
         if(eta.gt.6.0)then
          eef=p(4,lf) 
          if(idptl(j).ge.1000)eef=eef-prom
          if(idptl(j).le.-1000)eef=eef+prom
          pp=sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)
          if(pp.ne.0.)x=0.001*eef
         endif 
        endif
      elseif(inom.eq.105)then
        etot=maproj*emax(lf)+matarg*0.94  !nur richtig fur target frame!!!!!
        x=p(4,lf)/etot
      elseif(inom.eq.106)then
        x=isign(1,idptl(j))
      elseif(inom.eq.107)then                       !'ptevt'
        x=sqrt(p(2,lf)**2+p(1,lf)**2)             
      elseif(inom.eq.108)then                       !'pmxevt'
        x=pmxevt             
        
c--------------------------------- > 200 ----------------------------

      elseif(inom.eq.201)then
        x=1.
      elseif(inom.eq.202)then
        x=egyevt
      elseif(inom.eq.203)then
        x=bimevt
      elseif(inom.eq.204)then                       !'xbjevt'
        x=xbjevt
      elseif(inom.eq.205)then                       !'qsqevt'
        x=qsqevt
      elseif(inom.eq.206)then                       !'yevt'
        x=qsqevt/xbjevt/engy**2
      elseif(inom.eq.207)then                       !'eloevt'
        x=qsqevt/4./elepti+elepti*(1.-qsqevt/xbjevt/engy**2)
      elseif(inom.eq.208)then                       !nd1evt
        x=nsdiff(1,noweak(n))
      elseif(inom.eq.209)then                       !'nd2evt'
        x=nsdiff(2,noweak(n))
      elseif(inom.eq.284)then                       !'nd3evt'
        x=nsdiff(3,noweak(n))
      elseif(inom.eq.285)then                       !'nd4evt'
        x=nsdiff(4,noweak(n))
      elseif(inom.eq.287)then                       !'nd5evt'
        x=nsdiff(5,noweak(n))
      elseif(inom.eq.210)then                       !'theevt'
        eloevt=qsqevt/4./elepti+elepti*(1.-qsqevt/xbjevt/engy**2)
        x=acos(1-qsqevt/2./elepti/eloevt)/pi*180.
      elseif(inom.eq.211)then                       !'nspevt'
        x=0
        do i=1,nptl
         if((istptl(i).eq.30.or.istptl(i).eq.31)
     &      .and.int(idptl(i)/1000000).eq.1)x=x+1
        enddo
      elseif(inom.eq.212)then                       !'nhpevt'
        x=0
        do i=1,nptl
         if((istptl(i).eq.30.or.istptl(i).eq.31)
     &      .and.int(idptl(i)/1000000).eq.3)x=x+1
        enddo
      elseif(inom.eq.213)then                       !'sigtot'
        x=sigtot
      elseif(inom.eq.214)then                       !'sigela'
        x=sigela
      elseif(inom.eq.215)then                       !'sloela'
        x=sloela
      elseif(inom.eq.216)then                       !'nrgevt'
        x=0
        do i=1,nptl
          if(istptl(i).eq.31.and.int(idptl(i)/10000).eq.2)x=x+1
        enddo
      elseif(inom.eq.217)then                       !qevt
        x=sqrt(qsqevt)
      elseif(inom.eq.218)then   !qevt
        if(iappl.eq.8)then
          x=qtl
        else
          x=pprt(1,5)
        endif
      elseif(inom.eq.220)then!------------------------------------------
        x=sngl(bofra(1,lf))     !thrust
      elseif(inom.eq.221)then
        x=1.-sngl(bofra(1,lf))  !1-thrust
      elseif(inom.eq.222)then
        x=sngl(bofra(2,lf))     !major
      elseif(inom.eq.223)then
        x=sngl(bofra(3,lf))     !minor
      elseif(inom.eq.224)then
        x=sngl(bofra(2,lf)-bofra(3,lf)) !oblateness
      elseif(inom.eq.230)then!------------------------------------------
        x=1.5*(1.-sngl(bofra(1,lf))) !spherecity
      elseif(inom.eq.231)then
        x=1.5*sngl(bofra(3,lf)) !aplanarity  
      elseif(inom.eq.232)then
        x=3.*sngl(bofra(1,lf)*bofra(2,lf)+bofra(1,lf)*bofra(3,lf)
     &       +bofra(2,lf)*bofra(3,lf)) !c-parameter
      elseif(inom.eq.233)then
        x=27.*sngl(bofra(1,lf)*bofra(2,lf)*bofra(3,lf))   !d-parameter
      elseif(inom.eq.234)then                       !npoevt
        x=0
        do i=1,nptl
         if(istptl(i).eq.30.or.istptl(i).eq.31)x=x+1
        enddo
      elseif(inom.eq.235)then                       !npnevt
        x=npjevt+ntgevt    

   !....unused 236, 237

      elseif(inom.eq.238)then  !npxevt ... nr of pomerons, including absorbed
        x=0
        do i=1,nptl
         if(istptl(i).eq.30.or.istptl(i).eq.31)x=x+1
         if(mod(abs(idptl(i)),100).eq.94)x=x+0.5
        enddo
      elseif(inom.eq.240)then    !mu1evt ... charged ptl multipl for central rap
        x=multy1
      elseif(inom.eq.241)then    !muievt ... charged ptl multipl
        x=multyi
      elseif(inom.eq.242)then    !hgtevt ... charged ptl multipl for central eta
        x=multc05
      elseif(inom.eq.243)then                       !difevt
        npom=0
        do i=1,nptl
         if(istptl(i).eq.30.or.istptl(i).eq.31)npom=npom+1
        enddo
        x=0
        if(npom.eq.0)x=1
      elseif(inom.eq.244)then                       !dixevt
        zpom=0
        do i=1,nptl
         if(istptl(i).eq.30.or.istptl(i).eq.31)zpom=zpom+1
         if(mod(abs(idptl(i)),100).eq.94)zpom=zpom+0.5
        enddo
        x=0
        if(abs(zpom).lt.0.001)x=1
      elseif(inom.eq.250)then
        if(iappl.eq.8)then      !mass in
          x=-pptl(5,6)
        else
          x=pprt(5,3)
        endif
      elseif(inom.eq.251)then
        if(iappl.eq.8)then      !mass out
          x=pptl(5,7)
        else
          x=pprt(5,2)
        endif
      elseif(inom.eq.252)then
        if(iappl.eq.8)then
          x=-pptl(4,6)
        else
          x=pprt(4,2)
        endif
      elseif(inom.eq.253)then
        if(iappl.eq.8)then
          x=pptl(4,7)
        else
          x=pprt(4,3)
        endif
      elseif(inom.eq.254)then
        if(iappl.eq.8)then
          x=abs(pptl(3,6))
        else
          x=abs(pprt(3,2))
        endif
      elseif(inom.eq.255)then
        if(iappl.eq.8)then
          x=abs(pptl(3,7))
          do l=1,5
            p(l,lf)=pptl(l,7)
          enddo         
          if(imofra(1,lf).ne.0)then
            call utrota(imofra(1,lf),r1fra(1,lf),r1fra(2,lf),r1fra(3,lf)
     $           ,p(1,lf),p(2,lf),p(3,lf))
          endif
          if(imofra(2,lf).ne.0)then !the x-z exchanged is ok !!
            call utrota(imofra(2,lf),r2fra(3,lf),r2fra(2,lf),r2fra(1,lf)
     $           ,p(3,lf),p(2,lf),p(1,lf))
          endif
          if(imofra(3,lf).ne.0)then
            call utlob4(imofra(3,lf),bofra(1,lf),bofra(2,lf) 
     $           ,bofra(3,lf) ,bofra(4,lf),bofra(5,lf)     
     $           ,p(1,lf),p(2,lf),p(3,lf),p(4,lf))
          endif
          x=abs(p(3,lf))
        else
          x=abs(pprt(3,3))
        endif
      elseif(inom.eq.256)then  !pxfevt: leading proton xf in cms
        x=-2
c       pmax=sqrt((engy/2.)**2-prom**2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          if(idptl(i).eq.1120.and.istptl(i).eq.0)then
            if(iframe.eq.11)then
              pz=pptl(3,i)
            else
              amt=sqrt(prom**2+pptl(1,i)**2+pptl(2,i)**2)
              rap=alog((pptl(4,i)+pptl(3,i))/amt)
     &           -alog((sqrt(pnll**2+engy**2)+pnll)/engy)
              pz=amt*sinh(rap)
            endif
            x=max(x,pz/pmax)
          endif
        enddo
      elseif(inom.eq.257)then  !  pi+xf: pi+ yield at cms xf>0.01
        x=0. 
c        pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          if(idptl(i).eq.120.and.istptl(i).eq.0)then
            if(iframe.eq.11)then
              pz=pptl(3,i)
            else
              amt=sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2)
              rap=alog((pptl(4,i)+pptl(3,i))/amt)
     &           -alog((sqrt(pnll**2+engy**2)+pnll)/engy)
              pz=amt*sinh(rap) 
            endif
            if(pz/pmax.gt.0.01)x=x+1.
          endif
        enddo
      elseif(inom.eq.258)then  !  pi-xf: pi- yield at cms xf>0.01 
        x=0. 
c        pmax=sqrt((engy/2)**2-prom*2)
          pmax=pnullx               !???????????????????
        do i=1,nptl
          if(idptl(i).eq.-120.and.istptl(i).eq.0)then
            if(iframe.eq.11)then
              pz=pptl(3,i)
            else
              amt=sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2)
              rap=alog((pptl(4,i)+pptl(3,i))/amt)
     &           -alog((sqrt(pnll**2+engy**2)+pnll)/engy)
              pz=amt*sinh(rap) 
            endif
            if(pz/pmax.gt.0.01)x=x+1.
          endif
        enddo
      elseif(inom.eq.260)then!------------------------------
        x=sigcut
      elseif(inom.eq.261)then
        x=keu
      elseif(inom.eq.262)then
        x=ked
      elseif(inom.eq.263)then
        x=kes
      elseif(inom.eq.265)then
        x=kolevt
      elseif(inom.eq.266)then
        x=sigsd
      elseif(inom.eq.267)then
        x=nglevt        
      elseif(inom.eq.268)then  ! kppevt : collision number per participant
        x=kolevt/float(npjevt+ntgevt)         
      elseif(inom.eq.269)then  ! npievt : pion + multi per event
        x=0
        do i=1,nptl
         if(idptl(i).eq.120)x=x+1
        enddo                        
      elseif(inom.eq.270)then  ! np2evt : pion + multi per event
        x=0
        do i=1,nptl
         if(idptl(i).eq.120)x=x+1
        enddo                
        x=x/float(npjevt+ntgevt) 
      elseif(inom.eq.271)then
        x=sigdif
      elseif(inom.eq.272)then  !number of inelastic collisions per event
        x=koievt
      elseif(inom.eq.273)then  ! inelasticity (energy loss of leading particle)
        x=0.
        do i=maproj+matarg+1,nptl
          if(istptl(i).eq.0)then
            if((((abs(idptl(i)).gt.1000.and.abs(idptl(i)).lt.10000)
     *           .and.idproj.gt.1000).or.(iabs(idptl(i)).gt.100
     *           .and.idproj.lt.1000)).and.pptl(4,i)
     *           .gt.x.and.pptl(3,i).gt.0.)x=pptl(4,i)
          endif 
        enddo
        Eini=pptl(4,1)
        if(Eini.gt.0.)x=(Eini-x)/Eini
      elseif(inom.eq.274)then  ! elasticity (energy of leading particle)
        x=0.
        do i=maproj+matarg+1,nptl
          if(istptl(i).eq.0)then
            if((((abs(idptl(i)).gt.1000.and.abs(idptl(i)).lt.10000)
     *           .and.idproj.gt.1000).or.(iabs(idptl(i)).gt.100
     *           .and.idproj.lt.1000)).and.pptl(4,i)
     *           .gt.x.and.pptl(3,i).gt.0.)x=pptl(4,i)
          endif 
        enddo
        Eini=pptl(4,1)
        if(Eini.gt.0.)x=x/Eini
      elseif(inom.eq.275)then         !'itgevt'
        x=0
        if(nint(ypara(1,n)).ne.0)x=1
      elseif(inom.eq.276)then         !'hrdevt' ......  1 = hard event
        x=0
        if(nint(ypara(1,n)).ne.0)x=1
      elseif(inom.eq.277)then         !'ej1evt' .... et of jet 1       
        x=0
        if(nint(ypara(1,n)).ne.0)
     &  x=ypara(2,n)
      elseif(inom.eq.278)then         !'pj1evt' .... phi of jet 1      
        x=1000
        if(nint(ypara(1,n)).ne.0)
     &  x=ypara(4,n)
      elseif(inom.eq.279)then         !'ej2evt' .... et of jet 2       
        x=0
        if(nint(ypara(6,n)).ne.0)
     &  x=ypara(7,n)
      elseif(inom.eq.280)then         !'pj2evt' .... delta_phi of jet 2 1      
        x=1000
        if(nint(ypara(6,n)).ne.0)then
          x=ypara(9,n)-ypara(4,n)
           if(x.lt.-2.5*pi)then
            x=x+4*pi
           elseif(x.lt.-0.5*pi)then
            x=x+2*pi
          elseif(x.gt.3.5*pi)then
            x=x-4*pi
          elseif(x.gt.1.5*pi)then
            x=x-2*pi
          endif
        endif
      elseif(inom.eq.281)then         !'zppevt'  
        x=zppevt    
      elseif(inom.eq.282)then         !'zptevt'  
        x=zptevt    
      elseif(inom.eq.283)then         
        stop '**********not used*********'        
      elseif(inom.eq.286)then         !'mubevt'
        x=multeb        
      elseif(inom.eq.298)then
        x=sval(1,n)
      elseif(inom.eq.299)then
        x=sval(2,n)
      elseif(inom.eq.301)then   !---------------------------------------------
        if(j.eq.0)then          !initialize
          do l=1,4
            aimuni(l,n)=0.
          enddo
        elseif(j.gt.nptl)then   !final calculation
          am2=aimuni(4,n)**2-aimuni(3,n)**2
     $         -aimuni(2,n)**2-aimuni(1,n)**2
          x=sign(sqrt(abs(am2)),am2)
c          print *, x
        else                    !routine work
          do l=1,4
            aimuni(l,n)=aimuni(l,n)+p(l,lf)
          enddo
c          print *, j,(p(l,lf),l=1,5)
        endif
      elseif(inom.ge.302.and.inom.le.305)then   !-----------------------
        if(j.eq.0)then          !initialize
          do l=1,4
            aimuni(l,n)=0.
          enddo
        elseif(j.gt.nptl)then   !final calculation
          if(inom.eq.302) x=max(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         ,aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.303) x=min(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         ,aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.304) x=abs(aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         -aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n)))
          if(inom.eq.305) x=aimuni(1,n)/2/(aimuni(2,n)+aimuni(4,n))
     $         +aimuni(3,n)/2/(aimuni(2,n)+aimuni(4,n))
        else                    !routine work
          l=0
          if(p(3,lf).lt.0.)l=2
          aimuni(1+l,n)=aimuni(1+l,n)+sqrt(p(1,lf)**2+p(2,lf)**2)
          aimuni(2+l,n)=aimuni(2+l,n)
     $         +sqrt(p(1,lf)**2+p(2,lf)**2+p(3,lf)**2)

        endif       
      elseif(inom.eq.306.or.inom.eq.307.or.inom.eq.308)then !---------
        if(j.eq.0)then          !initialize
          do ll=1,8
            aimuni(ll,n)=0.
          enddo
        elseif(j.gt.nptl)then   !final calculation
          am2a=aimuni(4,n)**2-aimuni(3,n)**2
     $         -aimuni(2,n)**2-aimuni(1,n)**2
          am2b=aimuni(8,n)**2-aimuni(7,n)**2
     $         -aimuni(6,n)**2-aimuni(5,n)**2
          if(inom.eq.306)x=(max(0.,am2a,am2b))/engy**2
          if(inom.eq.307)x=(max(0.,min(am2a,am2b)))/engy**2
          if(inom.eq.308)x=(abs(am2a-am2b))/engy**2
        else                    !routine work
          ll=0
          if(p(3,lf).lt.0.)ll=4
          do l=1,4
            aimuni(l+ll,n)=aimuni(l+ll,n)+p(l,lf)
          enddo
        endif
      elseif (inom.eq.310.or.inom.eq.311) then !---------
        if(j.eq.0)then          !initialize
          aimuni(1,n)=0
          aimuni(2,n)=0
          do i=1,nptl
c            charge=0.
             if(istptl(i).eq.0) then 
               if (idptl(i).eq.idcod(1,n)) aimuni(1,n)=aimuni(1,n)+1.
               if (idptl(i).eq.idcod(2,n)) aimuni(2,n)=aimuni(2,n)+1.
             endif
           enddo
        elseif(j.gt.nptl)then   !final calculation
          if(aimuni(1,n).eq.0.or.aimuni(2,n).eq.0) then
            ncevt(n)=ncevt(n)-1
          endif
          x=xmin(n)-100.
          do i=1,nbin(n)
            zcbin(i,nac(n),n)=abs(zcbin(i,nac(n),n))
         enddo
       else                    !routine work
          if( istptl(j).eq.0
     $         .and. aimuni(1,n).ne.0. .and. aimuni(2,n).ne.0. ) then
            id1=idptl(j)
            if(id1.eq.idcod(1,n) .or. id1.eq.idcod(2,n)) then
              y1=sign(1.,pptl(3,j))*alog((pptl(4,j)+abs(pptl(3,j)))
     *             /sqrt(pptl(5,j)**2+pptl(1,j)**2+pptl(2,j)**2))
              do i=1,nptl
                if(i.eq.j .or. istptl(i).ne.0) goto 124
                id2=idptl(i)
                if(id2.eq.idcod(1,n) .or. id2.eq.idcod(2,n)) then
                  y2=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))
     *                 /sqrt(pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2))
                  dy=(y2-y1)
                  if(inom.eq.311) dy=abs(dy)
                  ib=1+int((dy-xmin(n))*xinc(n))
                  if(dy.ge.xmin(n).and.dy.le.xmax(n)) then
                    if( id1.eq.idcod(1,n) ) then
                      if( id2.eq.idcod(2,n) ) then 
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)+.5/aimuni(2,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)+1
                      else
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)-.5/aimuni(1,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)-1
                      endif
                    else        !id1 is idcod(2,n)  
                      if(id2.eq.idcod(1,n)) then 
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)+.5/aimuni(1,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)+1
                      else
                        bin(ib,nac(n),n)=bin(ib,nac(n),n)-.5/aimuni(2,n)
                        zcbin(ib,nac(n),n)=zcbin(ib,nac(n),n)-1
                      endif
                    endif
                  endif
                endif
 124          enddo
            endif
          endif    
        endif
      elseif (inom.eq.312) then !---------
        x=sigine
      elseif (inom.eq.313) then !---------
        x=sigineaa
      elseif (inom.eq.314) then !---------
        x=alpD(idxD0,iclpro,icltar)
      elseif (inom.eq.315) then !---------
        x=alpD(1,iclpro,icltar)
      elseif (inom.eq.316) then !---------
        x=betD(idxD0,iclpro,icltar)
        if(x.lt.0.)x=-10.*x
      elseif (inom.eq.317) then !---------
        x=betD(1,iclpro,icltar)
      elseif (inom.eq.318) then !---------
        x=rexdif(iclpro)
      elseif (inom.eq.319) then !---------
        x=rexdif(icltar)
      elseif(inom.eq.320)then    !m14evt ... multipl |eta|<1, pt>0.4
        x=multc14
      elseif(inom.eq.321)then    !ht3evt ... height |eta|<3.15
        x=multc3/6.3
      elseif (inom.eq.322) then !---------
        x=sigineex
      elseif (inom.eq.323) then !---------
        x=sigdifex
      elseif (inom.eq.324) then !---------
        x=sigsdex
      elseif (inom.eq.501) then !---------
        x=sval(1,1)
      elseif (inom.eq.502) then !---------        
        x=sval(1,2)
      elseif (inom.eq.503) then !---------        
        x=sval(1,3)
      elseif (inom.eq.504) then !---------        
        x=sval(1,4)
      endif                     !---------------------------------------
      end

c----------------------------------------------------------------------
      subroutine StandardVariables
c----------------------------------------------------------------------
      include 'epos.inc'
      common/stavar/multc05,multy1,multc14,multyi,multc3,imulty1,multeb
     &     ,multc1
      common/crvar/idlead
      parameter(mxxhis=50)
      common/varhis/icorrtrig(0:mxxhis),ihardevent(0:mxxhis)
     &,ijetfind1(0:mxxhis),ijetfind2(0:mxxhis)

      Emax=0
      multy1=0
      multc05=0
      multc14=0
      multc1=0
      multyi=0
      multc3=0
      multeb=0

      
      do i=maproj+matarg+1,nptl
c---multy1---charged ptl multipl for central rap
        if(istptl(i).eq.0)then
          amt=pptl(5,i)**2+pptl(1,i)**2+pptl(2,i)**2
          pt=pptl(1,i)**2+pptl(2,i)**2
          pp=sqrt(pptl(1,i)**2+pptl(2,i)**2+pptl(3,i)**2)
          if(amt.gt.0..and.pptl(4,i).gt.0.)then
            amt=sqrt(amt)
            rap=sign(1.,pptl(3,i))*alog((pptl(4,i)+abs(pptl(3,i)))/amt)
          else
            rap=1000.          
          endif
          if(pt.gt.0.)then
            pt=sqrt(pt)
            eta=sign(1.,pptl(3,i))*alog((pp+abs(pptl(3,i)))/pt)
          else
            eta=1000.          
          endif
          if(abs(idptl(i)).ge.100
     $         .and.abs(idptl(i)).lt.10000)then
            call idchrg(idptl(i),ch)
            if(abs(ch).gt.0.1)then
c---multyi---charged ptl multipl 
              multyi=multyi+1
c---multy1---charged ptl multipl for central rap
              if(abs(rap).le.1.)multy1=multy1+1
              if(abs(eta).le.0.5)multc05=multc05+1
              if(abs(eta).le.1.and.pt.gt.0.4)multc14=multc14+1
              if(abs(eta).le.1)multc1=multc1+1
              if(abs(rap).le.3.15)multc3=multc3+1
c---multeb---charged ptl multipl for back rap
              if(eta.gt.-3.8.and.eta.lt.-2.8)multeb=multeb+1
            endif
          endif
          if((((abs(idptl(i)).gt.1000.and.abs(idptl(i)).lt.10000)
     *         .and.idproj.gt.1000).or.(iabs(idptl(i)).gt.100
     *         .and.idproj.lt.1000)).and.pptl(4,i)
     *         .gt.Emax.and.pptl(3,i).gt.0.)then
            Emax=pptl(4,i)
            idlead=i
          endif
        endif
      enddo
      end

c----------------------------------------------------------------------
      subroutine jetfind(m,n)
c----------------------------------------------------------------------
c   m = 1 ou 2 (two different definitions)
c   n = histogram
c input(jet definition):     
c   xpara(1,n) ... 1=used (m=1)
c   xpara(2,n) ... etamin (m=1)
c   xpara(3,n) ... etamax (m=1)
c   xpara(4,n) ... rmax   (m=1) (rmax defining the cone)     
c   xpara(5,n) ... ichd   (m=1) (1=charged, 0=all)
c   xpara(6,n) ... 1=used (m=2) 
c   xpara(7,n) ... etamin (m=2)
c   xpara(8,n) ... etamax (m=2)
c   xpara(9,n) ... rmax   (m=2)       
c   xpara(10,n) .. ichd   (m=2) 
c output (jet properties):
c   ypara(1,n) ... 1 (found) or 0 if not  (m=1)
c   ypara(2,n) ... et                     (m=1)
c   ypara(3,n) ... eta of center          (m=1)
c   ypara(4,n) ... phi of center          (m=1)
c   ypara(5,n) 
c   ypara(6,n) ... 1 (found) or 0 if not  (m=2) 
c   ypara(7,n) ... et                     (m=2)
c   ypara(8,n) ... eta of center          (m=2)
c   ypara(9,n) ... phi of center          (m=2)
c   ypara(10,n) 
c----------------------------------------------------------------------
      
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      parameter (mxval=5)
      real ptx(mxval),lst(mxval),etax(mxval),phix(mxval)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 

      if(m.ne.1.and.m.ne.2)stop'jetfind: value of m not valid.      '
      
      etamin=           xpara(2+5*(m-1),n)
      etamax=           xpara(3+5*(m-1),n)
      rmax  =           xpara(4+5*(m-1),n)
      ichd  = nint(xpara(5+5*(m-1),n))

      ifound=0  
      do l=1,mxval
        ptx(l)=0
        lst(l)=0
        etax(l)=0
        phix(l)=0
      enddo

ctp060829      pp1=0
ctp060829      pp2=0
ctp060829      pp3=0

      do i=maproj+matarg+1,nptl
        iok=0
        if(istptl(i).eq.0.and.abs(idptl(i)).lt.10000)iok=1
        if(iok.eq.1)call idchrg(idptl(i),ch)
        if(ichd.eq.1.and.nint(ch).eq.0)iok=0
        if(iok.eq.1)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)          
          pt=sqrt(p1**2+p2**2) 
                if(pt.gt.0)then
            eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
            phi=sign(1.,p2)*acos(p1/pt)
          else
            eta=10000
            phi=0
          endif    
          do k=1,mxval
            iok=1
            if(m.eq.2)then
              dphi=phi-ypara(4,n)
              if(dphi.lt.-pi)dphi=dphi+2*pi
              if(dphi.gt. pi)dphi=dphi-2*pi
              if(abs(dphi).lt.pi/2)iok=0
            endif
            if(iok.eq.1.and.pt.gt.ptx(k)
     &        .and.eta.le.etamax.and.eta.ge.etamin)then
              do l=mxval,k+1,-1
               ptx(l)=ptx(l-1)
               lst(l)=lst(l-1)
               etax(l)=etax(l-1)
               phix(l)=phix(l-1)
              enddo
               ptx(k)=pt
               lst(k)=i
               etax(k)=eta
               phix(k)=phi
              goto2
            endif        
          enddo
  2       continue  
        endif
      enddo        

      kk=0           
      etx=0

      do k=1,mxval
       if(lst(k).ne.0)then

        ifound=1
        et=0
        etaxx=etax(k)
        phixx=phix(k)
        do j=maproj+matarg+1,nptl
          iok=0
          if(istptl(j).eq.0.and.abs(idptl(j)).lt.10000)iok=1
          if(iok.eq.1)call idchrg(idptl(j),ch)
          if(ichd.eq.1.and.nint(ch).eq.0)iok=0
          if(iok.eq.1)then
            p1=pptl(1,j)
            p2=pptl(2,j)
            p3=pptl(3,j)              
            pt=sqrt(p1**2+p2**2)  
            am=pptl(5,j)
                  if(pt.gt.0)then
              eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
              phi=sign(1.,p2)*acos(p1/pt)
            else
              eta=-10000
              phi=0
            endif         
            if(eta.le.etamax.and.eta.ge.etamin)then
              deta=eta-etaxx
              dphi=phi-phixx
              if(dphi.lt.-pi)dphi=dphi+2*pi
              if(dphi.gt. pi)dphi=dphi-2*pi
              if(deta**2+dphi**2.lt.rmax**2)et=et+sqrt(pt**2+am**2)
            endif
          endif
        enddo 
        if(et.gt.etx)then
          etx=et
          kk=k
        endif 

       endif          
      enddo

      ypara(1+5*(m-1),n)=ifound
      ypara(2+5*(m-1),n)=etx
      if(kk.gt.0)then
       ypara(3+5*(m-1),n)=etax(kk)
       ypara(4+5*(m-1),n)=phix(kk)
      endif
      return
      end


c----------------------------------------------------------------------
      subroutine hardevent(n)
c----------------------------------------------------------------------
c   n = histogram
c input(jet event conditions):     
c   xpara(2,n) ... pt1
c   xpara(3,n) ... pt2
c   xpara(4,n) ... absetamax
c   xpara(5,n) ... rmax        (r=sqrt(deltaeta**2+deltaphi**2))    
c   xpara(6,n) ... Et_min  
c output (jet event found or not):
c   ypara(1,n) ... 1 (found) or 0 if not
c----------------------------------------------------------------------
      
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 

      ypara(1,n)=0
      do i=maproj+matarg+1,nptl
       if(abs(idptl(i)).ge.100.and.abs(idptl(i)).lt.10000.
     $  and.istptl(i).eq.0)then
        call idchrg(idptl(i),ch)
        if(abs(ch).gt.0.1)then
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)          
          pt=sqrt(p1**2+p2**2) 
          if(pt.gt.0)then
            eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
            phi=sign(1.,p2)*acos(p1/pt)
          else
            eta=10000
            phi=0
          endif        
          if(pt.ge.xpara(2,n).and.abs(eta).lt.xpara(4,n))then
            et1=pptl(4,i)*pt/sqrt(p3**2+pt**2)        
            do j=maproj+matarg+1,nptl
              if(j.ne.i
     $        .and.abs(idptl(j)).ge.100.and.abs(idptl(j)).lt.10000.
     $        .and.istptl(j).eq.0)then
                call idchrg(idptl(j),ch)
                if(abs(ch).gt.0.1.and.abs(idptl(j)).ge.100
     $          .and.abs(idptl(j)).lt.10000.and.istptl(j).eq.0)then
                  p1=pptl(1,j)
                  p2=pptl(2,j)
                  p3=pptl(3,j)          
                  pt=sqrt(p1**2+p2**2)  
                        if(pt.gt.0)then
                   etax=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
                   phix=sign(1.,p2)*acos(p1/pt)
                  else
                    etax=-10000
                    phix=0
                  endif    
                  if(pt.ge.xpara(3,n).and.abs(etax).lt.xpara(4,n))then
                    deta=eta-etax
                    dphi=phi-phix
                    if(dphi.lt.-pi)dphi=dphi+2*pi
                    if(dphi.gt. pi)dphi=dphi-2*pi
                    if(deta**2+dphi**2.lt.xpara(5,n)**2)then
                    et2=pptl(4,j)*pt/sqrt(p3**2+pt**2)        
                     if(et1+et2.gt.xpara(6,n))then
                      ypara(1,n)=1
                      goto1
                     endif 
                    endif        
                  endif
                endif        
              endif         
            enddo
          endif
        endif
       endif          
      enddo

   1  continue                   
      return
      end

c----------------------------------------------------------------------
      subroutine corrtrig(n)
c----------------------------------------------------------------------
c   n = histogram
c input(trigger conditions):     
c   xpara(2,n) ... ptmin      
c   xpara(3,n) ... ptmax      
c   xpara(4,n) ... etamin     
c   xpara(5,n) ... etamax
c   xpara(6,n) ... 
c   xpara(7,n) ... 
c output (triggered particle (the one with highest pt if there are several)):
c   ypara(1,n) ... iptl or 0 if no particle found
c   ypara(2,n) ... phi of particle 
c----------------------------------------------------------------------
      
      include 'epos.inc'
      parameter (mxhis=500,mxcontr=300,mxidcd=60,mxtri=20,mxbin=400)
      parameter (mypara=10,mxpara=10)
      logical ilog,icnx,itrevt,idmod
      double precision bin,bbin,zcbin,zbbin
      common/bins/bin(mxbin,2,mxhis),zcbin(mxbin,2,mxhis)
     $     ,bbin(mxbin,2,mxcontr),itrevt(mxhis),zbbin(mxbin,2,mxcontr)
     $     ,nac(mxhis),ilog(mxhis),icnx(mxhis),xinc(mxhis),ncevt(mxhis)
     $     ,sval(2,mxhis),valtri(mxtri,mxhis),ntrc(mxtri,mxhis)
     $     ,xmin(mxhis),xmax(mxhis),nhis,noweak(mxhis)
     $     ,ivar(2,mxhis),inorm(mxhis),nbin(mxhis),nidcod(mxhis)
     $     ,idcod(mxidcd,mxhis),idmod(mxidcd,mxhis),ntri(mxhis)
     $     ,itri(mxtri,mxhis),xmitri(mxtri,mxhis),xmatri(mxtri,mxhis)
     $     ,xmitrp(mxtri,mxhis),xmatrp(mxtri,mxhis),xpara(mxpara,mxhis)  
     $     ,ypara(mypara,mxhis),lookcontr(mxhis) 
     $     ,lookcontrx(mxhis),ncontrall,icontrtyp(mxhis),nccevt(mxcontr) 

      pt0=xpara(2,n)

      do i=maproj+matarg+1,nptl
       if(abs(idptl(i)).ge.100.and.abs(idptl(i)).lt.10000.
     $  and.istptl(i).eq.0)then
        call idchrg(idptl(i),ch)
        if(abs(ch).gt.0.1)then   
          p1=pptl(1,i)
          p2=pptl(2,i)
          p3=pptl(3,i)
          pt=sqrt(p1**2+p2**2)  
          pt=max(pt,1e-20)
          eta=sign(1.,p3)*alog((sqrt(p3**2+pt**2)+abs(p3))/pt)
          phi=sign(1.,p2)*acos(p1/pt)
          if(pt.ge.pt0.and.pt.le.xpara(3,n).and.eta.gt.xpara(4,n)
     $      .and.eta.lt.xpara(5,n))then
            pt0=pt
            ypara(1,n)=i
            ypara(2,n)=phi
          endif            
        endif
       endif          
      enddo

      return
      end
      
c----------------------------------------------------------------------
