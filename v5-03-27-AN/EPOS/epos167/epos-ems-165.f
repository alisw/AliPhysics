c-----------------------------------------------------------------------
      subroutine emsaa(iret) 
c-----------------------------------------------------------------------
c  energy-momentum sharing
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/cwzero/wzero,wzerox
      double precision omega,omlog,oma,omb,wab,wba,wmatrix,wzero,nbar
     *,wzerox,rrr,eps,xprem,xmrem,om1intgck
      parameter(eps=1.d-30)
      common/col3/ncol,kolpt
c      logical modu
      common/cems5/plc,s
      double precision s,px,py,pomass,plc!,PhiExpo
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common/cncl3/iactpr(mamx),iacttg(mamx)
      common/nucl3/phi,bimp
      logical vpom
      dimension ishuff(500),icp(2),ict(2)
      call utpri('emsaa ',ish,ishini,4)
      
      ntry2=0
      irea=iret

 0001 continue
      iret=0
      iret2=0

c     initialize
c     ----------

      if(ntry2.gt.0)then
        call conre
        call conwr
      endif
      call emsipt   !initialize projectile and target
      call emsigr   !initialize grid



           if((iactpr(1).eq.1.and.iacttg(1).eq.1.and.maproj+matarg.eq.2)
     &           
     &.or.(iactpr(1).eq.1.and.maproj.eq.1.and.maproj+matarg.gt.2)
     &.or.(iacttg(1).eq.1.and.matarg.eq.1.and.maproj+matarg.gt.2)
     &.or.(maproj.gt.1.and.matarg.gt.1))then !if not nothing


        ievt0=0
      
        nSprmx=0
        do k=1,koll
          nSprmx=nSprmx+nprmx(k)
        enddo
     
        omlog=0
        nemsi=nemsi+1
        if(nemsi.le.4.and.iemsi1.eq.1)call xEmsI1(1,0,omlog) 
        if(ish.ge.6)write (ifch,*)'after xEmsI1'
        if(nemsi.le.4.and.iemsi2.eq.1)call xEmsI2(1,0)
        if(ish.ge.6)write (ifch,*)'after xEmsI2'
        if(ish.ge.6)call XPrint('Before Markov:&')
      

c     Markov
c     ------
     
      if(ish.ge.4)write(ifch,*)'Markov Process'
      kint=int(max(15.,2.*engy**0.2))
      if(koll.gt.50)kint=3*kint/int(log(float(koll)))
      kmcmx=nSprmx*kint        !50*kint  !100*kint  


      do kmc=1,kmcmx               !-----> start Metropolis 

       knprmx=0
       rrr=dble(rangen())
       do ik=1,koll
         knprmx=knprmx+nprmx(ik)
         if(rrr.le.dble(knprmx)/dble(nSprmx))then ! k-th pair
           k=ik
           goto 10
         endif
       enddo
 10    continue

       ip=iproj(k)
       it=itarg(k)
       n=1+int(rangen()*float(nprmx(k)))  ! n-th spot for k-th pair
       nbar=dble(npr(0,k))
       if(idpr(n,k).eq.0)nbar=nbar-1d0

       xprem=1.d0!xpp(ip)+xppr(n,k)        !consistently, it should be 1.
       xmrem=1.d0!xmt(it)+xmpr(n,k)
       wzerox=(nbar+1d0)
       wzero=wzerox    / ( wzerox
     &                    +om1intgck(k,xprem,xmrem)*gammaV(k) )

       if(ish.ge.8)write(ifch,*)'wzero',k,n,wzero,wzerox,gammaV(k)
     &                          ,om1intgck(k,xprem,xmrem)
       if(ish.ge.1.and.100000*(kmc/100000).eq.kmc)
     & write(ifmt,*)'kmc',kmc,kmcmx

       call StoCon(1,k,n) 
       call RemPom(k,n)
       call ProPo(k,n) 
       call ProXY(k,n)

       call StoCon(2,k,n) 

       if(idpr(n,k).eq.0.and.idx0.eq.0)then
         accept=accept+1.
       else   

         omb=omega(n,k)
         if(omb.le.0.d0)then
           reject=reject+1.      
           call RemPom(k,n)
           call StoCon(-1,k,n) 
         else

           wab=wmatrix(k,n)
           if(ish.ge.8)write(ifch,*)'omb',omb,wab,k,n
           if(wab.le.0.d0)then 
             write (ifmt,*)'wab,kmc',wab,omb,kmc,k,n,xpr(n,k),ypr(n,k)
     &  ,xppr(n,k),xmpr(n,k),xpp(ip),xmt(it),ip,it,idpr(n,k)
             write(ifmt,'(a,i12,d25.15)')'ems,seedf',nrevt+1,seedc
             iret=1
             goto 1000
           endif 
           call RemPom(k,n)
           call StoCon(-1,k,n)
           oma=omega(n,k)
           wba=wmatrix(k,n)
           if(oma.ge.0.d0.and.oma.le.eps*omb*wba/wab)then
             accept=accept+1.
             call RemPom(k,n)
             call StoCon(-2,k,n) 
             omlog=omlog+dlog(omb)
             goto 500
           elseif(oma.le.1.d-300.or.oma.ne.oma.or.omb.ne.omb)then 
             write (ifmt,*)'oma,kmc',oma,omb,kmc,k,n,xpr(n,k),ypr(n,k)
     &  ,xppr(n,k),xmpr(n,k),idpr(n,k),npr(1,k),xpp(ip),xmt(it),ip,it
             write(ifmt,'(a,i12,d25.15)')'ems,seedf',nrevt+1,seedc
             iret=1
             goto 1000
           endif

           z=sngl(omb/oma*wba/wab)
           if(ish.ge.8)write(ifch,*)'z,oma',z,oma,wba,k,n
           if(rangen().gt.z)then
             reject=reject+1.      
           else
             accept=accept+1.
             call RemPom(k,n)
             call StoCon(-2,k,n) 
             omlog=omlog-dlog(oma)+dlog(omb)
           endif 

 500       continue

         endif 

         endif 

       if(nemsi.le.4)then
         kplot=int(float(kmc)/float(kmcmx)*100.)
         if(iemsi1.eq.1)call xEmsI1(1,kplot,omlog) 
         if(iemsi2.eq.1)call xEmsI2(1,kplot) 
       endif

      enddo  !-----> end Metropolis


           elseif(iokoll.eq.1)then
             
             ievt0=0
             n=1

             do k=1,koll

       call ProPo(k,n) 
       call ProXY(k,n)

             enddo

           else
        
        ievt0=1
      
           endif


c --- Plot Pomeron b-distributions ---

      if(ish.ge.6)call XPrint('After Markov :&')

           if(ntry2.eq.0)then

      if(iemsb.eq.1)then ! plot
       do k=1,koll
        call xEmsB(1,1,k)
        if(nprt(k).gt.0)call xEmsB(1,2,k)
       enddo
      endif

      if(iemsbg.eq.1.and.ievt0.eq.0)then ! plot
        call xEmsBg(3,0,0)
        do k=1,koll
          call xEmsBg(1,0,k)
          if(nprt(k).gt.0)then
            call xEmsBg(1,-1,k)
            do n=1,nprmx(k)
              if(idpr(n,k).ne.0)call xEmsBg(1,idpr(n,k),k)
            enddo
          endif
        enddo
      endif

c --- Plot distr of pomeron number ---


      if(iemspm.eq.1.and.ievt0.eq.0)then
       do k=1,koll
           call xEmsPm(1,k,nprt(k))
       enddo
      endif
      
           endif

c -- Split Enhanced Pomerons and fix their nature ---

      do k=1,koll
        do n=1,nprmx(k)
          if(idfpr(n,k).eq.0)call ProPoTy(k,n)
        enddo
      enddo

      
c --- Count real interactions ---

      ncol=0
      do k=1,koll
       if(nprt(k).gt.0)then
        ncol=ncol+1
        itpr(k)=1
        ip=iproj(k)
        it=itarg(k)
        kolp(ip)=kolp(ip)+nprt(k)        !number of cut Pomerons
        kolt(it)=kolt(it)+nprt(k)        !on remnants
c        kolp(ip)=kolp(ip)+1
c        kolt(it)=kolt(it)+1
c       else
c        itpr(k)=0
       endif
      enddo
      if(ish.ge.5)write(ifch,*)'ncol:',ncol
     

c --- Calculate Z (written to zzremn)


      do ip=1,maproj
       call CalcZZ(1,ip)
      enddo
      do it=1,matarg
       call CalcZZ(-1,it)
      enddo


c ---  recalculate Zptn


      if(irzptn.eq.1)call recalcZPtn
      

c --- fix all variables


      if(ish.ge.4)write(ifch,*)'fix all variables'
      
      do k=1,koll
        if(itpr(k).ge.2)call ProDiSc(k)
      enddo
      
      do ip=1,maproj
       if(lproj(ip).ne.0)call ProReEx( 1,ip) 
      enddo 
      do it=1,matarg
       if(ltarg(it).ne.0)call ProReEx(-1,it)
      enddo 


      if(isigma.eq.1)then
      if(koll.eq.1)then
        if(itpr(1).ne.0)then
          anintine=anintine+1.
          if(itpr(1).eq.2)then
            anintdiff=anintdiff+1.
            if((iep(1).eq.0.and.iet(1).eq.2).or.
     &           (iet(1).eq.0.and.iep(1).eq.2))anintsdif=anintsdif+1.
          endif
        endif
      else
        aidif=0.
        aiine=0.
        do k=1,koll
          if(aidif.ge.0..and.itpr(k).eq.2)then
            aidif=aidif+iep(k)+iet(k)
          elseif(itpr(k).eq.1)then
            aiine=aiine+1.
            aidif=-0.5
          endif
        enddo
        if(aidif.gt.0.)anintdiff=anintdiff+1.
        if(aidif.eq.2.)anintsdif=anintsdif+1.
        if(aiine+aidif.gt.0.)anintine=anintine+1.
      endif
      endif

      if(ish.ge.6)call XPrint('After fixing:&')


c --- Plot MC pomeron number ---

      if(nemsi.le.4.and.ntry2.eq.0.and.irea.ge.0)then
       if(iemsi1.eq.1)call xEmsI1(1,100,omlog) 
       if(iemsi2.eq.1)call xEmsI2(1,100) 
       if(iemsi1.eq.1.and.ncol.gt.0)call xEmsI1(2,0,omlog)
       if(iemsi2.eq.1.and.ncol.gt.0)call xEmsI2(2,0)
       if((iemsi1.eq.1.or.iemsi2.eq.1).and.ncol.eq.0)nemsi=nemsi-1
      endif
      
      if(ntry2.eq.0)then
        if(iemsb.eq.1)then      ! plot
          do k=1,koll
            if(itpr(k).eq.0)call xEmsB(1,3,k)
            if(itpr(k).eq.1)call xEmsB(1,4,k)
            if(itpr(k).eq.2)call xEmsB(1,5,k)
            if(itpr(k).ne.0)call xEmsB(1,6,k)
          enddo
        endif
      endif



      if(ncol.eq.0)goto 998

c --- Treat Pomerons ---------------------------------------


c --- Check minimum mass ---
 
      do k=1,koll
      do n=1,nprmx(k)
        if(xpr(n,k).lt.cumpom**2/engy**2)then
          nnb=nbkpr(n,k)
          nnv=nvpr(n,k)
          if(nnv.ne.0)then
            nbkpr(nnv,k)=0                  !if bckp Pomeron 
          endif
          if(nnb.ne.0)then
            call VirPom(k,nnb,1)            !if hard backup exist
            nbkpr(n,k)=0                    !remove it
          endif
          call VirPom(k,n,2)
        endif
      enddo
      enddo

c --- Set String End Type and Pt
     
      do k=1,koll
        ip=iproj(k)
        it=itarg(k)
        do n=1,nprmx(k)

          if(idpr(n,k).gt.0)then

          ntry=0
          vpom=.false.
          do i=1,2
            icp(i)=icproj(i,ip)
            ict(i)=ictarg(i,it)
          enddo
 100      ntry=ntry+1
          iret=0
          if(ntry.ge.200)vpom=.true.
          if(ntry.gt.1)then
       if(ish.ge.4)write(ifch,*)'Try again setting string ends for k,n'
     &                               ,k,n,ntry
            do i=1,2
              icproj(i,ip)=icp(i)
              ictarg(i,it)=ict(i)
            enddo
            call RmPt(k,n)
          endif

        call ProSeTy(k,n)
        call ProSePt(k,n)
c      enddo
c      enddo
      
c --- Check Pomeron mass

c      do k=1,koll
c      do n=1,nprmx(k)
       if(idpr(n,k).ne.0.and.ivpr(n,k).ne.0)then
        px=xxp1pr(n,k)+xxp2pr(n,k)+xxm1pr(n,k)+xxm2pr(n,k)
        py=xyp1pr(n,k)+xyp2pr(n,k)+xym1pr(n,k)+xym2pr(n,k)
        pomass=xpr(n,k)*plc*plc-px*px-py*py
        if(pomass.le.0.d0)then
          nnv=nvpr(n,k)
          nnb=nbkpr(n,k)
          idfpom=iabs(idfpr(n,k))
          if(vpom)then
            call VirPom(k,n,3)  !call RmPt(k,n)
            if(nnv.ne.0)then    !bckp Pomeron
              nbkpr(nnv,k)=0
            endif
            if(nnb.ne.0)then    !big Pomeron with bckp one
              ivpr(nnb,k)=1
              nvpr(nnb,k)=0
              idfpr(nnb,k)=idfpom
              npr(1,k)=npr(1,k)+1
              npr(3,k)=npr(3,k)-1
            endif
          else
            goto 100
          endif
        endif
       endif
c      enddo
c      enddo
      
c --- Define String ends for "backup" Pomerons ---

c      do k=1,koll
c      do n=1,nprmx(k)
        if(nvpr(n,k).ne.0)call ProSeX(k,n,iret)
        if(iret.eq.1)then
          if(vpom)then
            nn=nvpr(n,k)
            call VirPom(k,n,7)
            nbkpr(nn,k)=0
          else
            goto 100
          endif
        endif
        iret=0
        iret2=0
c      enddo
c      enddo

c --- Define String ends for "normal" Pomerons ---

c      do k=1,koll
c      do n=1,nprmx(k)
        if(nvpr(n,k).eq.0)call ProSeX(k,n,iret)
        if(iret.eq.1)then
          if(vpom)then
            call VirPom(k,n,12)
          else
            goto 100
          endif
        endif
        iret=0
        iret2=0

      endif

      enddo
      enddo


c --- Write ---

 998  call emszz
      if(ncol.eq.0)goto 1000

      do k=1,koll
       if(itpr(k).eq.1)call emswrpom(k,iproj(k),maproj+itarg(k))
      enddo


c --- Treat hard Pomeron

      do k=1,koll
        do n=1,nprmx(k)
          if(idpr(n,k).eq.3)then
            if(ishpom.eq.1)then
              call psahot(k,n,iret)
              if(iret.eq.1)then
                if(nbkpr(n,k).ne.0)then
                  nn=nbkpr(n,k)
                  call ProSeX(k,nn,iret2)
                  if(iret2.eq.1)then
                    call VirPom(k,nn,7)
                    istptl(nppr(nn,k))=32
                    nbkpr(n,k)=0
                  else
                    ivpr(nn,k)=1
                    nvpr(nn,k)=0
                    idfpr(nn,k)=idfpr(n,k)
                    npr(1,k)=npr(1,k)+1
                    npr(3,k)=npr(3,k)-1
                    ansff=ansff+1 !counters
                    anshf=anshf-1
                  endif
                endif
                call VirPom(k,n,16)
                istptl(nppr(n,k))=32
              elseif(nbkpr(n,k).ne.0)then
                nn=nbkpr(n,k)
                call VirPom(k,nn,17)
                istptl(nppr(nn,k))=32
                nbkpr(n,k)=0
              endif
              iret=0
            else
              istptl(nppr(n,k))=32
              if(nbkpr(n,k).ne.0)then
                nn=nbkpr(n,k)
                istptl(nppr(nn,k))=32
              endif
            endif
          endif
        enddo
      enddo


c --- Treat "normal" soft Pomerons ---

      do k=1,koll
        do n=1,nprmx(k)
          if(nvpr(n,k).eq.0)then
            if(isopom.eq.1)then
              call ProSeF(k,n,iret)
              if(iret.eq.1)then 
                call VirPom(k,n,18)
                istptl(nppr(n,k))=32
              endif
              iret=0
            else
              istptl(nppr(n,k))=32
            endif
          endif
        enddo
      enddo


c --- Treat Remnants -----------------------------------------


c --- Diffractive Pt

      do k=1,koll
        call ProDiPt(k)
      enddo

      do ip=1,maproj
c Here and later "kolp(ip).ne.0" replaced by "iep(ip).ne.-1" to count
c projectile and target nucleons which are counted in paires but are not used
c in collision (no diffractive or inelastic interaction) as slow particles
c at the end. Then we can use them in ProRem to give mass to all other nucleons
c and avoid energy conservation violation that utrescl can not treat 
c (and it gives a reasonnable number of grey particles even if distributions 
c are not really reproduced).
c       if(kolp(ip).ne.0)call ProCop(ip,ip)
       if(iep(ip).ne.-1)call ProCop(ip,ip)
      enddo
      do it=1,matarg
       if(iet(it).ne.-1)call ProCot(it,maproj+it)
c       if(kolt(it).ne.0)call ProCot(it,maproj+it)
      enddo

     
c ---- Remnant Masses (ProReM)


      if(ish.ge.6)call XPrint('Before  ProReM:&')
      ntry=0
      call StoRe(1)             !Store Remnant configuration
 123  ntry=ntry+1
      nshuffi=maproj+matarg
      if(nshuffi.gt.500)
     &call utstop('ems: increase dimension for ishuff&')
      do ip=1,maproj
        ishuff(ip)=ip           !positive for projectile 
      enddo
      do it=1,matarg
        ishuff(maproj+it)=-it   !negative for target
      enddo

      nshuff=maproj+matarg

      do while(nshuff.gt.0)

        if(nshuff.gt.1.and.koll.eq.1.and.iep(1).ne.iet(1))then
c to set the mass of diffractive or not excited remnant first 
c (to avoid unlimited mass of inelastic remants)
          if(iep(1).ne.1.and.iep(1).lt.iet(1))then
            indx=1
          else
            indx=2
          endif
        else
          indx=1+int(rangen()*float(nshuff))
        endif
        if(ishuff(indx).gt.0)then
          ip=ishuff(indx)
c          if(kolp(ip).ne.0)call ProReM( 1,ip,iret) 
          if(iep(ip).ne.-1)call ProReM( 1,ip,iret) 
        else
          it=-ishuff(indx)
c          if(kolt(it).ne.0)call ProReM(-1,it,iret)  
          if(iet(it).ne.-1)call ProReM(-1,it,iret)  
        endif

        if(iret.eq.1)then
          !----------------------------------------
          !If there is a problem, try again shuffle (30 times), 
          !if it doesn't work, for pp, try 10 times with the same type 
          !of event and if doesn't work redo event; 
          !for pA redo event ; and for AB (with A or B >10)
          !continue with some ghosts ...
          !----------------------------------------
          if(ntry.lt.30)then
            if(ish.ge.3)write(ifch,*)'shuffle, try again',ntry         
            call StoRe(-1)         !Restore Remnant configuration
            iret=0
            goto 123          
          elseif(maproj+matarg.eq.2.and.ntry2.lt.10)then
            ntry2=ntry2+1
            if(ish.ge.2)write(ifch,*)'ProRem, try again ! ntry=',ntry2
            if(ish.ge.2)write(ifmt,*)'ProRem, try again ! ntry=',ntry2
            goto 0001
          elseif(maproj.le.10.or.matarg.le.10)then
            if(ish.ge.2)write(ifch,*)'ProRem, redo event ! ntry=',ntry
            if(ish.ge.2)write(ifmt,*)'ProRem, redo event ! ntry=',ntry
            goto 1000
          else
            iret=10
          endif
        endif

        ishuff(indx)=ishuff(nshuff)
        nshuff=nshuff-1

      enddo

      iret=0
      if(ish.ge.6)call XPrint('After ProReM:&')


c --- Write Z into zpaptl for connected strings


      if(isplit.eq.1)then
        do ip=1,maproj
          if(kolp(ip).ne.0)call WriteZZ(1,ip)
        enddo
        do it=1,matarg
          if(kolt(it).ne.0)call WriteZZ(-1,it)
        enddo
      endif


c --- Write Remnants


      do ip=1,maproj
c       if(kolp(ip).ne.0)call emswrp(ip,ip)
       if(iep(ip).ne.-1)call emswrp(ip,ip)
      enddo
      do it=1,matarg
c       if(kolt(it).ne.0)call emswrt(it,maproj+it)
       if(iet(it).ne.-1)call emswrt(it,maproj+it)
      enddo


c --- Remnant Flavors (ProReF)


      do ip=1,maproj
        call ProReF(1,ip)
      enddo 
      do it=1,matarg
        call ProReF(-1,it)
      enddo 
        
999   continue 

c     plot
c     ----

       if(iemspx.eq.1)then
       do ko=1,koll
        if(nprt(ko).gt.0)then
         do np=1,nprmx(ko)
          if(idpr(np,ko).gt.0)then
           call xEmsPx(1,sngl(xpr(np,ko)),sngl(ypr(np,ko)),nprt(ko))
          endif
         enddo
        endif
       enddo
      endif

      if(iemspbx.eq.1)then
       do k=1,koll
        if(nprt(k).gt.0)then
         do n=1,nprmx(k)
          if(idpr(n,k).eq.3)then
            je1=min(1,nemispr(1,n,k))
            je2=min(1,nemispr(2,n,k))
            jex=1+je1+2*je2
            call xEmsP2(1,1+idhpr(n,k),jex
     *            ,sngl(xppr(n,k))
     *            ,sngl(xmpr(n,k))
     *            ,sngl(xpprbor(n,k)),sngl(xmprbor(n,k))
     *            ,ptprboo(1,n,k),ptprboo(2,n,k)  )
          endif
         enddo
        endif
       enddo
      endif

      
      if(iemsse.eq.1)then
       do ko=1,koll
        if(nprt(ko).gt.0)then
         do np=1,nprmx(ko)
          if(idpr(np,ko).gt.0)then
           ptp1=sngl(xxp1pr(np,ko)**2+xyp1pr(np,ko)**2)
           ptp2=sngl(xxp2pr(np,ko)**2+xyp2pr(np,ko)**2)
           ptm1=sngl(xxm1pr(np,ko)**2+xym1pr(np,ko)**2)
           ptm2=sngl(xxm2pr(np,ko)**2+xym2pr(np,ko)**2)
           call xEmsSe(1,sngl(xp1pr(np,ko)),ptp1,1,1)
           call xEmsSe(1,sngl(xp2pr(np,ko)),ptp2,1,1)
           call xEmsSe(1,sngl(xm1pr(np,ko)),ptm1,-1,1)
           call xEmsSe(1,sngl(xm2pr(np,ko)),ptm2,-1,1)
           call xEmsSe(1,sngl(xp1pr(np,ko)),sngl(xm1pr(np,ko)),1,2)
           call xEmsSe(1,sngl(xm2pr(np,ko)),sngl(xp2pr(np,ko)),1,2)
          endif
         enddo
        endif
       enddo
      endif

      if(iemsdr.eq.1)then
       do i=maproj+matarg+1,nptl
        if(istptl(iorptl(i)).eq.41)then
          xpdr=(pptl(4,i)+pptl(3,i))/sngl(plc)
          xmdr=(pptl(4,i)-pptl(3,i))/sngl(plc)
          if(ityptl(i).eq.41)call xEmsDr(1,xpdr,xmdr,1)
          if(ityptl(i).eq.51)call xEmsDr(1,xpdr,xmdr,2)
          if(ityptl(i).eq.42)call xEmsDr(1,xpdr,xmdr,3)
          if(ityptl(i).eq.52)call xEmsDr(1,xpdr,xmdr,4)
        endif
       enddo
      endif

      if(iemsrx.eq.1)then
       do i=1,maproj
        if(kolp(i).gt.0)call xEmsRx(1,1,sngl(xpp(i)),sngl(xmp(i)))
       enddo
       do j=1,matarg
        if(kolt(j).gt.0)call xEmsRx(1,2,sngl(xmt(j)),sngl(xpt(j)))
       enddo
      endif


c     exit
c     ----

 1000 continue
c      write(*,*)'emsaa-iret',iret
      if(ish.ge.1.and.iret.ne.0)write(ifch,*)'iret not 0 (ems)=> redo' 
      call utprix('emsaa ',ish,ishini,4)
      return
      end       


c----------------------------------------------------------------------
      subroutine StoCon(mode,k,n)
c----------------------------------------------------------------------
c store or restore configuration
c   mode = 1 (store) or -1 (restore)
c   k = collision index
c   n = pomeron index
c----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      ip=iproj(k)
      it=itarg(k)

      if(mode.eq.1)then 

       do i=0,3
        nprx0(i)=npr(i,k)
       enddo
       nprtx0=nprt(k)
       idx0=idpr(n,k)
       xxpr0=xpr(n,k)
       yx0=ypr(n,k)
       xxppr0=xppr(n,k)
       xxmpr0=xmpr(n,k)
       nppx0=npp(ip)
       nptx0=npt(it)
       xppx0=xpp(ip)
       xppstx0=xppmx(ip)
       xmpstx0=xppmn(ip)
       xmtx0=xmt(it)
       xptstx0=xmtmx(it)
       xmtstx0=xmtmn(it)
      
      elseif(mode.eq.2)then 

       do i=0,3
        nprx(i)=npr(i,k)
       enddo
       nprtx=nprt(k)
       idx=idpr(n,k)
       xxpr=xpr(n,k)
       yx=ypr(n,k)
       xxppr=xppr(n,k)
       xxmpr=xmpr(n,k)
       nppx=npp(ip)
       nptx=npt(it)
       xppx=xpp(ip)
       xppstx=xppmx(ip)
       xmpstx=xppmn(ip)
       xmtx=xmt(it)
       xptstx=xmtmx(it)
       xmtstx=xmtmn(it)
      
      elseif(mode.eq.-1)then 

       do i=0,3
        npr(i,k)=nprx0(i)
       enddo
       nprt(k)=nprtx0
       idpr(n,k)=idx0
       xpr(n,k)=xxpr0
       ypr(n,k)=yx0
       xppr(n,k)=xxppr0
       xmpr(n,k)=xxmpr0
       npp(ip)=nppx0
       npt(it)=nptx0
       xpp(ip)=xppx0
       xppmx(ip)=xppstx0
       xppmn(ip)=xmpstx0
       xmt(it)=xmtx0
       xmtmx(it)=xptstx0
       xmtmn(it)=xmtstx0
      
      elseif(mode.eq.-2)then 

       do i=0,3
        npr(i,k)=nprx(i)
       enddo
       nprt(k)=nprtx
       idpr(n,k)=idx
       xpr(n,k)=xxpr
       ypr(n,k)=yx
       xppr(n,k)=xxppr
       xmpr(n,k)=xxmpr
       npp(ip)=nppx
       npt(it)=nptx
       xpp(ip)=xppx
       xppmx(ip)=xppstx
       xppmn(ip)=xmpstx
       xmt(it)=xmtx
       xmtmx(it)=xptstx
       xmtmn(it)=xmtstx
      
      else 
      call utstop('mode should integer from -2 to 2 (without 0)&') 
      endif 
      return
      end

c-------------------------------------------------------------------------
      subroutine RemPom(k,n)
c-------------------------------------------------------------------------
c remove pomeron  
c-------------------------------------------------------------------------
      include 'epos.incems'
      include 'epos.inc'

      ip=iproj(k)
      it=itarg(k)
      npr(idpr(n,k),k)=npr(idpr(n,k),k)-1  !nr of pomerons
      nprt(k)=npr(1,k)+npr(3,k)
      if(idpr(n,k).gt.0)then
       npp(ip)=npp(ip)-1                     !nr of pomerons per proj
       npt(it)=npt(it)-1                     !nr of pomerons per targ
       idpr(n,k)=0
       xpp(ip)=xpp(ip)+xppr(n,k)
       xmt(it)=xmt(it)+xmpr(n,k)
       xpr(n,k)=0.d0
       ypr(n,k)=0.d0
       xppr(n,k)=0.d0
       xmpr(n,k)=0.d0
       


      endif

      end

c-------------------------------------------------------------------------
      subroutine ProPo(k,n)
c-------------------------------------------------------------------------
c propose pomeron type = idpr(n,k
c-------------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      double precision wzero,wzerox
      common/cwzero/wzero,wzerox
      
      ip=iproj(k)
      it=itarg(k)

      idpr(n,k)=0

      if(dble(rangen()).gt.wzero)then
        idpr(n,k)=1


c nbr of pomerons per proj
       npp(ip)=npp(ip)+1 
c nbr of pomerons per targ
       npt(it)=npt(it)+1 

      endif

      npr(idpr(n,k),k)=npr(idpr(n,k),k)+1 !nr of pomerons
      nprt(k)=npr(1,k)+npr(3,k)
      

      end


c-------------------------------------------------------------------------
      subroutine ProXY(k,n)
c-------------------------------------------------------------------------
c propose pomeron x,y 
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      double precision xp,xm,om1xprk,om1xmrk,anip,anit,eps
     &,xprem,xmrem
      parameter (eps=1.d-30)
      

      ip=iproj(k)
      it=itarg(k)

      
      xpr(n,k)=0.d0
      ypr(n,k)=0.d0

      
      if(idpr(n,k).ne.0)then
          xprem=xpp(ip)
          xmrem=xmt(it)
c because of fom, it's not symetric any more if we choose always xp first 
c and then xm ... so choose it randomly.
          if(rangen().lt.0.5)then 
            xp=om1xprk(k,xprem,xmrem,1)
            xm=om1xmrk(k,xp,xprem,xmrem,1)
          else
            xm=om1xprk(k,xmrem,xprem,-1)
            xp=om1xmrk(k,xm,xmrem,xprem,-1)
          endif
          xpr(n,k)=xp*xm
          ypr(n,k)=0.d0
          if(xm.gt.eps.and.xp.gt.eps)then
            ypr(n,k)=0.5D0*dlog(xp/xm)
            xppr(n,k)=xp
            xmpr(n,k)=xm
          else
            if(ish.ge.1)write(ifmt,*)'Warning in ProXY ',xp,xm
            npr(idpr(n,k),k)=npr(idpr(n,k),k)-1 
            idpr(n,k)=0
            npr(idpr(n,k),k)=npr(idpr(n,k),k)+1 
            xpr(n,k)=0.d0
            ypr(n,k)=0.d0
            xppr(n,k)=0.d0
            xmpr(n,k)=0.d0
            nprt(k)=npr(1,k)+npr(3,k)
            npp(ip)=npp(ip)-1   !nr of pomerons per proj
            npt(it)=npt(it)-1   !nr of pomerons per targ
          endif
      
c Update xp and xm of remnant, and change the limit to have big enought mass.

        anip=dble(npp(ip))
        anit=dble(npt(it))
        xpp(ip)=xpp(ip)-xppr(n,k)
        xppmn(ip)=min(1.d0,anip*xpmn(ip)/xmpmx(ip))
        xmt(it)=xmt(it)-xmpr(n,k)
        xmtmn(it)=min(1.d0,anit*xtmn(it)/xptmx(it))

      endif

      end
      
c-------------------------------------------------------------------------
      double precision function wmatrix(k,n)
c-------------------------------------------------------------------------
c proposal matrix w(a->b), considering pomeron type, x, y 
c-------------------------------------------------------------------------

      include 'epos.incems'
      include 'epos.incsem'
      double precision wzero,wzerox,Womegak,xprem,xmrem,om1intgck
      common/cwzero/wzero,wzerox


c      ip=iproj(k)
c      it=itarg(k)

      if(idpr(n,k).eq.0)then
        wmatrix=wzero
      else
          xprem=1.d0!xpp(ip)+xppr(n,k)
          xmrem=1.d0!xmt(it)+xmpr(n,k)
          wmatrix=(1d0-wzero)/om1intgck(k,xprem,xmrem)
     *           *Womegak(xppr(n,k),xmpr(n,k),xprem,xmrem,k)
      endif
      
      
      end 

c-------------------------------------------------------------------------
      double precision function omega(n,k)
c-------------------------------------------------------------------------
c calculates partial omega for spot (k,n)
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      common/cwzero/wzero,wzerox
      double precision wzero,wzerox,eps
      parameter(eps=1.d-15)
      double precision PhiExpoK,omGamk,xp,xm,fom
      double precision plc,s
      common/cems5/plc,s
      common/nucl3/phi,bimp
      
      omega=0.d0

      ip=iproj(k)
      it=itarg(k)

      if(xpp(ip).lt.xppmn(ip)+eps.or.xpp(ip).gt.1.d0+eps)goto 1001
      if(xmt(it).lt.xmtmn(it)+eps.or.xmt(it).gt.1.d0+eps)goto 1001 
          
      omega=xpp(ip)**dble(alplea(iclpro))
     &     *xmt(it)**dble(alplea(icltar))
        
c      ztg=0
c      zpj=0
c      nctg=0
c      ncpj=0
c      zsame=nprt(k)
c      if(idpr(n,k).gt.0)then
c        if(nprt(k).le.0)stop'omega: nprt(k) should be positive !!!!    '
c        zsame=zsame-1 
c      endif
c      nlpop=nint(zsame)
c      nlpot=nint(zsame)
c      bglaub2=sigine/10./pi        !10= fm^2 -> mb
c      bglaub=sqrt(bglaub2)
c      b2x=epscrp*epscrp*bglaub2  
c      b2=bk(k)**2
c      ztgx=epscrw*exp(-b2/2./b2x)*fscra(engy/egyscr)   
c      zpjx=epscrw*exp(-b2/2./b2x)*fscra(engy/egyscr)  
c         
c      if(koll.gt.1)then
c        do li=1,lproj(ip)
c          kk=kproj(ip,li)
c          if(kk.ne.k)then
c            b2=bk(kk)**2
c            if(b2.le.bglaub2)nctg=nctg+1
c            ztg=ztg+epscrw*exp(-b2/2./b2x)*fscro(engy/egyscr) 
c            nlpop=nlpop+nprt(kk)
c          endif
c        enddo
c        do li=1,ltarg(it)
c          kk=ktarg(it,li)
c          if(kk.ne.k)then
c            b2=bk(kk)**2
c            if(b2.le.bglaub2)ncpj=ncpj+1
c            zpj=zpj+epscrw*exp(-b2/2./b2x)*fscro(engy/egyscr) 
c            nlpot=nlpot+nprt(kk)
c          endif
c        enddo
c      endif
      !  zpjx+zpj is equal to zparpro(k)
      !  ztgx+ztg is equal to zpartar(k)
      zprj=zparpro(k)  !zsame+zpj 
      ztgt=zpartar(k)  !zsame+ztg
      if(idpr(n,k).eq.0)then
        omega=omega*wzerox
      else
        xp=xppr(n,k)
        xm=xmpr(n,k)
c        !-------------------------------------------------------------------------
c        ! fom : part of Phi regularization; Phi -> Phi^(n) (n = number of Poms)
c        ! Phi^(0) relevant for Xsect unchanged, apart of (maybe) normalization (Z)
c        !-------------------------------------------------------------------------
        omega=omega*omGamk(k,xp,xm)*gammaV(k)*fom(zprj,xm,bk(k))
     &                                       *fom(ztgt,xp,bk(k))
      endif

      omega=omega*PhiExpoK(k,xpp(ip),xmt(it))


      if(omega.le.0.d0)goto 1001

      if(koll.gt.1)then
        do li=1,lproj(ip)
          kk=kproj(ip,li)
          if(itarg(kk).ne.it)then
            ipl=iproj(kk)
            itl=itarg(kk)
            omega=omega*PhiExpoK(kk,xpp(ipl),xmt(itl))
            if(omega.le.0.d0)goto 1001
          endif
        enddo
        do li=1,ltarg(it)
          kk=ktarg(it,li)
          if(iproj(kk).ne.ip)then
            ipl=iproj(kk)
            itl=itarg(kk)
            omega=omega*PhiExpoK(kk,xpp(ipl),xmt(itl))
            if(omega.le.0.d0)goto 1001
          endif
        enddo
      endif
      
      if(omega.lt.1.d-100)then
        if(ish.ge.6)write(*,*)'omega-exit',omega
        omega=0.d0
      elseif(omega.gt.1.d100)then
        if(ish.ge.6)write(*,*)'omega-exit',omega
        omega=0.d0
      endif
      
      return
      
 1001 continue

      omega=0.d0      
      return
      
      end

c-------------------------------------------------------------------------
      double precision function fom(z,x,b)
c-------------------------------------------------------------------------
      include 'epos.inc'
      double precision x,u,w,z0
      !----------------------------------------------------------------
      ! part of Phi regularization; Phi -> Phi^(n) (n = number of Poms)
      ! Phi^(0) relevant for Xsect unchanged
      !----------------------------------------------------------------
      fom=1d0
      if(z.gt.0..and.alpfomi.gt.0.)then
       z0=dble(alpfomi)
       u=dble(z**gamfom) 
c       u=z0*dble(z/z0)**2. 
       w=u/z0*exp(-dble(b*b/delD(1,iclpro,icltar)))
c       w=10.d0*u
       !---------------------------------------------------
       !e=exp(-0.05*u)  !analytic function with e(0)=1 
       !fom=((1-u)+(u+w)*sqrt(x**2+((u-1+e)/(u+w))**2)) 
       !     fom(z=0)=1  fom(x=0)=e  fom(x=1)~w
       !---------------------------------------------------
       fom=1.d0+w*x**betfom
       !---------------------------------------------------
      endif
      end

c-------------------------------------------------------------------------
      subroutine ProPoTy(k,n)  
c-------------------------------------------------------------------------
c propose pomeron type 
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      double precision ww,w0,w1,w2,w3,w4,w5,w(0:7),aks,eps
     *,xh,yp!,xp,xm
      parameter(eps=1.d-10)
      logical cont
      dimension nnn(3),kkk(3)

      if(idpr(n,k).eq.0)return
      if(ish.ge.4)write(ifch,*)'ProPoTy:k,n,idpr,x',k,n,idpr(n,k)
     *                                              ,xpr(n,k)
      if(idpr(n,k).ne.1)call utstop('ProPoTy: should not happen&')
            
      cont=.true.
      do i=1,3
        nnn(i)=0
        kkk(i)=0
      enddo

      idfpr(n,k)=1
      ip=iproj(k)
      it=itarg(k)
      xh=xpr(n,k)
      yp=ypr(n,k)
c      xp=xppr(n,k)
c      xm=xmpr(n,k)
      nnn(3)=n
      kkk(3)=k
      
      if(iep(ip).ne.-1)iep(ip)=-1
      if(iet(it).ne.-1)iet(it)=-1
    

      idpr(n,k)=1
      
        w0=0.d0
        w1=0.d0
        w2=0.d0
        w3=0.d0
        w4=0.d0
        w5=0.d0

        call WomTy(w,xh,yp,k)


        if(w(0).gt.0.d0)w0=w(0)
        if(w(1).gt.0.d0)w1=w(1)
        if(w(2).gt.0.d0)w2=w(2)
        if(w(3).gt.0.d0)w3=w(3)
        if(w(4).gt.0.d0)w4=w(4)
        if(w(5).gt.0.d0)w5=w(5)

        ww=w0+w1+w2+w3+w4+w5
        if(ish.ge.4)write(ifch,*)'ProPoTy:ww,ww_i'
     *       ,ww,w0/ww*100.d0,w1/ww*100.d0,w2/ww*100.d0
     *       ,w3/ww*100.d0,w4/ww*100.d0,w5/ww*100.d0


        aks=dble(rangen())*ww

        if(ww.lt.eps.or.aks.le.w0)then            !soft pomeron

          if(ish.ge.8)write(ifch,*)'ProPoTy:idpr',idpr(n,k)

        elseif(aks.ge.ww-w5)then !diffractive interaction
          itpr(k)=itpr(k)+2
          if(ish.ge.8)write(ifch,*)'ProPoTy:itpr',itpr(k)
          call RemPom(k,n)
          npr(0,k)=npr(0,k)+1 !nr of empty cells
 
        else

          idpr(n,k)=3
          if(ish.ge.8)write(ifch,*)'ProPoTy:idpr',idpr(n,k)
          npr(3,k)=npr(3,k)+1 
          npr(1,k)=npr(1,k)-1
          bhpr(n,k)=bk(k)
       
          aks=aks-w0
          if(aks.le.w1)then                             !gg-pomeron
            idhpr(n,k)=0
          elseif(aks.le.w1+w2)then                      !qg-pomeron
            idhpr(n,k)=1
          elseif(aks.le.w1+w2+w3)then                   !gq-pomeron
            idhpr(n,k)=2
          elseif(aks.le.w1+w2+w3+w4)then                !qq-pomeron
            idhpr(n,k)=3
          else
            call utstop('ems-unknown pomeron&')
          endif
        
        endif

        if(xpr(n,k).gt.0.d0)then
          antot=antot+1
          antotf=antotf+1
          if(idpr(n,k).eq.1)then
            ansf=ansf+1
            ansff=ansff+1
          endif
          if(idpr(n,k).eq.3)then
            ansh=ansh+1
            anshf=anshf+1
          endif
        endif

      do i=3,1,-1

        if(nnn(i).ne.0.and.kkk(i).ne.0.and.cont)then

          if(idpr(nnn(i),kkk(i)).eq.3)then  

                       !Backup soft Pomeron if sh not possible later

            kb=kkk(i)
            nb=nnn(i)
            ip=iproj(kb)
            it=itarg(kb)        
            do nn=1,nprmx(kb)
              if(idpr(nn,kb).eq.0)then !empty spot
                nbkpr(nb,kb)=nn
                nvpr(nn,kb)=nb
                idpr(nn,kb)=1
                ivpr(nn,kb)=2
                xpr(nn,kb)=xpr(nb,kb)
                ypr(nn,kb)=ypr(nb,kb)
                xppr(nn,kb)=xppr(nb,kb)
                xmpr(nn,kb)=xmpr(nb,kb)
                idfpr(nn,kb)=-idfpr(nb,kb)
                bhpr(nn,kb)=bhpr(nb,kb)
                idp1pr(nn,kb)=0
                idp2pr(nn,kb)=0
                idm1pr(nn,kb)=0
                idm2pr(nn,kb)=0
                xm1pr(nn,kb)=0.d0
                xp1pr(nn,kb)=0.d0
                xm2pr(nn,kb)=0.d0
                xp2pr(nn,kb)=0.d0
                xxm1pr(nn,kb)=0.d0
                xym1pr(nn,kb)=0.d0
                xxp1pr(nn,kb)=0.d0
                xyp1pr(nn,kb)=0.d0
                xxm2pr(nn,kb)=0.d0
                xym2pr(nn,kb)=0.d0
                xxp2pr(nn,kb)=0.d0
                xyp2pr(nn,kb)=0.d0
                goto 10
              endif
            enddo
      if(ish.ge.2)write(ifmt,*)'no empty lattice site, backup lost'
            
 10         continue
          endif
        endif
      enddo

      return
      end   

c-------------------------------------------------------------------------
      subroutine ProDiSc(k)
c-------------------------------------------------------------------------
c propose diffractive scattering
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incsem'
      include 'epos.incems'
      common/col3/ncol,kolpt
      
      ncol=ncol+1
      ip=iproj(k)
      it=itarg(k)
      kolp(ip)=kolp(ip)+itpr(k)/2   !number of cut on remnants
      kolt(it)=kolt(it)+itpr(k)/2
      itpr(k)=2

      
      end
     
c-------------------------------------------------------------------------
      subroutine ProReEx(ir,ii)
c-------------------------------------------------------------------------
c propose remnant excitation 
c for proj (iep) if ir=1 or target (iet) if ir=-1: 
c 0 = no,  1 = inel excitation,  2 = diffr excitation
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incsem'
      include 'epos.incems'
      common/cncl3/iactpr(mamx),iacttg(mamx)

       
      if(ir.eq.1)then                   !proj

        ip=ii
        if(iactpr(ip).eq.0)stop'ProReEx: should not happen (1)'
        mine=0
        mdif=0
        do l=1,lproj(ip)
          kp=kproj(ip,l)
          if(itpr(kp).eq.1)mine=1
          if(itpr(kp).eq.2)mdif=1
        enddo  
        iep(ip)=0
        r=rangen()
        if(mine.eq.1)then   !inelastic
c          if(r.lt.1.-(1.-rexndi(iclpro))**kolp(ip))iep(ip)=1   
          if(r.lt.1.-(1.-rexndi(iclpro)))iep(ip)=1   
        elseif(mdif.eq.1)then        !diffr
          if(r.lt.1.-(1.-rexdif(iclpro))**kolp(ip))iep(ip)=2   
c          if(r.lt.1.-(1.-rexdif(iclpro)))iep(ip)=2   
        endif
      elseif(ir.eq.-1)then                !targ
      
        it=ii
        if(iacttg(it).eq.0)stop'ProReEx: should not happen (2)'
        mine=0
        mdif=0
        do l=1,ltarg(it)
          kt=ktarg(it,l)
          if(itpr(kt).eq.1)mine=1
          if(itpr(kt).eq.2)mdif=1
        enddo  
        iet(it)=0
        r=rangen()
        if(mine.eq.1)then   !inelastic
c          if(r.lt.1.-(1.-rexndi(icltar))**kolt(it))iet(it)=1   
          if(r.lt.1.-(1.-rexndi(icltar)))iet(it)=1   
        elseif(mdif.eq.1)then        !diffr
          if(r.lt.1.-(1.-rexdif(icltar))**kolt(it))iet(it)=2   
c          if(r.lt.1.-(1.-rexdif(icltar)))iet(it)=2   
        endif
        
      endif

      end
      

c-------------------------------------------------------------------------
      subroutine ProDiPt(k)
c-------------------------------------------------------------------------
c propose transverse momentum for diffractive interaction
c-------------------------------------------------------------------------

      include 'epos.incems'
      include 'epos.inc'
      double precision xxe,xye,p5sqpr,p5sqtg
      double precision plc,s
      common/cems5/plc,s

      ip=iproj(k)
      it=itarg(k)
      

c generate p_t for diffractive 

       if(ptdiff.ne.0.)then
         if(itpr(k).eq.2)then
           ptd=ptdiff  
c           ad=pi/4./ptd**2
c           r=rangen()
           pt=ranpt()*ptd !sqrt(-alog(r)/ad)
         elseif(itpr(k).eq.0)then   !pt for non-wounded nucleon (usefull in ProRem to avoid problem in utrescl)
           ptnw=0.005
           pt=ranptd()*ptnw
         else
           xxe=0d0
           xye=0d0
           goto 10
         endif
         phi=2.*pi*rangen()
         xxe=dble(pt*cos(phi))
         xye=dble(pt*sin(phi))
       else
         xxe=0d0
         xye=0d0
       endif

c update remnant p_t 

 10    xxp(ip)=xxp(ip)-xxe
       xyp(ip)=xyp(ip)-xye
       xxt(it)=xxt(it)+xxe
       xyt(it)=xyt(it)+xye
      
       if(iep(ip).eq.6.and.iet(it).eq.6)then             !to simulate the fact that originally 
         if(itpr(k).eq.3)then
           call StoCon(-k,k,1)  !to fixe mass of corresponding remnants
           xpp(ip)=xpp(ip)-xppr(1,k)
           xpt(it)=xpt(it)+xppr(1,k)
           xmt(it)=xmt(it)-xmpr(1,k)
           xmp(ip)=xmp(ip)+xmpr(1,k)           
           idpr(1,k)=0
           xpr(1,k)=0.d0
           ypr(1,k)=0.d0
           xppr(1,k)=0.d0
           xmpr(1,k)=0.d0
         endif
         p5sqpr=xpp(ip)*plc*xmp(ip)*plc-dble(amproj*amproj)
         p5sqtg=xpt(it)*plc*xmt(it)*plc-dble(amtarg*amtarg)
         phi=2.*pi*rangen()
         ntry=0
 20      ntry=ntry+1
         pt=ranptcut(ptsemx)*ptsend**2     
         if(ntry.lt.100.and.(p5sqpr-dble(pt*pt).lt.0.d0
     &                   .or.p5sqtg-dble(pt*pt).lt.0.d0))then
             goto 20
         else
           pt=ranptcut(ptsemx)*ptsendi 
         endif
         xxe=dble(pt*cos(phi))
         xye=dble(pt*sin(phi))
         xxp(ip)=xxp(ip)-xxe
         xyp(ip)=xyp(ip)-xye
         xxt(it)=xxt(it)+xxe
         xyt(it)=xyt(it)+xye
       endif

       
       end
      
c-------------------------------------------------------------------------
      subroutine ProSePt(k,n)
c-------------------------------------------------------------------------
c propose transverse momentum for string ends
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      if(ivpr(n,k).eq.2)return            !Backup Pomeron

      ip=iproj(k)
      it=itarg(k)
      amk0=(qmass(1)+qmass(2)+qmass(3))     !mass for mt distribution but spoil <pt> of anti-proton at low energy 

ctp060829      id=idpr(n,k) 
ctp060829      ih=0
ctp060829      if(id.eq.3)ih=1
      
c generate p_t for string ends  (proj)

c      nph=0
c      do l=1,lproj(ip)
c        kk=kproj(ip,l)
c        nph=nph+npr(3,kk)        
c      enddo
c      
c      !---proj-----
c        zz=0
c        if(isplit.eq.1)then
c         if(lproj(ip).ge.1)then
c          do l=1,lproj(ip)
c           kpair=kproj(ip,l)
c           if(itpr(kpair).eq.1)then
c            zz=zz+zparpro(kpair)
c           endif
c          enddo 
c         endif
c        endif  
c      !------   
        ptsef=ptsemx
c        if(iep(ip).eq.0)ptsef=ptsef*ptsendi
        ptsendx=ptsend
        if(iep(ip).eq.0)ptsendx=ptsendi
        ptsendy = ptsendi 
        
      if(idp1pr(n,k).gt.0)then
         if(idp1pr(n,k).eq.4.or.idp1pr(n,k).eq.5)then   !diquarks
c          pt=ranptd()*ptsendx
c           if(iep(ip).eq.0)then
c            pt=ranpt()*ptsendy
c          else
c            pt=ranptd()*ptsendy
c          endif
           pt=ranptcut(ptsef)*ptsendy
           amk1=amk0!+qmass(0)*0.5 !mass for mt distribution with bounding energy for diquark
         else
c          pt=ranptd()*ptsendx
c           if(iep(ip).eq.0)then
c             pt=ranpt()*ptsendx
c           else
c             pt=ranptd()*ptsendx
c           endif
          pt=ranptcut(ptsef)*ptsendx
           amk1=amk0
         endif
         pt=sqrt(pt*pt+2.*pt*amk1) !sample mt-m0 instead of pt ...
         phi=2.*pi*rangen()
         xxp1pr(n,k)=dble(pt*cos(phi))
         xyp1pr(n,k)=dble(pt*sin(phi))
      else
         xxp1pr(n,k)=0d0
         xyp1pr(n,k)=0d0
      endif
      if(idp2pr(n,k).gt.0)then
         if(idp2pr(n,k).eq.4.or.idp2pr(n,k).eq.5)then
c           pt=ranptd()*ptsendy     
c           if(iep(ip).eq.0)then
c             pt=ranpt()*ptsendy
c           else
c             pt=ranptd()*ptsendy
c           endif
           pt=ranptcut(ptsef)*ptsendy
           amk1=amk0!+qmass(0)*0.5 !mass for mt distribution with bounding energy for diquark
         else
c           pt=ranptd()*ptsendx 
c           if(iep(ip).eq.0)then
c             pt=ranpt()*ptsendx
c           else
c             pt=ranptd()*ptsendx
c           endif
           pt=ranptcut(ptsef)*ptsendx
           amk1=amk0
         endif
         pt=sqrt(pt*pt+2.*pt*amk1) !sample mt-m0 instead of pt ...
         phi=2.*pi*rangen()
         xxp2pr(n,k)=dble(pt*cos(phi))
         xyp2pr(n,k)=dble(pt*sin(phi))
      else
         xxp2pr(n,k)=0d0
         xyp2pr(n,k)=0d0
      endif
c generate p_t for string ends  (targ)


c      nph=0
c      do l=1,ltarg(it)
c        kk=ktarg(it,l)
c        nph=nph+npr(3,kk)        
c      enddo
c
c      !---targ-----
c        zz=0
c        if(isplit.eq.1)then
c         if(ltarg(it).ge.1)then
c          do l=1,ltarg(it)
c           kpair=ktarg(it,l)
c           if(itpr(kpair).eq.1)then
c            zz=zz+zpartar(kpair)
c           endif
c          enddo 
c         endif  
c        endif
c      !---------
        ptsef=ptsemx
c        if(iet(it).eq.0)ptsef=ptsef*ptsendi
        ptsendx=ptsend
        if(iet(it).eq.0)ptsendx=ptsendi
        ptsendy = ptsendx 

      if(idm1pr(n,k).gt.0)then
         if(idm1pr(n,k).eq.4.or.idm1pr(n,k).eq.5)then
c           pt=ranptd()*ptsendy  
c           if(iet(it).eq.0)then
c             pt=ranpt()*ptsendy
c           else
c             pt=ranptd()*ptsendy
c           endif
           pt=ranptcut(ptsef)*ptsendy
           amk1=amk0!+qmass(0)*0.5 !mass for mt distribution with bounding energy for diquark
         else
c           pt=ranptd()*ptsendx 
c           if(iet(it).eq.0)then
c             pt=ranpt()*ptsendx
c           else
c             pt=ranptd()*ptsendx
c           endif
           pt=ranptcut(ptsef)*ptsendx
           amk1=amk0
         endif
         pt=sqrt(pt*pt+2.*pt*amk1) !sample mt-m0 instead of pt ...
         phi=2.*pi*rangen()
         xxm1pr(n,k)=dble(pt*cos(phi))
         xym1pr(n,k)=dble(pt*sin(phi))
      else
         xxm1pr(n,k)=0d0
         xym1pr(n,k)=0d0
      endif
      if(idm2pr(n,k).gt.0)then
         if(idm2pr(n,k).eq.4.or.idm2pr(n,k).eq.5)then
c           pt=ranptd()*ptsendy     
c           if(iet(it).eq.0)then
c             pt=ranpt()*ptsendy
c           else
c             pt=ranptd()*ptsendy
c           endif
           pt=ranptcut(ptsef)*ptsendy
           amk1=amk0!+qmass(0)*0.5 !mass for mt distribution with bounding energy for diquark
         else
c           pt=ranptd()*ptsendx 
c           if(iet(it).eq.0)then
c             pt=ranpt()*ptsendx
c           else
c             pt=ranptd()*ptsendx
c           endif
           pt=ranptcut(ptsef)*ptsendx
           amk1=amk0
         endif
         pt=sqrt(pt*pt+2.*pt*amk1) !sample mt-m0 instead of pt ...
         phi=2.*pi*rangen()
         xxm2pr(n,k)=dble(pt*cos(phi))
         xym2pr(n,k)=dble(pt*sin(phi))
      else
         xxm2pr(n,k)=0d0
         xym2pr(n,k)=0d0
      endif

c update backup soft pomeron p_t if exist

        if(nbkpr(n,k).ne.0)then
          nn=nbkpr(n,k)
          xxp1pr(nn,k)=xxp1pr(n,k)
          xyp1pr(nn,k)=xyp1pr(n,k)
          xxp2pr(nn,k)=xxp2pr(n,k)
          xyp2pr(nn,k)=xyp2pr(n,k)
          xxm1pr(nn,k)=xxm1pr(n,k)
          xym1pr(nn,k)=xym1pr(n,k)
          xxm2pr(nn,k)=xxm2pr(n,k)
          xym2pr(nn,k)=xym2pr(n,k)
        endif



c update remnant p_t (pomeron)
        xxp(ip)=xxp(ip)-xxp1pr(n,k)-xxp2pr(n,k)
        xyp(ip)=xyp(ip)-xyp1pr(n,k)-xyp2pr(n,k)
        xxt(it)=xxt(it)-xxm1pr(n,k)-xxm2pr(n,k)
        xyt(it)=xyt(it)-xym1pr(n,k)-xym2pr(n,k)
        
        if(ish.ge.6)then
        write(ifch,*) 'ProSePt'
        write(ifch,'(4i14/4d14.3/4d14.3/)')
     * idp1pr(n,k),idp2pr(n,k),idm1pr(n,k),idm2pr(n,k)
     *,xxp1pr(n,k),xxp2pr(n,k),xxm1pr(n,k),xxm2pr(n,k)
     *,xyp1pr(n,k),xyp2pr(n,k),xym1pr(n,k),xym2pr(n,k)
        endif
     
        end

c-----------------------------------------------------------------------
      subroutine ProSeX(k,n,iret)
c-----------------------------------------------------------------------
c calculates x of string ends
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/cems5/plc,s
      double precision s,plc
      common/cems10/a(0:ntypmx),b(0:ntypmx),d(0:ntypmx)
      double precision a,b,d
     *,xp,xm,ap1,ap2,am1,am2,aamin1,aamin2,u
     *,xmn1,xmn2
     
      iret=0
     
      if(itpr(k).ne.1)return
      if(idpr(n,k).ne.1.or.ivpr(n,k).eq.0)return
      
      if(idp1pr(n,k).eq.0.and.idp2pr(n,k).eq.0
     * .and.idm1pr(n,k).eq.0.and.idm2pr(n,k).eq.0)
     *call utstop('no Pomeron in ProSex&')
      
      xp=xppr(n,k)
      xm=xmpr(n,k)
      ap1=a(idp1pr(n,k))
      ap2=a(idp2pr(n,k))
      am1=a(idm1pr(n,k))
      am2=a(idm2pr(n,k))
      aamin1=ammn(idp1pr(n,k)+idm2pr(n,k))
      aamin2=ammn(idp2pr(n,k)+idm1pr(n,k))
      xmn1=(aamin1**2+(xxp1pr(n,k)+xxm2pr(n,k))**2
     &               +(xyp1pr(n,k)+xym2pr(n,k))**2)/s
      xmn2=(aamin2**2+(xxp2pr(n,k)+xxm1pr(n,k))**2
     &               +(xyp2pr(n,k)+xym1pr(n,k))**2)/s
    
      ntry=0
 999  ntry=ntry+1
      if(ntry.gt.100)then
        iret=1
        if(ish.ge.5)write(ifch,*)'Problem in ProSex(k,n)',k,n
        return
      endif        
           
    1 u=dble(rangen())**(1d0/(1d0+ap1)) 
      if(dble(rangen()).gt.(1d0-u)**ap2)goto1 
      xp1pr(n,k)=u*xp
      xp2pr(n,k)=(1-u)*xp
    2 u=dble(rangen())**(1d0/(1d0+am1)) 
      if(dble(rangen()).gt.(1d0-u)**am2)goto2 
      xm1pr(n,k)=u*xm
      xm2pr(n,k)=(1-u)*xm
      
      if(xp1pr(n,k)*xm2pr(n,k).lt.xmn1)then
      goto 999
c       fc=xp1pr(n,k)*xm2pr(n,k)/xmn1   !avoid virpom
c       if(fc.eq.0.)goto 999
c       xp1pr(n,k)=xp1pr(n,k)/sqrt(fc)
c       xm2pr(n,k)=xm2pr(n,k)/sqrt(fc)
      endif
      if(xp2pr(n,k)*xm1pr(n,k).lt.xmn2)then
      goto 999
c       fc=xp2pr(n,k)*xm1pr(n,k)/xmn2   !avoid virpom
c       if(fc.eq.0.)goto 999
c       xp2pr(n,k)=xp2pr(n,k)/sqrt(fc)
c       xm1pr(n,k)=xm1pr(n,k)/sqrt(fc)
      endif
      
      if(ish.ge.6)then
       write(ifch,*) 'ProSeX'
       write(ifch,'(2d28.3,i8)') xp,xm,ntry
       write(ifch,'(4d14.3)')xp1pr(n,k),xp2pr(n,k),xm1pr(n,k),xm2pr(n,k)
       write(ifch,'(4d14.3/)')xp1pr(n,k)*xm2pr(n,k)
     *                   ,xp2pr(n,k)*xm1pr(n,k),  xmn1, xmn2
      endif

      end
c-------------------------------------------------------------------------
      subroutine RmPt(k,n)
c-------------------------------------------------------------------------
c remove pt from pomeron
c-------------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      ip=iproj(k)
      it=itarg(k)
      xxp(ip)=xxp(ip)+xxp1pr(n,k)+xxp2pr(n,k)
      xyp(ip)=xyp(ip)+xyp1pr(n,k)+xyp2pr(n,k)
      xxt(it)=xxt(it)+xxm1pr(n,k)+xxm2pr(n,k)
      xyt(it)=xyt(it)+xym1pr(n,k)+xym2pr(n,k)
      xp1pr(n,k)=0d0
      xp2pr(n,k)=0d0
      xm1pr(n,k)=0d0
      xm2pr(n,k)=0d0
      xxm1pr(n,k)=0d0
      xym1pr(n,k)=0d0
      xxp1pr(n,k)=0d0
      xyp1pr(n,k)=0d0
      xxm2pr(n,k)=0d0
      xym2pr(n,k)=0d0
      xxp2pr(n,k)=0d0
      xyp2pr(n,k)=0d0
      idp1pr(n,k)=0
      idm2pr(n,k)=0
      idp2pr(n,k)=0
      idm1pr(n,k)=0
      end
      
c-------------------------------------------------------------------------
      subroutine VirPom(k,n,id)
c-------------------------------------------------------------------------
c create virtual pomeron
c virtual pomeron: ivpr(n,k)=0, otherwise ivpr(n,k)=1
c-------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/col3/ncol,kolpt
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision plc,s
      common/cems5/plc,s
c      data nvir/0/
c      save nvir
      
      call utpri('VirPom',ish,ishini,3)

      if(idpr(n,k).eq.0)return

      ip=iproj(k)
      it=itarg(k)

      nnv=nvpr(n,k)
      nnb=nbkpr(n,k)

c                        nvir=nvir+1
c                   print *,'  ',id,'   ',nvir

      if(ish.ge.3)then
      write(ifch,*)"virpom ",id," (n,k)",n,k,nnb,nnv,nppr(n,k)
      if(ish.ge.5)write(ifch,*)"remnant in",xpp(ip),xmt(it)
      endif

      if(nnv.ne.0)then
        nn=nnv
        kk=k
        if(idpr(nn,kk).eq.0)then
          nvpr(n,k)=0
        endif
      endif

      if(nnb.ne.0)then
        nn=nnb
        kk=k
        if(idpr(nn,kk).eq.0)then
          nbkpr(n,k)=0
         endif
      endif

 
      if(nbkpr(n,k).eq.0.and.nvpr(n,k).eq.0)then     !normal Pomeron

      npr(0,k)=npr(0,k)+1
      npp(ip)=npp(ip)-1  
      npt(it)=npt(it)-1
      npr(idpr(n,k),k)=npr(idpr(n,k),k)-1 
      nprt(k)=npr(1,k)+npr(3,k)
      antotf=antotf-1
      if(idpr(n,k).eq.1)ansff=ansff-1
      if(idpr(n,k).eq.3)anshf=anshf-1
      kolp(ip)=kolp(ip)-1  
      kolt(it)=kolt(it)-1
      xpp(ip)=xpp(ip)+xppr(n,k)
      xmt(it)=xmt(it)+xmpr(n,k)
      xxp(ip)=xxp(ip)+xxp1pr(n,k)+xxp2pr(n,k)
      xyp(ip)=xyp(ip)+xyp1pr(n,k)+xyp2pr(n,k)
      xxt(it)=xxt(it)+xxm1pr(n,k)+xxm2pr(n,k)
      xyt(it)=xyt(it)+xym1pr(n,k)+xym2pr(n,k)

      if(itpr(k).eq.1.and.nprt(k).eq.0)then  !no more Pomeron on this pair
        if(kolp(ip).eq.0)then
          kolp(ip)=1        !excite nucleon (remnant get pt in ProDiPt) 
          iep(ip)=6         !with inel mass and inverted string 
        endif
        if(kolt(it).eq.0)then
          kolt(it)=1
          iet(it)=6    
        endif
        if(koll.le.2.and.iep(ip).eq.6.and.iet(it).eq.6)then  !for small systems we can store lost informations
          itpr(k)=3             !to use it to define remnant mass
          call StoCon(k,k,n)    !store information on lost Pomeron
        endif
      endif  

      endif


      ivpr(n,k)=0
      nbkpr(n,k)=0
      nvpr(n,k)=0
      idpr(n,k)=0
      idfpr(n,k)=0                     
      xpr(n,k)=0d0
      ypr(n,k)=0d0
      xppr(n,k)=0d0
      xmpr(n,k)=0d0
      idp1pr(n,k)=0
      idp2pr(n,k)=0
      idm1pr(n,k)=0
      idm2pr(n,k)=0
      xm1pr(n,k)=0d0
      xp1pr(n,k)=0d0
      xm2pr(n,k)=0d0
      xp2pr(n,k)=0d0
      xxm1pr(n,k)=0d0
      xym1pr(n,k)=0d0
      xxp1pr(n,k)=0d0
      xyp1pr(n,k)=0d0
      xxm2pr(n,k)=0d0
      xym2pr(n,k)=0d0
      xxp2pr(n,k)=0d0
      xyp2pr(n,k)=0d0

       if(ish.ge.5)write(ifch,*)"remnant out",xpp(ip),xmt(it)
     
      call utprix('VirPom',ish,ishini,3)

      end   

c-----------------------------------------------------------------------
      subroutine StoRe(imod)
c-----------------------------------------------------------------------
c Store Remnant configuration (imod=1) before shuffle  to restore the 
c initial configuration (imod=-1) in case of problem.
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      if(imod.eq.1)then

c       initialize projectile

        do i=1,maproj
          xppst(i)=xpp(i)
          xmpst(i)=xmp(i)
          xposst(i)=xpos(i)
        enddo

c       initialize target

        do j=1,matarg
          xmtst(j)=xmt(j)
          xptst(j)=xpt(j)
          xtosst(j)=xtos(j)
        enddo
        
      elseif(imod.eq.-1)then

c       restore projectile

        do i=1,maproj
          xpp(i)=xppst(i)
          xmp(i)=xmpst(i)
          xpos(i)=xposst(i)
        enddo

c       restore target

        do j=1,matarg
          xmt(j)=xmtst(j)
          xpt(j)=xptst(j)
          xtos(j)=xtosst(j)
        enddo
        
      else
        
        call utstop('Do not know what to do in StoRe.&')

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine CalcZZ(ir,m)
c-----------------------------------------------------------------------
C Calculates zz for remnant m for proj (ir=1) or target (ir=-1)
c   writes it to zzremn(m, 1 or 2)
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      if(ir.eq.1)then
        if(kolp(m).eq.0)then
          zzremn(m,1)=0
          return
        endif    
      elseif(ir.eq.-1)then 
        if(kolt(m).eq.0)then
          zzremn(m,2)=0
          return
        endif
      endif 
      if(isplit.eq.1)then
        if(ir.eq.1)then
          zz=0
          if(lproj(m).ge.1)then
           do l=1,lproj(m)
            kpair=kproj(m,l)
            if(itpr(kpair).eq.1)then
             zz=zz+zparpro(kpair)
            endif
           enddo 
          endif  
          zzremn(m,1)=zz
        elseif(ir.eq.-1)then
          zz=0
          if(ltarg(m).ge.1)then
           do l=1,ltarg(m)
            kpair=ktarg(m,l)
            if(itpr(kpair).eq.1)then
             zz=zz+zpartar(kpair)
            endif
           enddo 
          endif  
          zzremn(m,2)=zz
        else
          stop'CalcZZ: invalid option.          '        
        endif
      else
        if(ir.eq.1) zzremn(m,1)=0
        if(ir.eq.-1)zzremn(m,2)=0
      endif
      end

c-----------------------------------------------------------------------
      subroutine WriteZZ(ir,irem)
c-----------------------------------------------------------------------
c Write Z into zpaptl(K) for connected strings
c                 K is the index for the string end
c                 on the corresponding remnant side
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)

      if(ir.eq.1)then
        jrem=1
      elseif(ir.eq.-1)then
        jrem=2
      endif

      do li=1,lremn(irem,jrem)
        kkk=kremn(irem,li,jrem)
         do n=1,nprmx(kkk)
          if(idpr(n,kkk).ne.0)then
            npom=nppr(n,kkk)
c              write(ifch,*)'remn',irem,' (',jrem,' )     pom',npom
c     &            ,'    ',zzremn(irem,jrem)
            ie=0
            do is=ifrptl(1,npom),ifrptl(2,npom)
              if(ie.eq.0)is1=is
              if(idptl(is).ne.9)ie=ie+1
              if(ie.eq.2)then
               is2=is
               ie=0
               if(ir.eq. 1)zpaptl(is1)=zzremn(irem,jrem)
               if(ir.eq.-1)zpaptl(is2)=zzremn(irem,jrem)
               do isi=is1,is2
c                     write(ifch,*)'            ',isi,idptl(isi),zpaptl(isi)
               enddo
              endif
            enddo
          endif
        enddo
      enddo
      
      end

c-----------------------------------------------------------------------
      subroutine ProReM(ir,irem,iret)
c-----------------------------------------------------------------------
c propose remnant mass of remnant irem in case of proj (ir=1) 
c or target (ir=-1)
c   (-> xmp, xpt)
c iret : input : if iret=10 force to give mass even if no more energy, 
c        when input not 10 : output = error if 1
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      double precision rr,xxx,xmin,xmax,msmin,xmmin,xpt2rem,xtest0
      double precision at,alp,xi,xii,eps,sx,xmin0,xtest(mamx),fxtest
      parameter(eps=1.d-20)
      common/cemsr5/at(0:1,0:6)
      double precision plc,s,p5sq,aremn,aremnex
      common/cems5/plc,s
      integer icrmn(2)
      logical cont,force
      character cremn*4
      
      call utpri('ProReM',ish,ishini,5)

      if(iret.eq.10)then
        force=.true.
      else
        iret=0
        force=.false.
      endif
      ntrymx=50       

c uncomment the following two lines to force the excitation 
     
ccc      force=.true.   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ccc      ntrymx=1       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c initial definitions

      ntry=0
      xxx=0.d0
      if(ir.eq.1)then
        cremn='proj'
        jrem=1
        jremo=2
        masso=matarg
        do j=1,masso
          xme(j)=0.d0
        enddo
        amremn=amproj
           !idx=isign(iabs(idproj)/10*10+1,idproj)
           !call idmass(idx,amremn)
        iremo1=itarg(1)
ctp        kolzi=kolp(irem)
        noevt=ntgevt
        msmin=dble(amremn*amremn)
        if(iep(irem).eq.6)goto 678
      elseif(ir.eq.-1)then
        cremn='targ'
        jrem=2
        jremo=1
        masso=maproj
        do j=1,masso
          xme(j)=0.d0
        enddo
        amremn=amtarg
           !idx=isign(iabs(idtarg)/10*10+1,idtarg)
           !call idmass(idx,amremn)
        iremo1=iproj(1)
ctp        kolzi=kolt(irem)
        noevt=npjevt
        msmin=dble(amremn*amremn)
        if(iet(irem).eq.6)goto 678
      endif 
      
ctp   noevt replace noxevt
ctp if iez=0, 5% energy violation allowed to give mass to the other side
ctp      noxevt=0      !?????? otherwise, energy is strongly not conserved
ctp      do i=1,masso
ctp       if(iez(i,jremo).gt.0)noxevt=noxevt+1
ctp      enddo
       

c ntry
      
    1 ntry=ntry+1   
      if(ntry.gt.ntrymx)then   
        if(ish.ge.5)then
          call utmsg('ProReM')
          write(ifch,*)'Remnant mass assignment not possible (ntry)'
          if(force)write(ifch,*)'Ignore p4 conservation'
          call utmsgf
        endif
        if(.not.force)then
          iret=1
        else
          if(ir.eq.1)then
            xmp(irem)=xxx
          else
            xpt(irem)=xxx
          endif
        endif
        goto 1000
      endif

c check

      if(xpz(irem,jrem).le.0.d0)then
        write(ifch,*)'ProRem ipp',xpz(irem,jrem),irem,lremn(irem,jrem)
        do li=1,lremn(irem,jrem)
          kkk=kremn(irem,li,jrem)
          write(ifch,*)'kkk',kkk
        enddo
        call XPrint('ProRem pro:&')
        call utstop('Big problem in ProRem pro!&')
      endif

c xtest = xminus-max,  corresponding mostly to a remnant mass 0.2

      xtest0=0.d0
      fxtest=0.2d0
      do j=1,masso
        cont=.false.
        xme(j)=0.d0
ctp        if(xmz(j,jremo).gt.eps.and.iez(j,jrem).gt.0)then !xmz(,jremo)=xplus
ctp060824        if(xmz(j,jremo).gt.eps.and.iez(j,jrem).ge.0)then !xmz(,jremo)=xplus
c        if(iez(j,jremo).gt.0.or.koll.eq.1)then !xmz(,jremo)=xplus
          if(xmz(j,jremo).gt.eps)then !xmz(,jremo)=xplus
            cont=.true.
            xmmin=xzos(j,jremo)/xmz(j,jremo)
          else
            xmmin=xzos(j,jremo)
          endif
          xtest(j)=xpz(j,jremo)-xmmin !maximal momentum available
!this term is very important for non excited remnants in pp, it changes the xf 
! distribution of proton and the multiplicity at low energy. Fxtest should not
! be to close to 0. otherwise it makes a step in xf distribution of p at
! 1-fxtest but if fxtest=1, multiplicity at low energy is too high ...
          if(.not.cont)then
            if(xtest(j).gt.0d0)then
              xtest(j)=min(xtest(j),fxtest/xpz(irem,jrem))
            else
              xtest(j)=min(1.d0,fxtest/xpz(irem,jrem))
            endif
          endif
c        else
c          xtest(j)=0.01d0 !maximal momentum available for non exited state
c        endif
         xtest0=max(xtest0,xtest(j))
c        print *,iep(1),iet(1),iez(irem,jrem),xtest(j),xpz(j,jremo),xmmin
c     & ,xzos(j,jremo),xmz(j,jremo)
      enddo
ctp060824      if(.not.cont)xtest=min(1.d0,0.2d0/xpz(irem,jrem))  

       
      cont=.true.

c defs

      sx=s*xpz(irem,jrem)
      icrmn(1)=icremn(1,irem,jrem)
      icrmn(2)=icremn(2,irem,jrem)
      xpt2rem=xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0

c  fremnux (+) or fremnux2 (-) ?

      aremn=dble(fremnux2(icrmn)) !dble(max(amremn,fremnux2(icrmn)))
c      if(iez(j,jrem).eq.2)then
c        aremnex=aremn   
c      else
c        aremnex=max(amzmn(idz(irem,jrem),jrem)   !makes remnant to heavy at low energy
c     &   +amemn(idz(irem,jrem),iez(irem,jrem)) 
c     &           ,dble(fremnux(icrmn)))     
        aremnex=aremn+amemn(idz(irem,jrem),iez(irem,jrem))
c      endif

c determine xminus 

c      xmin0=1.05*(aremn**2d0+xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0)/sx
c      xmin=1.1*(aremnex**2d0+xxz(irem,jrem)**2d0+xyz(irem,jrem)**2d0)/sx
      xmin0=1.01d0*(aremn**2d0+xpt2rem)/sx
      xmin=1.01d0*(aremnex**2d0+xpt2rem)/sx
      xmax=min(1.d6/s,xtest0)             !to avoid ultra high mass remnants
c for diffractive remnant, mass should never exceed 5% of the proj or targ energy
      if(iez(irem,jrem).eq.1)then
        xmax=min(xmax,max(dble(xmaxremn),xmin))
      elseif(iez(irem,jrem).eq.2)then
        xmax=min(xmax,max(dble(xmaxdiff),xmin))
      endif
      if(koll.eq.1)xmax=min(xmax,xpz(iremo1,jremo))
      if(xmin.ge.xmax-eps)then
        xmin=xmin0
        if(koll.ne.1)xmax=1.d0
        if(xmin.ge.xmax-eps)then
        if(.not.force)then
          iret=1
        else
          xmz(irem,jrem)=xmin
        endif
        goto 1000
        endif
      endif
      xmin0=xmin
      rr=dble(rangen())
      if(iez(irem,jrem).gt.0)then
c        xmin=xmin-xpt2rem/sx                     !no pt
c        xmax=xmax-xpt2rem/sx                     !no pt
        alp=at(idz(irem,jrem),iez(irem,jrem))
c        print *,'alp',iez(irem,jrem),alp,xmin,xmax
        if(dabs(alp-1.d0).lt.eps)then
          xxx=xmax**rr*xmin**(1d0-rr)
        else
          xxx=(rr*xmax**(1d0-alp)+(1d0-rr)*xmin**(1d0-alp))
     &                                             **(1d0/(1d0-alp))
        endif
c        xxx=xxx+xpt2rem/sx                       !no pt
      else
c        xmin=dble(amremn)**2d0/sx                !no pt
c        xxx=xmin+xpt2rem/sx                      !no pt
        xmin=(dble(amremn)**2d0+xpt2rem)/sx
        xxx=xmin
        if(xmin.gt.xmax+eps)then
          if(ish.ge.6)write(ifch,*)'xmin>xmax for proj not possible (2)'
          if(.not.force)then
            iret=1
          else
            xmz(irem,jrem)=xxx
          endif
          goto 1000
        endif
        !to have nice diffractive pic, do not allow too much fluctuation
c        xmin0=0.92d0*xxx  
c        xmin0=0.9d0*xxx 
        xmin0=min(0.99d0,1d0-fxtest*dble(1.-rangen()))*xxx
c        xmin0=dble(0.9+0.09*rangen())*xxx
      endif
      xzos(irem,jrem)=xmin0*xpz(irem,jrem)
      msmin=xmin*sx
c      msmin=xmin*sx+xpt2rem                      !no pt

c partition xminus between nucleons of the other side

      xii=1d0
      iimax=noevt                 !number of opposite side participants
ctp      iimax=noxevt                 !number of opposite side participants
      ii=iimax
      iro=int(rangen()*masso)+1   ! choose ramdomly a nucleon to start

      do while(ii.gt.0)

        cont=iez(iro,jremo).lt.0.or.xme(iro).lt.-0.99
        do while(cont) 
          iro=iro+1
          if(iro.gt.masso)iro=iro-masso
          ii=ii-1
          if(ii.lt.1)goto 1
          cont=iez(iro,jremo).lt.0.or.xme(iro).lt.-0.99
        enddo

        if(ii-1.gt.0)then
         xi=xii*dble(rangen())**(1.d0/dble(ii-1))
        else
         xi=0d0
        endif
        xme(iro)=xxx*(xii-xi)

        xmmin=xzos(iro,jremo)
        if(xmz(iro,jremo).gt.eps)then
          xmmin=xmmin/xmz(iro,jremo)
        elseif(koll.eq.1.and.xtest(iro).gt.eps)then
          xmmin=xmmin/min(xpz(irem,jrem),xtest(iro))
        elseif(xtest(iro).gt.eps)then
          xmmin=xmmin/xtest(iro)
        endif
        if((xpz(iro,jremo)-xme(iro)).lt.xmmin)then
c          write(ifch,*)'     skip ',cremn,' ',ii,iimax,ntry,xxx
c     &    ,xpz(iro,jremo)-xme(iro),xmmin
          if(ii.le.1)goto1
          xme(iro)=-1.d0
        else 
          xii=xi
c          write(ifch,*)'       ok ',cremn,' ',ii,iimax,ntry,xme(iro)/xxx
        endif
        iro=iro+1
        if(iro.gt.masso)iro=iro-masso
        ii=ii-1
        
      enddo

c check xmz(irem,jrem)

      xmz(irem,jrem)=xxx

 678  p5sq=xpz(irem,jrem)*plc*xmz(irem,jrem)*plc
c      write(ifch,*)'final mass',p5sq,msmin,xpz(irem,jrem),xmz(irem,jrem)
c     &,force
      if(p5sq.lt.msmin)then
        if(ish.ge.5)then
          call utmsg('ProReM')
          write(ifch,*)'Remnant mass assignment not possible (M<Mmin)!'
          if(force)write(ifch,*)'Ignore p4 conservation'
          call utmsgf
        endif
        if(.not.force)then
          iret=1
          goto 1000
        elseif(xpz(irem,jrem).gt.0.d0)then
          xmz(irem,jrem)=msmin/(plc*plc*xpz(irem,jrem))
        endif
      endif

c subtract xme

      do iro=1,masso
        if(xme(iro).gt.0.d0)then
          xpz(iro,jremo)=xpz(iro,jremo)-xme(iro)  !xpz(,jremo)=xminus
        endif
      enddo
        
 1000 continue
 
      call utprix('ProReM',ish,ishini,5)

      end   
        
c-----------------------------------------------------------------------
      subroutine ProSeTy(k,n)
c-----------------------------------------------------------------------
c creates proposal for string ends, idp., idm.
c updates quark counters
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'

      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      double precision pes,xfqp,xfqt   !so01
      parameter(eps=1.e-6)
      common/ems9/xfqp(0:9),xfqt(0:9)
      common/emsx3/pes(0:3,0:6)
      
      if(idpr(n,k).eq.2)stop'no Reggeons any more'

      ip=iproj(k)
      it=itarg(k)
      
      if(idpr(n,k).eq.3)then
       pssp=0.
       pvsp=0.
       psap=0.
       pddp=0.
       psvvp=0.
       paasp=0.
       psst=0.
       pvst=0.
       psat=0.
       pddt=0.
       psvvt=0.
       paast=0.
       if(idhpr(n,k).eq.3)then  !so01
        idp1pr(n,k)=2
        idp2pr(n,k)=8
        idm1pr(n,k)=2
        idm2pr(n,k)=8
        ivp(ip)=ivp(ip)          !-1
        ivt(it)=ivt(it)          !-1
       elseif(idhpr(n,k).eq.2)then
        idp1pr(n,k)=1
        idp2pr(n,k)=1
        idm1pr(n,k)=2
        idm2pr(n,k)=8
        ivt(it)=ivt(it)          !-1
       elseif(idhpr(n,k).eq.1)then
        idp1pr(n,k)=2
        idp2pr(n,k)=8
        idm1pr(n,k)=1
        idm2pr(n,k)=1
        ivp(ip)=ivp(ip)          !-1
       elseif(idhpr(n,k).eq.0)then
        idp1pr(n,k)=1
        idp2pr(n,k)=1
        idm1pr(n,k)=1
        idm2pr(n,k)=1
       else
        call utstop('ProSeTy-idhpr????&')
       endif


      elseif(idpr(n,k).eq.1)then
       
c    projectile

       if(iabs(idfpr(n,k)).eq.1)then  

       ntry=0
  1    ntry=ntry+1
       if(ntry.gt.10)call utstop('something goes wrong in sr ProSeTy&')
       pss=wgtsea
       pvs=wgtval      !  *ivp(ip)
       psa=wgtval      !  *iap(ip)
       pdd=wgtdiq
       psvv=wgtqqq(iclpro)     !  *ivp(ip)*(ivp(ip)-1)/2.
       paas=wgtqqq(iclpro)    !  *iap(ip)*(iap(ip)-1)/2.
       su=pss+pvs+psa+pdd+psvv+paas
       pssp = pss /su
       pvsp = pvs /su
       psap = psa /su
       pddp = pdd /su
       psvvp= psvv/su
       paasp= paas/su
       r=rangen()
       if(r.gt.(pssp+pvsp+psap+pddp+psvvp).and.paasp.gt.eps)then
        idp1pr(n,k)=5
        idp2pr(n,k)=1
        idsppr(n,k)=6
c        iap(ip)=iap(ip)-2
       elseif(r.gt.(pssp+pvsp+psap+pddp).and.psvvp.gt.eps)then
        idp1pr(n,k)=1
        idp2pr(n,k)=5
        idsppr(n,k)=5
c        ivp(ip)=ivp(ip)-2
       elseif(r.gt.(pssp+pvsp+psap).and.pddp.gt.eps)then
        idp1pr(n,k)=4
        idp2pr(n,k)=4
        idsppr(n,k)=4
       elseif(r.gt.(pssp+pvsp).and.psap.gt.eps)then
        idp1pr(n,k)=1
        idp2pr(n,k)=2
        idsppr(n,k)=2
c        iap(ip)=iap(ip)-1
       elseif(r.gt.pssp.and.pvsp.gt.eps)then
        idp1pr(n,k)=2
        idp2pr(n,k)=1
        idsppr(n,k)=1
c        ivp(ip)=ivp(ip)-1
       elseif(pssp.gt.eps)then
        idp1pr(n,k)=1
        idp2pr(n,k)=1
        idsppr(n,k)=0
       else 
        goto1
       endif

       else     
        idp1pr(n,k)=1
        idp2pr(n,k)=1
        idsppr(n,k)=0
       endif

c    target 

       if(iabs(idfpr(n,k)).eq.1)then  


       ntry=0
  2    ntry=ntry+1
       if(ntry.gt.10)call utstop('something goes wrong in sr ProSeTy&')
       pss=wgtsea
       pvs=wgtval      !   *ivt(it)
       psa=wgtval      !   *iat(it)
       pdd=wgtdiq
       psvv=wgtqqq(icltar)     !   *ivt(it)*(ivt(it)-1)/2.
       paas=wgtqqq(icltar)    !   *iat(it)*(iat(it)-1)/2.
       su=pss+pvs+psa+pdd+psvv+paas
       psst = pss /su
       pvst = pvs /su
       psat = psa /su
       pddt = pdd /su
       psvvt= psvv/su
       paast= paas/su
       r=rangen()
       if(r.gt.(psst+pvst+psat+pddt+psvvt).and.paast.gt.eps)then
        idm1pr(n,k)=5
        idm2pr(n,k)=1
        idstpr(n,k)=6
c        iat(it)=iat(it)-2
       elseif(r.gt.(psst+pvst+psat+pddt).and.psvvt.gt.eps)then
        idm1pr(n,k)=1
        idm2pr(n,k)=5
        idstpr(n,k)=5
c        ivt(it)=ivt(it)-2
       elseif(r.gt.(psst+pvst+psat).and.pddt.gt.eps)then
        idm1pr(n,k)=4
        idm2pr(n,k)=4
        idstpr(n,k)=4
       elseif(r.gt.(psst+pvst).and.psat.gt.eps)then
        idm1pr(n,k)=1
        idm2pr(n,k)=2
        idstpr(n,k)=2
c        iat(it)=iat(it)-1
       elseif(r.gt.psst.and.pvst.gt.eps)then
        idm1pr(n,k)=2
        idm2pr(n,k)=1
        idstpr(n,k)=1
c        ivt(it)=ivt(it)-1
       elseif(psst.gt.eps)then
        idm1pr(n,k)=1
        idm2pr(n,k)=1
        idstpr(n,k)=0
       else 
        goto2
       endif

       else    
        idm1pr(n,k)=1
        idm2pr(n,k)=1
        idstpr(n,k)=0
       endif
       
      elseif(idpr(n,k).eq.0)then
      
        idp1pr(n,k)=0
        idm2pr(n,k)=0
        idp2pr(n,k)=0
        idm1pr(n,k)=0
       
      endif

        if(ish.ge.6)then
      write(ifch,'(a,2(6(f3.2,1x),2x),$)')'ProSeTy ',
     * pssp,pvsp,psap,pddp,psvvp,paasp, psst,pvst,psat,pddt,psvvt,paast
      write(ifch,'(2x,3i3,2x,2(i2,1x,2i2,1x,2i2,2x))')idpr(n,k),n,k
     * ,idsppr(n,k),idp1pr(n,k),idp2pr(n,k),ivp(ip),iap(ip)
     * ,idstpr(n,k),idm1pr(n,k),idm2pr(n,k),ivt(it),iat(it)
        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine ProSeF(k,n,iret)
c-----------------------------------------------------------------------
c starting from string properties as already determined in EMS,
c one determines string end flavors 
c by checking compatibility with remnant masses.
c strings are written to /cems/ and then to /cptl/
c remnant ic is updated (icproj,ictarg)
c------------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      
      double precision plc,s,pstg,pend
      common/cems5/plc,s
      common/cems/pstg(5,2),pend(4,4),idend(4)      
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icp(2),ict(2),ic(2),icp1(2),icp2(2),icm1(2),icm2(2)
      integer icini(2)
      integer jcp1(nflav,2),jcp2(nflav,2),jcm1(nflav,2),jcm2(nflav,2)
      common/col3/ncol,kolpt /cfacmss/facmss /cts/its

      call utpri('ProSeF',ish,ishini,6)
      
c     entry
c     -----

      iret=0
      
      if(ncol.eq.0)return
      if(itpr(k).ne.1)return
      
      ip=iproj(k)
      it=itarg(k)
       
      if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)return
      if(idpr(n,k).eq.2)stop'Reggeon'
      if(idpr(n,k).eq.3)goto 1000  
      if(ish.ge.5)then
          write(ifch,*)'soft Pomeron'
          write(ifch,*)'k:',k,'  n:',n,'  ip:',ip,'  it:',it
      endif
      np=nppr(n,k)
        
c         string ends

          pend(1,1)=xxp1pr(n,k)
          pend(2,1)=xyp1pr(n,k)
          pend(3,1)=xp1pr(n,k)*plc/2d0
          pend(4,1)=dsqrt(pend(1,1)**2+pend(2,1)**2+pend(3,1)**2)
          pend(1,2)=xxp2pr(n,k)
          pend(2,2)=xyp2pr(n,k)
          pend(3,2)=xp2pr(n,k)*plc/2d0
          pend(4,2)=dsqrt(pend(1,2)**2+pend(2,2)**2+pend(3,2)**2)
          pend(1,4)=xxm1pr(n,k)
          pend(2,4)=xym1pr(n,k)
          pend(3,4)=-xm1pr(n,k)*plc/2d0
          pend(4,4)=dsqrt(pend(1,4)**2+pend(2,4)**2+pend(3,4)**2)
          pend(1,3)=xxm2pr(n,k)
          pend(2,3)=xym2pr(n,k)
          pend(3,3)=-xm2pr(n,k)*plc/2d0
          pend(4,3)=dsqrt(pend(1,3)**2+pend(2,3)**2+pend(3,3)**2)
          
c         strings

          pstg(1,1)=xxp1pr(n,k)+xxm2pr(n,k)
          pstg(2,1)=xyp1pr(n,k)+xym2pr(n,k)
          pstg(3,1)=(xp1pr(n,k)-xm2pr(n,k))*plc/2d0
          pstg(4,1)=(xp1pr(n,k)+xm2pr(n,k))*plc/2d0
          pstg(5,1)=dsqrt((pstg(4,1)-pstg(3,1))*(pstg(4,1)+pstg(3,1))
     &                   -pstg(1,1)**2-pstg(2,1)**2)
          pstg(1,2)=xxp2pr(n,k)+xxm1pr(n,k)
          pstg(2,2)=xyp2pr(n,k)+xym1pr(n,k)
          pstg(3,2)=(xp2pr(n,k)-xm1pr(n,k))*plc/2d0
          pstg(4,2)=(xp2pr(n,k)+xm1pr(n,k))*plc/2d0
          pstg(5,2)=dsqrt((pstg(4,2)-pstg(3,2))*(pstg(4,2)+pstg(3,2))
     &                   -pstg(2,2)**2-pstg(1,2)**2)

c         initialize

          ntry=0
  777     ntry=ntry+1
          if(ntry.gt.100)goto1001
          
          do i=1,2
           icp(i)=icproj(i,ip)
           ict(i)=ictarg(i,it)
           icp1(i)=0
           icp2(i)=0
           icm1(i)=0
           icm2(i)=0
           do j=1,nflav
            jcp1(j,i)=0
            jcp2(j,i)=0
            jcm1(j,i)=0
            jcm2(j,i)=0
           enddo
          enddo
          idpj0=idtr2(icp)
          idtg0=idtr2(ict)
          do j=1,4
           idend(j)=0
          enddo

          if(ish.ge.5)write(ifch,'(a,3x,2i7,i9)')' proj: '
     *     ,(icp(l),l=1,2),idpj0
          if(ish.ge.5)write(ifch,'(a,3x,2i7,i9)')' targ: '
     *    ,(ict(l),l=1,2),idtg0

c         determine string flavors

          call fstrfl(icp,ict,icp1,icp2,icm1,icm2
     *               ,idp1pr(n,k),idp2pr(n,k),idm1pr(n,k),idm2pr(n,k)
     *               ,iabs(idfpr(n,k)),iret)
          if(iret.ne.0)then
            jerr(1)=jerr(1)+1  ! > 9 quarks per flavor attempted. 
                               !    OK when happens rarely.
            goto 1001
          endif

c         check mass string 1

          ic(1)=icp1(1)+icm2(1)
          ic(2)=icp1(2)+icm2(2)
          if(ic(1).gt.0.or.ic(2).gt.0)then
           am=sngl(pstg(5,1))
           call iddeco(icp1,jcp1)
           call iddeco(icm2,jcm2)
           ammns=utamnx(jcp1,jcm2)
           if(ish.ge.7)write(ifch,'(a,2i7,2e12.3)')
     *           ' string 1 - ic,mass,min.mass:',ic,am,ammns
           if(am.lt.ammns*facmss)then
             goto 777   !avoid virpom
           endif
           idend(1)=idtra(icp1,0,0,3)
           idend(3)=idtra(icm2,0,0,3)
           if(ish.ge.7)write(ifch,'(a,2i4)') ' string 1 - SE-ids:'
     *      ,idend(1),idend(3)
          endif

c         check mass string 2

          ic(1)=icp2(1)+icm1(1)
          ic(2)=icp2(2)+icm1(2)
          if(ic(1).gt.0.or.ic(2).gt.0)then
           am=sngl(pstg(5,2))
           call iddeco(icp2,jcp2)
           call iddeco(icm1,jcm1)
           ammns=utamnx(jcp2,jcm1)
           if(ish.ge.7)write(ifch,'(a,2i7,2e12.3)')
     *           ' string 2 - ic,mass,min.mass:',ic,am,ammns
           if(am.lt.ammns*facmss)then
             goto 777  !avoid virpom
           endif
           idend(2)=idtra(icp2,0,0,3)
           idend(4)=idtra(icm1,0,0,3)
           if(ish.ge.7)write(ifch,'(a,2i4)') ' string 2 - SE-ids:'
     *      ,idend(2),idend(4)
          endif

          if(ish.ge.5)then
          write(ifch,'(a,i10)')' pom:   '
     *    ,idptl(np)
          write(ifch,'(a,2i5)')' str 1: ' 
     *    ,idend(1),idend(3)
          write(ifch,'(a,2i5)')' str 2: '
     *    ,idend(2),idend(4)
          write(ifch,'(a,2i7,1x,a)')' proj:  '
     *    ,(icp(l),l=1,2)
          write(ifch,'(a,2i7,1x,a)')' targ:  '
     *    ,(ict(l),l=1,2)
          endif

c         update remnant ic

          do i=1,2
           icproj(i,ip)=icp(i)
           ictarg(i,it)=ict(i)
          enddo

          call idtr4(idptl(ip),icini)   !excited remnant ?
          if(ish.ge.5)write(ifch,*)'icini proj',icini
     &    ,(icp(1)-icini(1)),(icp(2)-icini(2))
          if((icp(1)-icini(1))+(icp(2)-icini(2)).ne.0)iep(ip)=1 
          call idtr4(idptl(maproj+it),icini)
          if(ish.ge.5)write(ifch,*)'icini targ',icini
     &    ,(ict(1)-icini(1)),(ict(2)-icini(2))
          if((ict(1)-icini(1))+(ict(2)-icini(2)).ne.0)iet(it)=1 
          if(ish.ge.5)write(ifch,*)'iep,iet ',iep(ip),iet(it)

c         write strings to /cptl/

          its=idp1pr(n,k)+idm2pr(n,k)
          call fstrwr(1,1,3,k,n)
          its=idp2pr(n,k)+idm1pr(n,k)
          call fstrwr(2,2,4,k,n)
          
c     exit
c     ----

1000  continue
      call utprix('ProSeF',ish,ishini,6)
      return

1001  iret=1
      goto1000

      end

c-----------------------------------------------------------------------
      subroutine fstrfl(icp,ict,icp1,icp2,icm1,icm2
     *                         ,idp1,idp2,idm1,idm2,idfp,iret)
c-----------------------------------------------------------------------
c knowing the string end types (idp1,idp2,idm1,idm2) 
c               and remnant flavors (icp,ict)
c               and remnant link of the string (idfp)
c one determines quark flavors of string ends (icp1,icp2,icm1,icm2)
c               and updates remnant flavors (icp,ict)   
c iret=0   ok
c iret=1   problem, more than 9 quarks per flavor attempted
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer icp(2),ict(2),icp1(2),icp2(2),icm1(2),icm2(2)
      integer jcp(6,2),jct(6,2),jcpi(6,2),jcti(6,2)
      integer iq(2,4)
c      data neuz/0/proz/0/dtaz/0/
c      save neuz,proz,dtaz

      call utpri('fstrfl',ish,ishini,7)
      
c     entry
c     -----

      iret=0
      iret1=0
      iret2=0
      iret3=0
      iret4=0
      
      if(idfp.eq.2)stop'fstrfl: should not happen (2).      '
      if(idfp.eq.3)stop'fstrfl: should not happen (3).      '
      if(idp1.eq.4)stop'fstrfl: diq code 4 not used any more'
      if(idm1.eq.4)stop'fstrfl: diq code 4 not used any more'
      if(idp2.eq.4)stop'fstrfl: diq code 4 not used any more'
      if(idm2.eq.4)stop'fstrfl: diq code 4 not used any more'
      if(idp1.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idp2.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idm1.eq.8)stop'fstrfl: fragm quarks not used any more'
      if(idm2.eq.8)stop'fstrfl: fragm quarks not used any more'

c determine flavors of string ends (u,d,s)

      call iddeco(icp,jcpi)
      call iddeco(ict,jcti)
      call iddeco(icp,jcp)
      call iddeco(ict,jct)
      if(ish.ge.7)then
       write(ifch,'(a,2i7,5x,6i2,3x,6i2,3x,i1)')' proj:',icp,jcp
       write(ifch,'(a,2i7,5x,6i2,3x,6i2,3x,i1)')' targ:',ict,jct
      endif

c empty

      if(idp1.eq.0)then
       iq(1,1)=0
       iq(2,1)=0
      endif
      if(idp2.eq.0)then
       iq(1,2)=0
       iq(2,2)=0
      endif
      if(idm1.eq.0)then
       iq(1,4)=0
       iq(2,4)=0
      endif
      if(idm2.eq.0)then
       iq(1,3)=0
       iq(2,3)=0
      endif

c valence quarks

      if(idp1.eq.2)then
       iq(1,1)=idrafl(iclpro,jcp,1,'s',iret)
       iq(2,1)=0
      endif
      if(idp2.eq.2)then
       iq(1,2)=idrafl(iclpro,jcp,2,'s',iret)
       iq(2,2)=0
      endif
      if(idm1.eq.2)then
       iq(1,4)=idrafl(icltar,jct,1,'s',iret)
       iq(2,4)=0
      endif
      if(idm2.eq.2)then
       iq(1,3)=idrafl(icltar,jct,2,'s',iret)
       iq(2,3)=0
      endif

c sea quarks      

      if(idp1.eq.1)then
       iq(1,1)=idrafl(iclpro,jcp,1,'s',iret1)  
       iq(2,1)=0
      endif
      if(idm1.eq.1)then
       iq(1,4)=idrafl(icltar,jct,1,'s',iret4)  
       iq(2,4)=0
      endif
      if(idp2.eq.1)then 
       iq(1,2)=idrafl(iclpro,jcp,2,'s',iret2)
       iq(2,2)=0
      endif
      if(idm2.eq.1)then
       iq(1,3)=idrafl(icltar,jct,2,'s',iret3)
       iq(2,3)=0
      endif

c diquarks, code 5 (former valence, but actually sea)

      if(idp1.eq.5)then
c       fc=puds
c       iq(1,1)=idraflx(fc,iclpro,jcp,2,'s',iret)
c       if(iq(1,1).eq.3)fc=fc*puds
c       iq(2,1)=idraflx(fc,iclpro,jcp,2,'s',iret)    
       iq(1,1)=idrafl(iclpro,jcp,2,'d',iret)
       iq(2,1)=idrafl(iclpro,jcp,2,'d',iret)    
      endif
      if(idm1.eq.5)then
c       fc=puds
c       iq(1,4)=idraflx(fc,icltar,jct,2,'s',iret)
c       if(iq(1,4).eq.3)fc=fc*puds
c       iq(2,4)=idraflx(fc,icltar,jct,2,'s',iret)
       iq(1,4)=idrafl(icltar,jct,2,'d',iret)
       iq(2,4)=idrafl(icltar,jct,2,'d',iret)
      endif
      if(idp2.eq.5)then
c       fc=puds
c       iq(1,2)=idraflx(fc,iclpro,jcp,1,'s',iret)
c       if(iq(1,2).eq.3)fc=fc*puds
c       iq(2,2)=idraflx(fc,iclpro,jcp,1,'s',iret)
       iq(1,2)=idrafl(iclpro,jcp,1,'d',iret)
       iq(2,2)=idrafl(iclpro,jcp,1,'d',iret)
      endif
      if(idm2.eq.5)then
c       fc=puds
c       iq(1,3)=idraflx(fc,icltar,jct,1,'s',iret)
c       if(iq(1,3).eq.3)fc=fc*puds
c       iq(2,3)=idraflx(fc,icltar,jct,1,'s',iret)
       iq(1,3)=idrafl(icltar,jct,1,'d',iret)
       iq(2,3)=idrafl(icltar,jct,1,'d',iret)
      endif
       
      if(iret.ne.0)goto 1000


c in case of saturated remnants, use the same flavor for quark and anti-quark
c at string-end
      if(iret1.ne.0.and.iret2.ne.0)then
        call iddeco(icp,jcp)
        if(rangen().gt.0.5)then
          iq(1,2)=iq(1,1)
        else
          iq(1,1)=iq(1,2)
        endif
      elseif(iret1.eq.0.and.iret2.ne.0.and.idp1.eq.1)then
        call iddeco(icp,jcp)
        iq(1,2)=iq(1,1)
        if(idp1.eq.4)iq(2,1)=iq(1,1)
c        if(idp2.eq.4)iq(2,2)=iq(1,2)
      elseif(iret2.eq.0.and.iret1.ne.0.and.idp2.eq.1)then
        call iddeco(icp,jcp)
        iq(1,1)=iq(1,2)
        if(idp1.eq.4)iq(2,1)=iq(1,1)
c        if(idp2.eq.4)iq(2,2)=iq(1,2)
      elseif(iret1.ne.0.or.iret2.ne.0)then
        iret=1
        goto 1000
      endif

      if(iret3.ne.0.and.iret4.ne.0)then
        call iddeco(ict,jct)
        if(rangen().gt.0.5)then
          iq(1,4)=iq(1,3)
        else
          iq(1,3)=iq(1,4)
        endif
      elseif(iret3.eq.0.and.iret4.ne.0.and.idm1.eq.1)then
        call iddeco(ict,jct)
        iq(1,4)=iq(1,3)
c        if(idm2.eq.4)iq(2,3)=iq(1,3)
        if(idm1.eq.4)iq(2,4)=iq(1,4)
      elseif(iret4.eq.0.and.iret3.ne.0.and.idm2.eq.1)then
        call iddeco(ict,jct)
        iq(1,3)=iq(1,4)
c        if(idm2.eq.4)iq(2,3)=iq(1,3)
        if(idm1.eq.4)iq(2,4)=iq(1,4)
      elseif(iret3.ne.0.or.iret4.ne.0)then
        iret=1
        goto 1000
      endif


c determine icp,ict      
 
      call idenco(jcp,icp,iret)
      if(iret.ne.0)goto 1000
      call idenco(jct,ict,iret)
      if(iret.ne.0)goto 1000

      ifla=iq(1,1)
      iflb=iq(2,1)
      iflc=iq(1,3)
      ifld=iq(2,3)
      if(ish.ge.7)write(ifch,'(a,2i3,4x,2i3)')
     *' string 1, string ends:',ifla,iflb,iflc,ifld

      if(ifla.gt.0)then
       if(iflb.eq.0)then
        icp1(1)=10**(6-ifla)
        icp1(2)=0
       else
        icp1(1)=0
        icp1(2)=10**(6-ifla)
        icp1(2)=icp1(2)+10**(6-iflb)
       endif
      endif

      if(iflc.gt.0)then
       if(ifld.eq.0)then
        icm2(1)=0
        icm2(2)=10**(6-iflc)
       else
        icm2(1)=10**(6-iflc)
        icm2(1)=icm2(1)+10**(6-ifld)
        icm2(2)=0
       endif
      endif

      ifla=iq(1,4)
      iflb=iq(2,4)
      iflc=iq(1,2)
      ifld=iq(2,2)
      if(ish.ge.7)write(ifch,'(a,2i3,4x,2i3)')
     *' string 2, string ends:',ifla,iflb,iflc,ifld

      if(ifla.gt.0)then
       if(iflb.eq.0)then
        icm1(1)=10**(6-ifla)
        icm1(2)=0
       else
        icm1(1)=0
        icm1(2)=10**(6-ifla)
        icm1(2)=icm1(2)+10**(6-iflb)
       endif
      endif

      if(iflc.gt.0)then
       if(ifld.eq.0)then
        icp2(1)=0
        icp2(2)=10**(6-iflc)
       else
        icp2(1)=10**(6-iflc)
        icp2(1)=icp2(1)+10**(6-ifld)
        icp2(2)=0
       endif
      endif

      if(ish.ge.7)then
        write(ifch,'(a,2i7,4x,2i7)')
     *  ' SE-forw:',icp1(1),icp1(2),icp2(1),icp2(2)
        write(ifch,'(a,2i7,4x,2i7)')
     *  ' SE-back:',icm1(1),icm1(2),icm2(1),icm2(2)
        write(ifch,'(a,2i7,5x,6i2,3x,6i2)')' proj:',icp,jcp
        write(ifch,'(a,2i7,5x,6i2,3x,6i2)')' targ:',ict,jct
      endif
      
c     exit
c     ----

1000  continue
      call utprix('fstrfl',ish,ishini,7)
      return
      end

c-----------------------------------------------------------------------
      integer function jdrafl(icl,jc,mod,iret)
c-----------------------------------------------------------------------
c mod=1
c returns random flavor of a quark
c
c mod=2
c jc : quark content of remnant
c returns random flavor and  update remant with corresponding q-qbar pair \
c if there is  enough place (else iret=1)
c
c     id=1 u, id=2 d, id=3 s 
c-----------------------------------------------------------------------
      include 'epos.inc'
      integer jc(nflav,2)

c        write(*,*)'entry jdrafl, j,c,jc: ',j,c,jc

       pu=rstrau(icl)
       pd=rstrad(icl)
       ps=rstras(icl)

       s=pu+pd+ps
       if(s.gt.0.)then
         r=rangen()*s
         if(r.gt.(pu+pd).and.ps.gt.1d-10)then
           i=3
         elseif(r.gt.pu.and.pd.gt.1d-10)then
           i=2
         else
           i=1
         endif
       else
         i=1+int((2.+rstras(icl))*rangen())
       endif
       jdrafl=i

c      write(*,*)'jc before updating',jc
c      write(*,*)'i,j,jc',i,j,jc

       if(mod.eq.2)then
         call idsufl2(i,1,jc,iret)
         call idsufl2(i,2,jc,iret)
       endif

      return
      end


cc-----------------------------------------------------------------------
c      subroutine fremfl(icp,ict,iret)
cc-----------------------------------------------------------------------
cc checks projectile and target flavor (icp,ict)
cc in case of reggeon exchange they do not correspond to hadrons.
cc one transfers therefore flavor from one side to the other in order
cc to have hadron flavor.
cc icp and ict are modified correspondingly
cc-----------------------------------------------------------------------
c      include 'epos.inc'
c      integer icp(2),ict(2),jcp(6,2),jct(6,2),kp(4),kt(4)
c
c      call utpri('fremfl',ish,ishini,7)
c      
cc     entry
cc     -----
c
c      iret=0
c
c      call iddeco(icp,jcp)
c      call iddeco(ict,jct)
c
c      iakp=0
c      iakt=0
c      ikp=0
c      ikt=0
c      do l=1,4
c       kp(l)=jcp(l,1)-jcp(l,2)
c       kt(l)=jct(l,1)-jct(l,2)
c       iakp=iakp+iabs(kp(l))
c       iakt=iakt+iabs(kt(l))
c       ikp=ikp+kp(l)
c       ikt=ikt+kt(l)
c      enddo
c      if(ish.ge.7)write(ifch,*)'iak_p:',iakp,' ik_p:',ikp
c      if(ish.ge.7)write(ifch,*)'iak_t:',iakt,' ik_t:',ikt
c
c      if(iakp.eq.4)then
c       if(ikp.eq.4.or.ikp.eq.-2)then
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       elseif(ikp.eq.-4.or.ikp.eq.2)then
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       else
c        call utstop('fremfl&')
c       endif
c      elseif(iakt.eq.4)then
c       if(ikt.eq.4.or.ikt.eq.-2)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       elseif(ikt.eq.-4.or.ikt.eq.2)then
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       else
c        call utstop('fremfl&')
c       endif
c      elseif(iakp.eq.3)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       else
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       endif
c      elseif(iakt.eq.3)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       endif
c      elseif(iakp.eq.2)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       endif
c      elseif(iakt.eq.2)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       else
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       endif
c      elseif(iakp.eq.1)then
c       if(ikp.gt.0)then
c        ifl=idrafl(jcp,2,'v',iret)
c        iqp=2      ! add antiquark
c        iqt=1      ! subtract antiquark
c       else
c        ifl=idrafl(jcp,1,'v',iret)
c        iqp=1      ! add quark
c        iqt=2      ! subtract quark
c       endif
c      elseif(iakt.eq.1)then
c       if(ikt.gt.0)then
c        ifl=idrafl(jct,2,'v',iret)
c        iqp=1      ! subtract antiquark
c        iqt=2      ! add antiquark
c       else
c        ifl=idrafl(jct,1,'v',iret)
c        iqp=2      ! subtract quark
c        iqt=1      ! add quark
c       endif
c      else
c       call utstop('fremfl: error&')
c      endif
c
c      if(ish.ge.7)write(ifch,*)'iq_p:',iqp,' iq_t:',iqt,' if:',ifl
c      call uticpl(icp,ifl,iqp,iret) 
c      if(iret.ne.0)goto1000
c      call uticpl(ict,ifl,iqt,iret)
c      if(iret.ne.0)goto1000
c
cc     exit
cc     ----
c
c1000  continue
c      call utprix('fremfl',ish,ishini,7)
c      return
c      end
c
c-----------------------------------------------------------------------
      subroutine fstrwr(j,ii,jj,k,n)
c-----------------------------------------------------------------------
c take pstg(5,j),pend(4,ii),idend(ii),pend(4,jj),idend(jj)  (/cems/)
c and write it to /cptl/
c-----------------------------------------------------------------------
c  j:     string 1 or 2
c  ii,jj: string end (1,2: proj; 3,4: targ) 
c  k:     current collision
c  n:     current pomeron
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      
      double precision pstg,pend
      common/cems/pstg(5,2),pend(4,4),idend(4)      
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision  pp(4)
      common/cts/its

      call utpri('fstrwr',ish,ishini,7) 

      if(idend(ii).ne.0.and.idend(jj).ne.0)then

c string

       call utlob2(1,pstg(1,j),pstg(2,j),pstg(3,j),pstg(4,j),pstg(5,j)
     * ,pend(1,ii),pend(2,ii),pend(3,ii),pend(4,ii),20)
       pp(1)=0d0
       pp(2)=0d0
       pp(3)=.5d0*pstg(5,j)
       pp(4)=.5d0*pstg(5,j)
       call utrot2
     * (-1,pend(1,ii),pend(2,ii),pend(3,ii),pp(1),pp(2),pp(3))
       call utlob2(-1,pstg(1,j),pstg(2,j),pstg(3,j),pstg(4,j),pstg(5,j)
     * ,pp(1),pp(2),pp(3),pp(4),21)

       npom=nppr(n,k)
       if(ifrptl(1,npom).eq.0)ifrptl(1,npom)=nptl+1
       ifrptl(2,npom)=nptl+2
       istptl(npom)=31
       
       nptl=nptl+1
       pptl(1,nptl)=sngl(pp(1))
       pptl(2,nptl)=sngl(pp(2))
       pptl(3,nptl)=sngl(pp(3))
       pptl(4,nptl)=sngl(pp(4))
       pptl(5,nptl)=0.
       istptl(nptl)=20
       iorptl(nptl)=npom
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coord(1,k)
       xorptl(2,nptl)=coord(2,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=xorptl(4,nptl)
       tivptl(2,nptl)=xorptl(4,nptl)
       idptl(nptl)=idend(ii)
       ityptl(nptl)=ityptl(npom)+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(nptl)=0.
       
       nptl=nptl+1
       do i=1,4
        pptl(i,nptl)=sngl(pstg(i,j))-pptl(i,nptl-1)
       enddo
       pptl(5,nptl)=0.

       istptl(nptl)=20
       iorptl(nptl)=nppr(n,k)
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coord(1,k)
       xorptl(2,nptl)=coord(2,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=xorptl(4,nptl)
       tivptl(2,nptl)=xorptl(4,nptl)
       idptl(nptl)=idend(jj)
       ityptl(nptl)=ityptl(npom)+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(nptl)=0.

       if(ish.ge.7)then
        write(ifch,100)' kink:',(pptl(l,nptl-1),l=1,4),idptl(nptl-1)
        write(ifch,100)' kink:',(pptl(l,nptl),l=1,4),idptl(nptl)
       endif

      elseif(idend(ii).ne.0.and.idend(jj).eq.0)then

c resonance

       npom=nppr(n,k)
       if(ifrptl(1,npom).eq.0)ifrptl(1,npom)=nptl+1
       ifrptl(2,npom)=nptl+1
       istptl(npom)=31

       nptl=nptl+1
       idptl(nptl)=idend(ii)
       pptl(1,nptl)=sngl(pstg(1,j))
       pptl(2,nptl)=sngl(pstg(2,j))
       pptl(3,nptl)=sngl(pstg(3,j))
       pptl(4,nptl)=sngl(pstg(4,j))
       pptl(5,nptl)=sngl(pstg(5,j))
       istptl(nptl)=0
       iorptl(nptl)=npom
       jorptl(nptl)=0
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coord(1,k)
       xorptl(2,nptl)=coord(2,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=coord(4,k)
       call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
       tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
       ityptl(nptl)=ityptl(npom)+2+j
       itsptl(nptl)=its
       rinptl(nptl)=-9999
       qsqptl(nptl)=0.
       zpaptl(nptl)=0.

       if(ish.ge.7)then
        write(ifch,100)'  res:',(pptl(l,nptl),l=1,4),idptl(nptl)
       endif
      elseif(idend(ii).eq.0.and.idend(jj).eq.0)then
       goto1000
      else
       call utstop('error in fstrwr&')
      endif

  100 format(a,4e9.3,i5)
      
1000  continue 
      call utprix('fstrwr',ish,ishini,7)
      return
      end 

c-----------------------------------------------------------------------
      subroutine ProReF(ir,m)
c-----------------------------------------------------------------------
c  proposes flavor for remnant m for proj (ir=1) or target (ir=-1) 
c  and writes remnant into /cptl/ as string or hadron 
c   ityptl definitions:
c      51  41  ...  rmn drop                        
c      52  42  ...  rmn str inel
c      53  43  ...  rmn str diff
c      54  44  ...  rmn str after droplet or hadron split
c      55  45  ...  rmn res
c      56  46  ...  rmn res after droplet or hadron split 
c      57  47  ...  rmn res after all Pomeron killed 
c      58  48  ...  rmn res from diff
c      59  49  ...  hadron split
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      
      double precision plc,s   !,ptt1,ptt2
      common/cems5/plc,s
      common/cdfptl/idfptl(mxptl)
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat,zor,tor
      common/cttaus/tpro,zpro,ttar,ztar,ttaus,detap,detat 
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)
      double precision amasmin,amasini,mdrmax
      integer icf(2),icb(2)
      integer jcf(nflav,2),jcdummy(nflav,2)
      logical gdrop, ghadr,gproj,gtarg
      double precision ept(5),ep(4),aa(5),am2t
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon    
                                                                       
      call utpri('ProReF',ish,ishini,3)

      if(ir.ne.1.and.ir.ne.-1)stop'ProReF: wrong ir'
      
      irmdropx=irmdrop
 55   idrop=0
      gdrop=.false.
      ghadr=.false.
      iret=0
      dens=0.15
      
      if(ir.eq.1)then
c        if(kolp(m).le.0)goto1000
        if(iep(m).le.-1)goto1000
        gproj=.true.
        gtarg=.false.
        mm=npproj(m)
        iept=iep(m)
        zz=zzremn(m,1)
        iclpt=iclpro
      elseif(ir.eq.-1)then
c        if(kolt(m).le.0)goto1000  
        if(iet(m).le.-1)goto1000  
        gproj=.false.
        gtarg=.true.
        mm=nptarg(m)
        iept=iet(m)
        zz=zzremn(m,2)
        iclpt=icltar
      else
        call utstop('ProReF: ir ???&')
      endif        
      if(ish.ge.3)write(ifch,*)'remnant particle index:',mm
      
      if(ish.ge.8)call alist('ProRef&',1,nptl)
      antotre=antotre+1

      mmini=mm
      nptlini=nptl
      minfra=min(minfra,nptlini)   !for trigger condition
      
      do l=1,5
       ept(l)=dble(pptl(l,mm))
      enddo

      ifrptl(1,mm)=0
      ifrptl(2,mm)=0

c  initialize forward and backward ic (to transform remnant into string)

      if(gproj)then
        icf(1)=icproj(1,m)
        icf(2)=icproj(2,m)
      else                     !gtarg
        icf(1)=ictarg(1,m)
        icf(2)=ictarg(2,m)
      endif
      icb(1)=0
      icb(2)=0

      call iddeco(icf,jcf)
      call idquacjc(jcf,nqu,naq)

c define masses

      amasmin=dble(fremnux2(icf))**2.d0
      if(ept(5).le.0.)then
        ept(5)=sqrt(2*amasmin)
        if(ish.ge.1)then
          call utmsg('ProReF') 
          write(ifch,*)'zero remnant mass -> amasmin'
          call utmsgf
        endif
      endif
      am2t=(ept(4)+ept(3))*(ept(4)-ept(3))-(ept(1)**2+ept(2)**2)
      if(ish.ge.1
     &   .and.(am2t.le.0d0.or.abs(am2t-ept(5)*ept(5)).gt.ept(5)))then
          write(ifch,*)'Precision problem in ProRef, p:',
     &             (ept(k),k=1,4),ept(5)*ept(5),am2t
      endif
      ept(4)=sqrt(ept(3)*ept(3)+ept(2)*ept(2)+ept(1)*ept(1)
     &           +ept(5)*ept(5))

      if(ish.ge.3)then
        if(gproj)
     *      write(ifch,'(a,5e11.3,2i7)')' proj:'
     *      ,(sngl(ept(k)) ,k=1,5),(icproj(k,m) ,k=1,2)
        if(gtarg)
     *      write(ifch,'(a,5e11.3,2i7)')' targ:'
     *      ,(sngl(ept(k)) ,k=1,5),(ictarg(k,m),k=1,2)
      endif

      amasini=ept(5)*ept(5)

c      mdrmax=amasmin+dble(amdrmax*amdrmax)
      mdrmax=dble(fremnux(icf)+amdrmax)**2.d0

      if(ish.ge.4)write(ifch,*)'remnant masses:',am2t,amasini,amasmin
     &                                          ,mdrmax

c.............................exotic ...................................

c      if(amasini.gt.amasmin.and.irmdropx.eq.1)then

c      if((iept.eq.6.or.
c     &   .not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
      if(.not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
     &           .or.(nqu.eq.1.and.naq.eq.1))
     &    .and.amasini.gt.amasmin.and.irmdropx.eq.1)then

c      if((
c     &   .not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
c     &           .or.(nqu.eq.1.and.naq.eq.1)).or.
c     &   (iept.ne.0.and.iept.le.2.and.reminv/ept(5).gt.rangen()))
c     &    .and.amasini.gt.amasmin.and.irmdropx.eq.1)then

         !print*,'-------------------------------------------' !!!
         !print*,jcf
         !print*,icf,sqrt(amasini),sqrt(amasmin),sqrt(mdrmax)  !!!
         !print*,nqu,naq                                      !!!
        if(amasini.gt.mdrmax.or.(jcf(4,1)+jcf(4,2).ne.0))then         !charm not possible in droplet
          call getdroplet(ir,icf,jcf,ept,aa,gdrop,mdrmax)
          !--------------------------------
          !emit a droplet, update the remnant string flavour and 5-momentum 
          ! input
          !     ir ......... 1  projectile, -1  target remnant
          !     ept ........ remnant  5-momentum 
          !     jcf ........ remnant jc
          ! output
          !     gdrop ...  .true. = successful droplet emission
          !                          icf, ept ....... droplet  ic and 5-momentum 
          !                          jcf, a ......... remnant string jc and 5-momentum
          !               .false. = unsuccessful
          !                          jcf, ept .... unchanged, 
          !                          emits hadrons instead of droplet
c          !                          considered as droplet jc and 5-momentum
          !-------------------------------------
          if(.not.gdrop)goto 500
        endif 

        !...........droplet
        !also in case of unsuccessful drop emission, then remnant = droplet ! 
        idrop=1
        nptl=nptl+1   
        t=xorptl(4,mm)
        istptl(mm)=41
        ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        tivptl(2,mm)=t
c            Remnant radius to have eps=dens GeV/fm3
        radptl(nptl)=(3.*sngl(ept(5))/4./pi/dens)**0.3333
        dezptl(nptl)=0.
        do l=1,5
          pptl(l,nptl)=sngl(ept(l))
        enddo
        idx=idtra(icf,0,0,3)
        if(idx.ne.0)then
         amx=sngl(ept(5))
         call idres(idx,amx,idrx,iadjx,1)
         idx=idrx
        endif
        if(idx.eq.0)then
          istptl(nptl)=10
          idptl(nptl)=8*10**8+icf(1)*100+icf(2)/100
          if(gproj)then
            ityptl(nptl)=40
          else  !gtarg
            ityptl(nptl)=50
          endif
        else  
          istptl(nptl)=0
          idptl(nptl)=idx
          pptl(5,nptl)=amx
          if(gproj)then
            ityptl(nptl)=45
            if(iept.eq.6)ityptl(nptl)=47
          else  !gtarg
            ityptl(nptl)=55
            if(iept.eq.6)ityptl(nptl)=57
          endif
        endif  
        iorptl(nptl)=mm
        jorptl(nptl)=0 
        ifrptl(1,nptl)=0 
        ifrptl(2,nptl)=0 
        xorptl(1,nptl)=xorptl(1,mm)
        xorptl(2,nptl)=xorptl(2,mm)
        xorptl(3,nptl)=xorptl(3,mm)
        xorptl(4,nptl)=t
        tivptl(1,nptl)=t
        call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
        tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
        do l=1,4
          ibptl(l,nptl)=0
        enddo
        andropl=andropl+1
        if(ish.ge.3)write(ifch,*)'Proref,ept(5),id',ept(5),idptl(nptl)
        !print*,nptl,idptl(nptl),sngl(ept(5)),pptl(5,nptl)  !!!

        !..........remnant update
        if(gdrop)then  !drop emission: new remnant -> ept, icf
          idrop=0
          do l=1,5
            ept(l)=aa(l)
          enddo
          call idquacjc(jcf,nqu,naq)
          call idenco(jcf,icf,iret)
          if(iret.eq.1)call utstop('Pb in ProRef in strg+drop process&')
          !!!  print*,'new remnant:',icf,ept(5)    !!!
          nptl=nptl+1   
          t=xorptl(4,mm)
          ifrptl(2,mm)=nptl
          do l=1,5
            pptl(l,nptl)=sngl(ept(l))
          enddo
          idptl(nptl)=idptl(mm)
          istptl(nptl)=40
          iorptl(nptl)=mm
          jorptl(nptl)=0 
          ifrptl(1,nptl)=0 
          ifrptl(2,nptl)=0 
          xorptl(1,nptl)=xorptl(1,mm)
          xorptl(2,nptl)=xorptl(2,mm)
          xorptl(3,nptl)=xorptl(3,mm)
          xorptl(4,nptl)=t
          tivptl(1,nptl)=t
          tivptl(2,nptl)=ainfin
          if(gproj)then
            ityptl(nptl)=40
          else   !gtarg
            ityptl(nptl)=50
          endif
          do l=1,4
            ibptl(l,nptl)=0
          enddo
        endif

        !........decay mini-droplet......
        mm=nptlini+1
        nptlb=nptl
        if(iabs(idptl(mm)).gt.10**8)then
          if(ish.ge.3)write(ifch,*)'Make droplet'
          if(nptlb.gt.mxptl-10)call utstop('ProRef: mxptl too small&')
          iret=0
          if(ifrade.gt.0.and.ispherio.eq.0)call hnbaaa(mm,iret)!Decay remn
          if(iret.ne.1.and.nptl.ne.nptlb)then ! ---successful decay---
            istptl(mm)=istptl(mm)+1
            ifrptl(1,mm)=nptlb+1
            ifrptl(2,mm)=nptl
            t=tivptl(2,mm)
            x=xorptl(1,mm)+(t-xorptl(4,mm))*pptl(1,mm)/pptl(4,mm)
            y=xorptl(2,mm)+(t-xorptl(4,mm))*pptl(2,mm)/pptl(4,mm)
            z=xorptl(3,mm)+(t-xorptl(4,mm))*pptl(3,mm)/pptl(4,mm)
            do 21 n=nptlb+1,nptl
              iorptl(n)=mm
              jorptl(n)=0
              istptl(n)=0
              ifrptl(1,n)=0
              ifrptl(2,n)=0
              if(idfptl(mm).eq.0)then
                idfptl(n)=0
              else
                idfptl(n)=1
              endif
              radius=0.8*sqrt(rangen())
              phi=2*pi*rangen()
              ti=t
              zi=z
              xorptl(1,n)=x + radius*cos(phi)
              xorptl(2,n)=y + radius*sin(phi)
              xorptl(3,n)=zi
              xorptl(4,n)=ti
              iioo=mm
              zor=dble(xorptl(3,iioo))
              tor=dble(xorptl(4,iioo))
              call idquac(iioo,nq,ndummy1,ndummy2,jcdummy)
              r=rangen()
              tauran=-taurea*alog(r)
              call jtaix(n,tauran,zor,tor,zis,tis)
              tivptl(1,n)=amax1(ti,tis)
              call idtau(idptl(n),pptl(4,n),pptl(5,n),taugm)
              r=rangen()
              tivptl(2,n)=t+taugm*(-alog(r))
              if(gproj)then
                ityptl(n)=41
                if(iept.eq.6)ityptl(n)=47
              else  !gtarg
                ityptl(n)=51
                if(iept.eq.6)ityptl(n)=57
              endif
              radptl(n)=0.
              dezptl(n)=0.
              itsptl(n)=0
              rinptl(nptl)=-9999
   21       continue
            if(iabs(idptl(nptlb+1)).le.6) then 
              call gakli2(0,0)
              if(ish.ge.1)write (ifmt,*)'string from drop:nptlb+1,nptl:'
     *                                 ,nptlb+1,nptl
              istptl(nptlb+1)=1
              do n=nptlb+2,nptl
                istptl(n)=20
                zpaptl(n)=0.
              enddo
              call gakfra(iret)
              call gakli2(0,0)
            endif
            jerr(4)=jerr(4)+1
          elseif(ifrade.gt.0.and.ispherio.eq.0)then ! Unsuccessful decay
            jerr(5)=jerr(5)+1
            if(ish.ge.4)write(ifch,*)
     *         '***** Unsuccessful remnant cluster decay'
     *             ,' --> do RemoveHadrons instead.'
            mm=mmini
            nptl=nptlini
            irmdropx=0
            goto 55
          endif
        endif
        
        if(idrop.eq.1)goto 1000   
        !successful drop decay, no additional string, nothing to do 

      endif

c...............................................................

 500  mm=mmini
      if(gdrop)mm=nptlini+2
      istptl(mm)=41
      ifrptl(1,mm)=nptl+1

c........................remove hadrons.........................
      
      nbar=0
      nmes=0

      if(.not.((nqu.eq.3.and.naq.eq.0).or.(nqu.eq.0.and.naq.eq.3)
     &          .or.(nqu.eq.1.and.naq.eq.1)))then
        if(irmdropx.eq.irmdrop)then
          jerr(6)=jerr(6)+1
             !call utmsg('ProReF')
             !write(ifch,*)'***** condition for droplet treatment: '
             !write(ifch,*)'*****  amasini.gt.amasmin.and.irmdropx.eq.1 = '        
             !*           ,amasini.gt.amasmin.and.irmdropx.eq.1
             !write(ifch,*)'***** amasini,amasmin,irmdropx:'
             !*                 ,amasini,amasmin,irmdropx
             !write(ifch,*)'***** nqu,naq:',nqu,naq 
             !write(ifch,*)'***** call RemoveHadrons'
             !call utmsgf
        endif
        call RemoveHadrons(gproj,gtarg,ghadr,m,mm,jcf,icf,ept)
      endif

c........................ determine idr (0=string, else=resonance).......

      if(icf(1).eq.0.and.icf(2).eq.0)then
        id=110
      else
        id=idtra(icf,0,0,3)
      endif  
      idr=0
      am=sngl(ept(5))
      call idres(id,am,idr,iadj,1)
      if(iabs(mod(idr,10)).le.2.and.idr.ne.0)then
       id=idr
      else
       idr=0
      endif                                !ckeck on-shell mass (see uti)
      if(iadj.ne.0.and.iept.gt.0.and.ept(5).gt.0.d0
     &     .and.(dabs((ept(4)+ept(3))*(ept(4)-ept(3))          
     $           -ept(2)**2-ept(1)**2-dble(am)**2).gt.0.3d0))idr=0

      if(ish.ge.3)then
        write(ifch,'(a,5e11.3)')' updt:',(sngl(ept(k)) ,k=1,5)
        write(ifch,*)'            icf: ',icf,' idr: ',idr,' iept: ',iept
      endif

      if(iept.eq.3)stop'ProReF: iept=3 ???'

c...........................................string...................
      if(iept.gt.0.and.idr.eq.0)then 

        !... nqu of remainder string

        anstrg0=anstrg0+1
        if(gdrop)anstrg1=anstrg1+1

        call iddeco(icf,jcf)
        nqu=0
        do l=1,4
          nqu=nqu+jcf(l,1)-jcf(l,2)
        enddo

        if(zbarfl.lt.0.)stop'ProReF: not supported any more.         '

        !......determine forward momentum ep   
          
           !ptt=0.5*min(zopmax,zopinc*zz)
           !phi=2.*pi*rangen()
           !ptt1=dble(ptt*cos(phi))
           !ptt2=dble(ptt*sin(phi))

        ep(1)=0
        ep(2)=0
        ep(3)=ir*0.5d0*ept(5)
        ep(4)=   0.5d0*ept(5)

        call utlob2(-1,ept(1),ept(2),ept(3),ept(4),ept(5)
     *     ,ep(1),ep(2),ep(3),ep(4),25)

        !....determine forward and backward flavor icf, icb
 
        ireminv=0
c        if(iept.le.2.and.ept(5)/reminv.lt.rangen())ireminv=1
c        if(iept.eq.2.and.ept(5).lt.reminv.and.rangen().lt.0.5)ireminv=1
        if(iept.eq.6.and.rangen().lt.0.25)ireminv=1
c        if(iept.le.2)then
c          if(reminv/ept(5).gt.rangen())ireminv=1
c        elseif(iept.eq.6)then
c          ireminv=1
c        endif
        if(nqu.eq.3)then      !---baryon---
          iq=idrafl(iclpt,jcf,1,'v',iret)
          call uticpl(icf,iq,2,iret)       ! antiquark
          call uticpl(icb,iq,1,iret)       ! quark
          if(ireminv.eq.1)then
           iq=idrafl(iclpt,jcf,1,'v',iret)
           call uticpl(icf,iq,2,iret)       ! antiquark
           call uticpl(icb,iq,1,iret)       ! quark
          endif
        elseif(nqu.eq.-3)then !---antibaryon---
          iq=idrafl(iclpt,jcf,2,'v',iret)
          call uticpl(icf,iq,1,iret)       ! quark
          call uticpl(icb,iq,2,iret)       ! antiquark
          if(ireminv.eq.1)then
           iq=idrafl(iclpt,jcf,2,'v',iret)
           call uticpl(icf,iq,1,iret)       ! quark
           call uticpl(icb,iq,2,iret)       ! antiquark
          endif
        elseif(nqu.eq.0)then !---meson---
           iq1=idrafl(iclpt,jcf,1,'v',iret)
           iq2=idrafl(iclpt,jcf,2,'v',iret)
           if(rangen().gt.0.5)then
             call uticpl(icf,iq1,2,iret) ! subtract quark
             call uticpl(icb,iq1,1,iret) ! add quark
           else
             call uticpl(icf,iq2,1,iret) ! subtract antiquark
             call uticpl(icb,iq2,2,iret) ! add antiquark
           endif
c        elseif(nqu.eq.0)then !---meson---
c          if(iept.ne.1.and.iept.ne.6.and.rangen().lt.0.5)then
c           iq=idrafl(iclpt,jcf,1,'v',iret)
c           call uticpl(icf,iq,2,iret)       ! subtract quark
c           call uticpl(icb,iq,1,iret)       ! add quark
c          else
cc put quark in forward direction always for inelastic
c           iq=idrafl(iclpt,jcf,2,'v',iret)
c           call uticpl(icf,iq,1,iret)       ! subtract antiquark
c           call uticpl(icb,iq,2,iret)       ! add antiquark
c          endif
        else
          call utmsg('ProReF')
          write(ifch,*)'***** neither baryon nor antibaryon nor meson.'
          write(ifch,*)'*****  number of net quarks:',nqu
          call utstop('ProRef&') 
        endif
c        if(nqu.eq.3)then      !---baryon---
c          iq1=idrafl(iclpt,jcf,1,'v',iret)
c          iq2=idrafl(iclpt,jcf,1,'v',iret)
c          iq3=idrafl(iclpt,jcf,1,'v',iret)
c          amdqa=qmass(iq2)+qmass(iq3)+qmass(0)
c          if(rangen().lt.qmass(iq1)/amdqa)ireminv=1
c          if(ireminv.ne.1)then
c            call uticpl(icf,iq1,2,iret) ! antiquark
c            call uticpl(icb,iq1,1,iret) ! quark
c          else
c            call uticpl(icf,iq2,2,iret) ! antiquark
c            call uticpl(icb,iq2,1,iret) ! quark
c            call uticpl(icf,iq3,2,iret) ! antiquark
c            call uticpl(icb,iq3,1,iret) ! quark
c          endif
c        elseif(nqu.eq.-3)then !---antibaryon---
c          iq1=idrafl(iclpt,jcf,2,'v',iret)
c          iq2=idrafl(iclpt,jcf,2,'v',iret)
c          iq3=idrafl(iclpt,jcf,2,'v',iret)
c          amdqa=qmass(iq2)+qmass(iq3)+qmass(0)
c          if(rangen().lt.qmass(iq1)/amdqa)ireminv=1
c          if(ireminv.ne.1)then
c            call uticpl(icf,iq1,1,iret) ! antiquark
c            call uticpl(icb,iq1,2,iret) ! quark
c          else
c            call uticpl(icf,iq2,1,iret) ! antiquark
c            call uticpl(icb,iq2,2,iret) ! quark
c            call uticpl(icf,iq3,1,iret) ! antiquark
c            call uticpl(icb,iq3,2,iret) ! quark
c          endif
c        elseif(nqu.eq.0)then !---meson---
c           iq1=idrafl(iclpt,jcf,1,'v',iret)
c           iq2=idrafl(iclpt,jcf,2,'v',iret)
c           if(rangen().lt.qmass(iq1)/qmass(iq2))then
cc           if(rangen().gt.0.5)then
c             call uticpl(icf,iq1,2,iret) ! subtract quark
c             call uticpl(icb,iq1,1,iret) ! add quark
c           else
c             call uticpl(icf,iq2,1,iret) ! subtract antiquark
c             call uticpl(icb,iq2,2,iret) ! add antiquark
c           endif
c        else
c          call utmsg('ProReF')
c          write(ifch,*)'***** neither baryon nor antibaryon nor meson.'
c          write(ifch,*)'*****  number of net quarks:',nqu
c          call utstop('ProRef&') 
c        endif

        !..... forward string end
         
        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        pptl(1,nptl)=sngl(ep(1))
        pptl(2,nptl)=sngl(ep(2))
        pptl(3,nptl)=sngl(ep(3))
        pptl(4,nptl)=sngl(ep(4))
        pptl(5,nptl)=0.
        istptl(nptl)=20
        iorptl(nptl)=mm
        if(.not.gdrop)istptl(mm)=41
        jorptl(nptl)=0
        if(nmes.eq.0.and.nbar.eq.0.and..not.gdrop)ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        xorptl(1,nptl)=xorptl(1,mm)
        xorptl(2,nptl)=xorptl(2,mm)
        xorptl(3,nptl)=xorptl(3,mm)
        xorptl(4,nptl)=xorptl(4,mm)
        tivptl(1,nptl)=xorptl(4,nptl)
        tivptl(2,nptl)=xorptl(4,nptl)
        idptl(nptl)=idtra(icf,0,0,3)
        if(gproj)then
          if(iep(m).lt.1)stop'ProReF: iep(m)<1     '
          ityptl(nptl)=41+iep(m)  ! =42 =43 =47
          if(gdrop.and.iep(m).ne.6)ityptl(nptl)=44
          if(ghadr)ityptl(nptl)=44
        else  !gtarg
          if(iet(m).lt.1)stop'ProReF: iet(m)<1     ' 
          ityptl(nptl)=51+iet(m)  !=52 =53 =57
          if(gdrop.and.iet(m).ne.6)ityptl(nptl)=54
          if(ghadr)ityptl(nptl)=54
        endif  
        itsptl(nptl)=1
        qsqptl(nptl)=0.
        rinptl(nptl)=-9999
        !write(6,'(a,i9,$)')'     ',idptl(nptl) !======================
        if(gproj)then
          zpaptl(nptl)=zz
        else  !gtarg
          zpaptl(nptl)=0
        endif
        if(ish.ge.3)then
          write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
          write(ifch,*)' id: ',idptl(nptl)
        endif
        !....... backward string end

        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        pptl2=0.
        do i=1,3
         pptl(i,nptl)=sngl(ept(i)-ep(i))
         pptl2=pptl2+pptl(i,nptl)*pptl(i,nptl)
        enddo
        pptl(4,nptl)=sqrt(pptl2)
        pptl2=sngl(ept(4)-ep(4))
        if(ish.ge.1.and.abs(pptl2-pptl(4,nptl)).gt.max(0.1,
     &                                         0.1*abs(pptl2)))then
          write(ifmt,*)
     &    'Warning in ProRef: inconsistent backward string end energy !'
     &    ,pptl(4,nptl),pptl2,abs(pptl2-pptl(4,nptl))
          if(ish.ge.2)write(ifch,*)
     &    'Warning in ProRef: inconsistent backward string end energy !'
     &    ,(pptl(kkk,nptl),kkk=1,4),pptl2,abs(pptl2-pptl(4,nptl))
        endif
        pptl(5,nptl)=0.
        istptl(nptl)=20
        iorptl(nptl)=mm
        jorptl(nptl)=0
        ifrptl(2,mm)=nptl
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        xorptl(1,nptl)=xorptl(1,mm)
        xorptl(2,nptl)=xorptl(2,mm)
        xorptl(3,nptl)=xorptl(3,mm)
        xorptl(4,nptl)=xorptl(4,mm)
        tivptl(1,nptl)=xorptl(4,nptl)
        tivptl(2,nptl)=xorptl(4,nptl)
        idptl(nptl)=idtra(icb,0,0,3)
        if(gproj)then
          ityptl(nptl)=41+iep(m)  ! =42 =43 =47
          if(gdrop.and.iep(m).ne.6)ityptl(nptl)=44
          if(ghadr)ityptl(nptl)=44
        else  !gtarg
          ityptl(nptl)=51+iet(m)  !=52 =53 =57
          if(gdrop.and.iep(m).ne.6)ityptl(nptl)=54
          if(ghadr)ityptl(nptl)=54
        endif
        itsptl(nptl)=1
        qsqptl(nptl)=0.
        rinptl(nptl)=-9999
        !write(6,'(a,i9)')'     ',idptl(nptl) 
        if(gtarg)then
          zpaptl(nptl)=zz
        else  !gproj
          zpaptl(nptl)=0
        endif
        if(ish.ge.3)then
          write(ifch,'(a,5e11.3,$)')' kink:',(pptl(k,nptl),k=1,5)
          write(ifch,*)' id: ',idptl(nptl)
        endif
         
c............................no string = resonance...................
      else 
        
        anreso0=anreso0+1
        if(gdrop)anreso1=anreso1+1

        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('ProRef: mxptl too small&')
        if(iept.eq.0)call idmass(id,am)
        idptl(nptl)=id
        pptl(1,nptl)=sngl(ept(1))
        pptl(2,nptl)=sngl(ept(2))
        pptl(3,nptl)=sngl(ept(3))
        pptl(4,nptl)=sngl(ept(4))
        pptl(5,nptl)=am 
        istptl(nptl)=0
        iorptl(nptl)=mm
        if(.not.gdrop)istptl(mm)=41
        jorptl(nptl)=0
        if(nmes.eq.0.and.nbar.eq.0.and..not.gdrop)ifrptl(1,mm)=nptl
        ifrptl(2,mm)=nptl
        ifrptl(1,nptl)=0
        ifrptl(2,nptl)=0
        xorptl(1,nptl)=xorptl(1,mm)
        xorptl(2,nptl)=xorptl(2,mm)
        xorptl(3,nptl)=xorptl(3,mm)
        xorptl(4,nptl)=xorptl(4,mm)
        tivptl(1,nptl)=xorptl(4,nptl)
        call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
        tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
        if(gproj)then
          ityptl(nptl)=45
          if(gdrop)ityptl(nptl)=46
          if(ghadr)ityptl(nptl)=46
          if(iept.eq.6)ityptl(nptl)=47
          if(iept.eq.2)ityptl(nptl)=48
        else   !gtarg
          ityptl(nptl)=55
          if(gdrop)ityptl(nptl)=56
          if(ghadr)ityptl(nptl)=56
          if(iept.eq.6)ityptl(nptl)=57
          if(iept.eq.2)ityptl(nptl)=58
        endif
        itsptl(nptl)=0
        qsqptl(nptl)=0.
        rinptl(nptl)=-9999

        if(ish.ge.3)write(ifch,'(a,5e10.3,i7)')' nucl:'
     *         ,(pptl(i,nptl),i=1,5),idptl(nptl)

      endif
c.......................................................................      
c      print *,iep(1),iet(1),ityptl(nptl)       
 1000 call utprix('ProReF',ish,ishini,3)
ctp060829        if(ityptl(nptl).gt.60)print*,ityptl(nptl)
      return

      end

c---------------------------------------------------------------------------------------
      subroutine RemoveHadrons(gproj,gtarg,ghadr,m,mm,jcf,icf,ept)
c---------------------------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      integer jcf(nflav,2),icf(2)
      double precision aa(5),ept(5)  
      logical ghadr,gproj,gtarg
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)

      if(gproj)then
        ir=1
      elseif(gtarg)then
        ir=-1
      else
        call utstop('RemoveHadron : neither proj or targ !&')
      endif
      
      call idquacjc(jcf,nqu,naq)
      if(nqu.eq.naq)then
        nmes=nqu-1
        nbar=0
      elseif(nqu.gt.naq)then
        nmes=naq
        nbar=(nqu-naq-3)/3     !nbar baryons
      else
        nmes=nqu
        nbar=(naq-nqu-3)/3    !nbar antibaryons
      endif
      if(nmes+nbar.gt.0)ghadr=.true.
c  remove mesons 
       if(nmes.gt.0)then
          do mes=1,nmes
            !write(ifch,*)'remove meson',mes,' / ',nmes
            call gethadron(1,idd,aa,jcf,ept,ir,iret)
            call idenco(jcf,icf,iret2)  
            if(iret.eq.0.and.iret2.eq.0)then
              nptl=nptl+1
              if(nptl.gt.mxptl)
     &             call utstop('RemoveHadrons: mxptl too small&')
              idptl(nptl)=idd
              do i=1,5
                pptl(i,nptl)=sngl(aa(i))
              enddo           
              iorptl(nptl)=mm
              jorptl(nptl)=0
              if(mes.eq.1)then
                ifrptl(1,mm)=nptl
                ifrptl(2,mm)=nptl
              else
                ifrptl(2,mm)=nptl
              endif
              ifrptl(1,nptl)=0
              ifrptl(2,nptl)=0
              istptl(nptl)=0
              if(gproj)then
                ityptl(nptl)=49
                xorptl(1,nptl)=xproj(m)
                xorptl(2,nptl)=yproj(m)
                xorptl(3,nptl)=zproj(m)
              elseif(gtarg)then
                ityptl(nptl)=59
                xorptl(1,nptl)=xtarg(m)
                xorptl(2,nptl)=ytarg(m)
                xorptl(3,nptl)=ztarg(m)
              endif
              xorptl(4,nptl)=xorptl(4,mm)
              tivptl(1,nptl)=xorptl(4,nptl)
              call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
              tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
              qsqptl(nptl)=0.
            endif
c           deleted: after abstracting a meson, 
c           check if the NEW remnant is a H-Dibaryon
          enddo           
        endif                           
c remove (anti)baryons 
        call idquacjc(jcf,nqu,naq)
        if(nbar.gt.0)then    
          do nb=1,nbar         
            !write(ifch,*)'remove baryon',nb,' / ',nbar
            if(nqu.gt.0)then            
              call gethadron(2,idd,aa,jcf,ept,ir,iret)
            else
              call gethadron(3,idd,aa,jcf,ept,ir,iret)
            endif               
            call idenco(jcf,icf,iret2)
            if(iret.eq.0.and.iret2.eq.0)then
              nptl=nptl+1
              if(nptl.gt.mxptl)
     &             call utstop('RemoveHadron: mxptl too small&')
              idptl(nptl)=idd
              do i=1,5
                pptl(i,nptl)=sngl(aa(i))
              enddo
              iorptl(nptl)=mm
              jorptl(nptl)=0
              if(nmes.eq.0.and.nb.eq.1)then
                ifrptl(1,mm)=nptl
                ifrptl(2,mm)=nptl
              else
                ifrptl(2,mm)=nptl
              endif
              ifrptl(1,nptl)=0
              ifrptl(2,nptl)=0
              istptl(nptl)=0
              if(gproj)then
                ityptl(nptl)=49
                xorptl(1,nptl)=xproj(m)
                xorptl(2,nptl)=yproj(m)
                xorptl(3,nptl)=zproj(m)
              elseif(gtarg)then
                ityptl(nptl)=59
                xorptl(1,nptl)=xtarg(m)
                xorptl(2,nptl)=ytarg(m)
                xorptl(3,nptl)=ztarg(m)
              endif
              xorptl(4,nptl)=xorptl(4,mm)
              tivptl(1,nptl)=xorptl(4,nptl)
              call idtau(idptl(nptl),pptl(4,nptl),pptl(5,nptl),taugm)
              tivptl(2,nptl)=tivptl(1,nptl)+taugm*(-alog(rangen()))
              qsqptl(nptl)=0.
            endif
c             deleted: after abstracting a (anti)baryon, 
c                                  check if the NEW remnant is a H-Dibaryon
          enddo                  
        endif
      end
      
c------------------------------------------------------------------
         subroutine gethadron(imb,idf,a,jc,ep,ir,iret)
c------------------------------------------------------------------  
c       goal:  emit a hadron (imb= 1 meson, 2 baryon, 3 antibaryon)
c              update the remnant flavour and 5-momentum 
c
c       idf ,a : hadron id and 5-momentum 
c       ir     : 1  projectile, -1  target remnant
c       jc, ep : remnant flavor and 5-momentum
c       iret   : in case of error, keep correct momentum in remnant
c                and lose the quarks of the (not) emitted hadron
c-----------------------------------------------------------------

        include 'epos.inc'
        include 'epos.incems'
        common/cems5/plc,s
        double precision s,plc
        double precision ep(5),a(5),re(5),p1(5)
        integer jc(nflav,2),ifh(3)!,ic(2)
        common /ems12/iodiba,bidiba  ! defaut iodiba=0. if iodiba=1, study H-Dibaryon 
        double precision ptm,qcm,u(3),utpcmd,ptt,phi,sxini,strmas
     &                  ,ampt2dro,ampt2str,p5sq,amasex,drangen

        call utpri('gethad',ish,ishini,5)

        iret=0
        do i=1,5
          a(i)=0.d0
          re(i)=ep(i)
        enddo
    
         if(ish.ge.5)then          
           write(ifch,*)'remnant flavour and 5-momentum:',jc, ep,  ir          
         endif       
         !write(*,'(/a,5f8.3)')'p before: ',ep
         
         if(ir.eq.1)then
           iclpt=iclpro
         else
           iclpt=icltar
         endif
  
c  get the id and mass of hadron, the remnant jc is updated
         
          if(imb.eq.1)then              ! a meson 
            ifq=idrafl(iclpt,jc,1,'v',iret)
            ifa=idrafl(iclpt,jc,2,'v',iret)
            if(ifq.le.ifa)then
                     idf=ifq*100+ifa*10
            else
              idf=-(ifq*10+ifa*100)
            endif            
            call idmass(idf,amss)    
            
          elseif(imb.eq.2)then            ! a baryon                       
            do ik=1,3
              ifh(ik)=idrafl(iclpt,jc,1,'v',iret)              
            enddo
            call neworder(ifh(1),ifh(2),ifh(3))
            idf=ifh(1)*1000+ifh(2)*100+ifh(3)*10
            if(ifh(1).ne.ifh(2).and.ifh(2).ne.ifh(3)
     $        .and.ifh(1).ne.ifh(3))  idf=2130  
            if(ifh(1).eq.ifh(2).and.ifh(2).eq.ifh(3))idf=idf+1
            call idmass(idf,amss)        
     
          elseif(imb.eq.3)then           ! an antibaryon                       
            do ik=1,3
              ifh(ik)=idrafl(iclpt,jc,2,'v',iret)
            enddo 
            call neworder(ifh(1),ifh(2),ifh(3))
            idf=ifh(1)*1000+ifh(2)*100+ifh(3)*10
            if(ifh(1).ne.ifh(2).and.ifh(2).ne.ifh(3)
     $        .and.ifh(1).ne.ifh(3))  idf=2130  
            if(ifh(1).eq.ifh(2).and.ifh(2).eq.ifh(3))idf=idf+1
            idf=-idf
            call idmass(idf,amss)               
          endif                      
             
          if(iret.ne.0)call utstop('Not enough quark in gethad ???&')

c boost remnant in rest frame
      if(ish.ge.6) write (ifch,*) 'on-shell check'
        do k=1,5
          p1(k)=ep(k)
        enddo
        p1(5)=(p1(4)-p1(3))*(p1(4)+p1(3))-p1(2)**2-p1(1)**2
        if(p1(5).gt.0d0.and.abs(p1(5)-ep(5)*ep(5)).lt.ep(5))then
          p1(5)=sqrt(p1(5))
        else
          if(ish.ge.1)write(ifch,*)'Precision problem in gethad, p:',
     &             (p1(k),k=1,5),ep(5)*ep(5)
          p1(5)=ep(5)
          p1(4)=sqrt(p1(3)*p1(3)+p1(2)*p1(2)+p1(1)*p1(1)+p1(5)*p1(5))
        endif

c initial limits

          ptm=p1(5)
          amasex=dble(amss)
          strmas=dble(2.*utamnz(jc,4))
          sxini=ptm*ptm   
c redo

       nredo=0  
 777   continue
       nredo=nredo+1
       if(nredo.gt.20)then
          !write(ifch,*)'nredo.gt.20 -> only drop'
         if(ish.ge.4)write(ifch,*)
     &  'Pb with hadron momentum in Gethad !'
         iret=1
         goto 999    
       endif

c fix pt

          ptt=dble(ranpt()*alpdro(2))**2         !pt+pl
          if(ptt.ge.sxini)goto 777
          sxini=sqrt(sxini-ptt)


          ampt2dro=amasex**2d0
          ampt2str=strmas**2d0

          a(5)=amasex
          re(5)=sxini-a(5)
          if(re(5).le.strmas)then
            if(ish.ge.6)write(ifch,*)
     &           'Pb with initial mass in Gethad, retry',ir
     &       ,amasex,strmas,sxini,ptm,ptt
            goto 777
          endif


c two body decay
          if(ish.ge.6)write(ifch,*)'2 body decay',ptm,a(5),re(5)
          qcm=utpcmd(ptm,a(5),re(5))
          u(3)=2.d0*drangen(qcm)-1.d0
          phi=2.d0*dble(pi)*drangen(u(3))     
          u(1)=sqrt(1.d0-u(3)**2)*cos(phi)
          u(2)=sqrt(1.d0-u(3)**2)*sin(phi)
          if(u(3).ge.0d0)then          !send always hadron backward
            do j=1,3
              re(j)=qcm*u(j)
              a(j)=-re(j)
            enddo
          else
            do j=1,3
              a(j)=qcm*u(j)
              re(j)=-a(j)
            enddo
          endif

          re(4)=sqrt(qcm**2+re(5)**2)
          a(4)=sqrt(qcm**2+a(5)**2)

          if(ish.ge.6)write(ifch,*)'momentum in rest frame : ',re,a


c Fix re of remnant

c boost string in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,re(1),re(2),re(3),re(4),81)

         p5sq=(re(4)+re(3))*(re(4)-re(3))-(re(1)*re(1)+re(2)*re(2))
         if(p5sq.gt.ampt2str)then
           re(5)=sqrt(p5sq) 
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with remnant mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2str
             write(ifch,*)'   momentum four vector:',(re(ii),ii=1,4)
           endif
           goto 777
         endif 

c Fix a of hadron

c boost droplet in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,a(1),a(2),a(3),a(4),82)

         p5sq=(a(4)+a(3))*(a(4)-a(3))-(a(1)**2.d0+a(2)**2.d0)
         if(abs(p5sq-ampt2dro).le.0.1)then
           a(5)=sqrt(p5sq) 
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with hadron mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2dro
             write(ifch,*)'   momentum four vector:',(a(ii),ii=1,4)
           endif
           goto 777
         endif 


 999    continue


        if(iret.eq.1)then      !If problem with momenta do not update remnant

          if(ish.ge.4)
     *    write(ifch,*)'no hadron emission in gethad'

        else     !update the 3-momentum and energy of remnant: ep

          if(ish.ge.1.and.abs(ep(4)-re(4)-a(4)).gt.1.e-2*ep(4))then
            write(ifmt,*)'Pb with energy conservation in gethad'
            if(ish.ge.6)then
              write(ifch,*)'Pb with energy conservation :'
              write(ifch,*)'   p1_ini:',ep(1),'  p1:',re(1)+a(1)
              write(ifch,*)'   p2_ini:',ep(2),'  p2:',re(2)+a(2)
              write(ifch,*)'   p3_ini:',ep(3),'  p3:',re(3)+a(3)
            endif
          endif
          
          do i=1,5
            ep(i)=re(i)
          enddo
          if(ish.ge.5)then
            write(ifch,*)'get hadron with id and 5-momentum:',idf, a
          endif       
        
        endif
        !do i=1,5
        !  sm(i)=ep(i)+a(i)
        !enddo
        !write(*,'(a,5f8.3,i5)')'p after:  ',sm,iret

c      ghost condition
c         if(abs((a(4)+a(3))*(a(4)-a(3))
c     $           -a(2)**2-a(1)**2-a(5)**2).gt.0.3
c     $      .and.  abs(1.-abs(a(3))/a(4)).gt.0.01)print*,iret,dd

c$$$        if(iodiba.eq.1)then  ! for H-dibaryon study ??????????
c$$$          call idenco(jc,ic,iret) 
c$$$          if(ic(1).eq.222000.and.ic(2).eq.0)ep(5)=ep(5)-bidiba
c$$$        endif  
         
        if(ish.ge.5)then
          write(ifch,*)'new remnant flavour and 5-momentum:',jc, ep 
        endif       
c          write(ifmt,*)'get hadron with id and 5-momentum:',idf, a
c          write(ifmt,*)'new remnant flavour and 5-momentum:',jc, ep 
        
        call utprix('gethad',ish,ishini,5)
        
      return
      end



c------------------------------------------------------------------
         subroutine getdroplet(ir,ic,jc,ep,a,pass,mdrmax)
c------------------------------------------------------------------  
c  emit a droplet, update the remnant string flavour and 5-momentum 
c
c input
c       ir ........ 1  projectile, -1  target remnant
c       ep ........ remnant  5-momentum 
c       jc ........ remnant jc
c output
c       pass ...  .true. = successful droplet emission
c                            ic, ep ....... droplet  ic and 5-momentum 
c                            jc, a ........ remnant string jc and 5-momentum
c                 .false. = unsuccessful
c                            jc, ep .... unchanged, 
c                            considered as droplet jc and 5-momentum
c-----------------------------------------------------------------

        include 'epos.inc'
        include 'epos.incems'
        double precision ep(5),a(5),p1(5),re(5),eps,amasex,mdrmax
        double precision xxx,rr,alp,p5sq,xmin,xmax,ampt2str
     &  ,sxini,strmas,xxxmax,xxxmin,ampt2dro,mdrmaxi
        parameter(eps=1.d-20)
        integer jc(nflav,2),jcini(nflav,2),jcfin(nflav,2),ifh(3),ic(2)
        integer icx(2)   !,icxx(2)
        logical pass
        common/cems5/plc,s
        double precision s,plc,ptm,qcm,u(3),utpcmd,ptt,drangen,phi
      
        call utpri('getdro',ish,ishini,4)
        
        iret=0
        iret2=0
        mdrmaxi=mdrmax
        pass=.true.
        idps=0
        idms=0
        do i=1,nflav
          jcini(i,1)=jc(i,1)
          jcini(i,2)=jc(i,2)
          jcfin(i,1)=0
          jcfin(i,2)=0
        enddo

        call idquacjc(jc,nqu,naq)

        do i=1,5
          a(i)=0.d0
          re(i)=0.d0
        enddo
        npart=nqu+naq
   
         if(ir.eq.1)then
           iclpt=iclpro
         else
           iclpt=icltar
         endif
  
         if(ish.ge.5)then          
           write(ifch,*)'remnant flavour and 5-momentum:',jc,ep,npart 
         endif       
           
c  get id of string ends, the remnant string jc is updated

         if(npart.lt.3.and.ep(5).lt.mdrmax.and.iclpt.ne.4)then !light droplet with few quarks
            pass=.false.
            goto 1000
         elseif(npart.lt.3)then    !few quarks but heavy, add some quarks to extract a q-qbar string (should not exit directly because of the large mass)
           ifq=jdrafl(iclpt,jcini,2,iret2)
           if(nqu.eq.1.and.naq.eq.1)then
             idps=1
             idms=1
             nqu=2
             naq=2
           else
             call utstop('This should not happen (getdrop) !&')
           endif
         elseif((nqu.eq.2.and.naq.le.2).or.(nqu.le.2.and.naq.eq.2))then
           idps=1
           idms=1
         elseif(naq.eq.0)then
           idps=5
           idms=1
         elseif(nqu.eq.0)then
           idps=1
           idms=5
         else                 !There is enough q or aq to do qq-q string
         
           
           if(jcini(4,1)-jcini(4,2).eq.0)then !if c-cbar

             idps=1
             idms=1

           else

c One chooses the first q or aq

           rrr=rangen()
           npart=nqu+naq
           if(jcini(4,1)+jcini(4,2).ne.0)then !if some charm take it out
             if(jcini(4,1).ne.0)then
               idps=1
               nqu=nqu-1
             else
               idms=1
               naq=naq-1
             endif
           elseif(rrr.gt.float(naq)/float(npart))then
             idps=1
             nqu=nqu-1
           else
             idms=1
             naq=naq-1
           endif
           
c One chooses the second one

           rrr=rangen()
           npart=nqu+naq
           if(idps.eq.1.and.jcini(4,1).ne.0)then !if some charm take it out
             idps=5
           elseif(idms.eq.1.and.jcini(4,2).ne.0)then !if some charm take it out
             idms=5
           elseif(rrr.gt.float(naq)/float(npart))then
             if(idps.eq.1.and.nqu.ge.2)then
               idps=5
             else
               idps=1
             endif
           else
             if(idms.eq.1.and.naq.ge.2)then
               idms=5
             else
               idms=1
             endif
           endif

c If there is already 2 q or 2 aq as string end, we know that we need 
c a third one to complete the string           

           if(idps.eq.5)idms=1
           if(idms.eq.5)idps=1
           if(idps.eq.1.and.idms.ne.5)idms=1
           if(idms.eq.1.and.idps.ne.5)idps=1

         endif

         endif

         if(ish.ge.5)then          
           write(ifch,*)'remnant string ends :',idps,idms  
         endif       

          if(idps.ne.5.and.idms.ne.5)then              ! q-aq string 
            if(jcini(4,1).eq.1)then
              ifq=idrafl(iclpt,jcini,1,'c',iret)
            else
              ifq=idrafl(iclpt,jcini,1,'v',iret)
            endif
            if(jcini(4,1).eq.1)then
              ifa=idrafl(iclpt,jcini,2,'c',iret)
            else
              ifa=idrafl(iclpt,jcini,2,'v',iret)
            endif
            jcfin(ifq,1)=1
            jcfin(ifa,2)=1
              
          elseif(idps.eq.5)then                       ! qq-q string
            do ik=1,3
              if(jcini(4,1).ne.0)then
                ifh(ik)=idrafl(iclpt,jcini,1,'c',iret)
              else
                ifh(ik)=idrafl(iclpt,jcini,1,'v',iret)
              endif
              jcfin(ifh(ik),1)=jcfin(ifh(ik),1)+1
            enddo
     
          elseif(idms.eq.5)then                        !aqaq-aq string 
            do ik=1,3
              if(jcini(4,2).ne.0)then
                ifh(ik)=idrafl(iclpt,jcini,2,'c',iret)
              else
                ifh(ik)=idrafl(iclpt,jcini,2,'v',iret)
              endif
              jcfin(ifh(ik),2)=jcfin(ifh(ik),2)+1
            enddo 
          endif                      

          if(iret.ne.0)call utstop('Not enough quark in getdro ???&')
          if(jcini(4,1)+jcini(4,2).ne.0)
     &         call utstop('There is sitll charm quark in getdro???&')

c droplet id
            
         call idenco(jcini,icx,iret)
         if(iret.eq.1)then
           call utstop('Exotic flavor in getdroplet !&')
         endif
         amx=0
         idx=idtra(icx,0,0,3)
         if(idx.ne.0)call idmass(idx,amx)
ccc         print*,idx,amx
         

c boost remnant in rest frame
      if(ish.ge.6) write (ifch,*) 'on-shell check'
        do k=1,5
          p1(k)=ep(k)
        enddo
        p1(5)=(p1(4)-p1(3))*(p1(4)+p1(3))-p1(2)**2-p1(1)**2
        if(p1(5).gt.0d0.and.abs(p1(5)-ep(5)*ep(5)).lt.ep(5))then
          p1(5)=sqrt(p1(5))
        else
          if(ish.ge.2)write(ifch,*)'Precision problem in getdro, p:',
     &             (p1(k),k=1,5),ep(5)*ep(5)
          p1(5)=ep(5)
          p1(4)=sqrt(p1(3)*p1(3)+p1(2)*p1(2)+p1(1)*p1(1)+p1(5)*p1(5))
        endif
      if(ish.ge.6) write (ifch,*) 'boost vector:',p1

c limits for momenta 

          mamod=4
          mamos=4
          fad=alpdro(1)
          fas=2
          ptm=p1(5)
          amasex=dble(fad*utamnz(jcini,mamod))  
          strmas=dble(fas*utamnz(jcfin,mamos))


c redo

       nredo=0  
 777   continue
       nredo=nredo+1
       if(nredo.eq.10)then
          amasex=dble(utamnz(jcini,mamod))  
          strmas=dble(utamnz(jcfin,mamos))
       elseif(nredo.gt.20)then
          !write(ifch,*)'nredo.gt.20 -> only drop'
         if(ish.ge.4)write(ifch,*)
     &     'Pb with string mass in Getdrop, continue with gethad'
          pass=.false.
         goto 1000       
       endif

c fix pt

          sxini=ptm*ptm    
          ptt=dble(ranpt()*alpdro(2))**2         !pt+pl
          if(ptt.ge.sxini)goto 777
          sxini=sqrt(sxini-ptt)


          ampt2dro=amasex**2d0
          ampt2str=strmas**2d0
          if(ampt2dro.gt.mdrmaxi)then
            mdrmaxi=2d0*ampt2dro
c            write(ifmt,*)'Warning Mmin>Mmax in Getdroplet'
          endif

          xxxmax=min(mdrmaxi,(sxini-strmas)**2)    !strmas/(strmas+ampt2)
          xxxmin=ampt2dro

          if(xxxmin.gt.xxxmax)then
            !write(ifch,*)'Warning Mmin>sxini -> only drop'
           if(ish.ge.4)write(ifch,*)
     &     'Pb with ampt2 in Getdrop, retry',nredo,ir
     &             ,ampt2dro,ampt2str,xxxmin,xxxmax,sxini,ptt,mdrmaxi
            goto 777
          endif



c fix mass

          rr=drangen(xxxmax)
          xmax=xxxmax
          xmin=xxxmin
          alp=dble(alpdro(3))
          if(dabs(alp-1.d0).lt.eps)then
            xxx=xmax**rr*xmin**(1d0-rr)
          else
            xxx=(rr*xmax**(1d0-alp)+(1d0-rr)*xmin**(1d0-alp))
     &                                                **(1d0/(1d0-alp))
          endif

c          write(ifch,*)'ini',xmin,xxx,xmax,rr,ampt2dro
c     &                   ,(sxini-sqrt(xxx)),ampt2str,p1(5)



          re(5)=sqrt(xxx)
          a(5)=sxini-re(5)
          if(a(5).le.strmas)then
            if(ish.ge.6)write(ifch,*)
     &           'Pb with initial mass in Getdrop, retry',ir
     &       ,xmin,xxx,xmax,rr,ampt2dro,ampt2str
            goto 777
          endif


c two body decay
          if(ish.ge.6)write(ifch,*)'2 body decay',ptm,re(5),a(5)
          qcm=utpcmd(ptm,re(5),a(5))
          u(3)=2.d0*drangen(qcm)-1.d0
          phi=2.d0*dble(pi)*drangen(u(3))     
          u(1)=sqrt(1.d0-u(3)**2)*cos(phi)
          u(2)=sqrt(1.d0-u(3)**2)*sin(phi)
          if(u(3).lt.0d0)then          !send always droplet backward
            do j=1,3
              re(j)=qcm*u(j)
              a(j)=-re(j)
            enddo
          else
            do j=1,3
              a(j)=qcm*u(j)
              re(j)=-a(j)
            enddo
          endif

          re(4)=sqrt(qcm**2+re(5)**2)
          a(4)=sqrt(qcm**2+a(5)**2)

          if(ish.ge.6)write(ifch,*)'momentum in rest frame : ',re,a



c Fix a of string

c boost string in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,a(1),a(2),a(3),a(4),71)

         p5sq=(a(4)+a(3))*(a(4)-a(3))-(a(1)**2.d0+a(2)**2.d0)
         if(p5sq.gt.ampt2str)then
           a(5)=sqrt(p5sq) 
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with string mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2str
             write(ifch,*)'   momentum four vector:',(a(ii),ii=1,4)
           endif
           goto 777
         endif 

c Fix ep of droplet

c boost droplet in collision frame
        call utlob2(-1,p1(1),p1(2),p1(3),p1(4),p1(5)
     $       ,re(1),re(2),re(3),re(4),72)

         p5sq=(re(4)+re(3))*(re(4)-re(3))-(re(1)*re(1)+re(2)*re(2))
         if(p5sq.gt.ampt2dro)then
           re(5)=sqrt(p5sq)
         else
           if(ish.ge.6)then
             write(ifch,*)'Pb with droplet mass -> retry'
             write(ifch,*)'   m^2:',p5sq,'  m_min^2:',ampt2dro
             write(ifch,*)'   momentum four vector:',(re(ii),ii=1,4)
           endif
           goto 777
         endif 

         
       if(ish.ge.1.and.abs(ep(4)-re(4)-a(4)).gt.1.e-2*ep(4))then
         write(ifmt,*)'Pb with energy conservation in getdro'
         if(ish.ge.6)then
           write(ifch,*)'Pb with energy conservation :'
           write(ifch,*)'   p1_ini:',ep(1),'  p1:',re(1)+a(1)
           write(ifch,*)'   p2_ini:',ep(2),'  p2:',re(2)+a(2)
           write(ifch,*)'   p3_ini:',ep(3),'  p3:',re(3)+a(3)
         endif
       endif
                               
c If OK, save flavours of droplet and string   
         do i=1,5
           ep(i)=re(i)
         enddo
         ic(1)=icx(1)
         ic(2)=icx(2)
         do i=1,nflav
           jc(i,1)=jcfin(i,1)
           jc(i,2)=jcfin(i,2)
         enddo

         if(ish.ge.6)then
           write(ifch,*)'droplet:'
           write(ifch,*)ic
           write(ifch,*)ep
           write(ifch,*)'string remnant:'
           write(ifch,*)jc 
           write(ifch,*)a 
         endif       
        
 1000    continue
         call utprix('getdro',ish,ishini,4)
         end

c-----------------------------------------------------
       subroutine neworder(n1, n2, n3) 
c-----------------------------------------------------
c make 3 integers ordered like 1 2 3
c------------------------------------------------------
            if(n2.lt.n1)then
              ifb=n2
              n2=n1
              n1=ifb
            endif     
            if(n3.lt.n1)then
              ifb=n3
              n3=n2
              n2=n1
              n1=ifb
            elseif(n3.lt.n2)then   
              ifb=n3
              n3=n2
              n2=ifb
            endif
         end

c-----------------------------------------------------------------------
      function idtr2(ic)
c-----------------------------------------------------------------------
c transforms ic to id such that only hadrons have nonzero id
c-----------------------------------------------------------------------
      parameter (nidt=30)
      integer idt(3,nidt),ic(2)
      data idt/
     * 100000,100000, 110   ,100000,010000, 120   ,010000,010000, 220
     *,100000,001000, 130   ,010000,001000, 230   ,001000,001000, 330
     *,100000,000100, 140   ,010000,000100, 240   ,001000,000100, 340
     *,000100,000100, 440
     *,300000,000000,1111   ,210000,000000,1120   ,120000,000000,1220
     *,030000,000000,2221   ,201000,000000,1130   ,111000,000000,1230
     *,021000,000000,2230   ,102000,000000,1330   ,012000,000000,2330
     *,003000,000000,3331   ,200100,000000,1140   ,110100,000000,1240
     *,020100,000000,2240   ,101100,000000,1340   ,011100,000000,2340
     *,002100,000000,3340   ,100200,000000,1440   ,010200,000000,2440
     *,001200,000000,3440   ,000300,000000,4441/

      idtr2=0
      if(ic(1).eq.0.and.ic(2).eq.0)then
       if(rangen().ge.0.5)then
        idtr2=110
        ic(1)=100000
        ic(2)=100000
       else
        idtr2=220
        ic(1)=10000
        ic(2)=10000
       endif
       return
      endif
      do 1 i=1,nidt
       if(ic(2).eq.idt(1,i).and.ic(1).eq.idt(2,i))idtr2=-idt(3,i)
       if(ic(1).eq.idt(1,i).and.ic(2).eq.idt(2,i))idtr2=idt(3,i)
1     continue
      return
      end

c----------------------------------------------------------------------
      subroutine emsini(e,idpj,idtg)
c----------------------------------------------------------------------
c  energy-momentum sharing initializations
c----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      common/cemsr5/at(0:1,0:6)
      common/cems5/plc,s
      common/cems10/a(0:ntypmx),b(0:ntypmx),d(0:ntypmx)
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      double precision d,a,b,plc,s,amd,dcel,xvpr,xdm,at,xdm2
      common/ems3/dcel,ad
      common/cems13/xvpr(0:3)

c abreviations

      plc=dble(e)
      s=plc**2
      amd=dble(delrex)
      

c alpha (0=0, 1=s, 2=v, 4=d, 8=f)

      a(0)=0d0      
      a(1)=dble(alpsea)    
      a(2)=dble(alpval)    
      a(3)= 0.0d0    
      a(4)=dble(alpdiq)    
      a(5)=dble(a(4))   
      a(6)= 0.0d0    
      a(7)= 0.0d0    
      a(8)=dble(a(2))    
      a(9)= 0.0d0   
       
c beta (0=0, 1=s, 2=v, 4=d, 8=f)

      b(0)=0.0d0     
      b(1)=dble(-alpqua)   
      b(2)=dble(-alpqua)  
      b(3)=0.0d0     
      b(4)=0.0d0     
      b(5)=0.0d0     
      b(6)=0.0d0     
      b(7)=0.0d0     
      b(8)=dble(-alpqua)
      b(9)=0.0d0  
         
      
c alpha_trailing and beta_trailing (0=meson, 1=baryon; 
c                                   0=no excit, 1=nondiffr, 2=diffr, 
c                                   6=nondiffr but no pomeron)

      at(0,0)=0.0d0
      at(0,1)=dble(alpndi)  
      at(0,2)=dble(alpdi)    
      at(0,3)=dble(alpndi)    
      at(0,6)=dble(alpndi)    
      at(1,0)=0.0d0
      at(1,1)=dble(alpndi)  
      at(1,2)=dble(alpdi)    
      at(1,3)=dble(alpndi)    
      at(1,6)=dble(alpndi)    
      
c minimal string masses ( i+j, each one: 0=0, 1=s, 2=v, 4=d, 8=f)

      ammn(0)=0d0
      ammn(1)=0d0
      ammn(2)=dble(ammsqq)+amd   
      ammn(3)=dble(ammsqq)+amd   
      ammn(4)=dble(ammsqq)+amd   
      ammn(5)=dble(ammsqd)+amd   
      ammn(6)=dble(ammsqd)+amd   
      ammn(7)=0d0
      ammn(8)=dble(ammsdd)+amd   
      ammn(9)=dble(ammsqd)+amd   
      ammn(10)=dble(ammsqd)+amd   
      ammn(12)=dble(ammsqd)+amd   
      ammn(16)=0.14d0
      
c minimal pomeron masses (0=0, 1=softPom, 2=regge, 3=hard)

      amprmn(0)=0d0
      amprmn(1)=0d0
      amprmn(2)=0d0
      amprmn(3)=dsqrt(4d0*dble(q2min))
      
c cutoff for virtual pomeron (0=0, 1=soft Pom, 2=regge, 3=hard)

      xvpr(0)=0d0
      xvpr(1)=dble(cumpom**2)/s  
      xvpr(2)=dble(cumpom**2)/s
      xvpr(3)=0.0d0**2/s
      
c minimal remnant masses (0=meson, 1=baryon)

      xdm=0.35d0                  !<pt>
      call idmass(idpj,ampj)
      if(iabs(idpj).gt.1000)then
       ampmn(0)=0.14d0+xdm
       ampmn(1)=dble(ampj)+xdm
      else
       ampmn(0)=dble(ampj)+xdm
       ampmn(1)=0.94d0+xdm
      endif
      call idmass(idtg,amtg)
      if(iabs(idtg).gt.1000)then
       amtmn(0)=0.14d0+xdm
       amtmn(1)=dble(amtg)+xdm
      else
       amtmn(0)=dble(amtg)+xdm
       amtmn(1)=0.94d0+xdm
      endif
      
c minimal excitation masses (0=meson, 1=baryon 
c                            0=no excit, 1=nondiffr, 2=diffr, 
c                                   6=nondiffr but no pomeron)

      xdm2=0.35d0
c to take into account increase of mean pt in inelastic remnants
c      if(isplit.eq.1)xdm=max(2d0*max(0d0,sqrt(log(s))-2.5d0),xdm)
      amemn(0,0)=0d0
      amemn(0,1)=xdm2+0.31d0
c      amemn(0,1)=0d0
      amemn(0,2)=xdm2+0.31d0
c      amemn(0,2)=0.d0
      amemn(0,3)=xdm2+0.31d0
c      amemn(0,6)=xdm2+0.31d0
      amemn(0,6)=0.d0
c      amemn(1,0)=0d0
      amemn(1,1)=xdm2+0.15d0 !+0.15d0
c      amemn(1,1)=0d0
      amemn(1,2)=xdm2+0.15d0
c      amemn(1,2)=0.d0
      amemn(1,3)=xdm2+0.15d0
c      amemn(1,6)=xdm2+0.15d0
      amemn(1,6)=0.d0

c maximal excitation masses (0=no excit, 1=nondiffr, 2=diffr)

      amemx(0)=2d0*xdm
      amemx(1)=plc
      amemx(2)=plc

      if(idpj.gt.1000)then     ! baryon

c initial quark configuration
       ivp0=3
       iap0=0
       idp0=1
       isp0=0
       
      elseif(idpj.lt.-1000)then     ! antibaryon

c initial quark configuration
       ivp0=0
       iap0=3
       idp0=1
       isp0=0
       
      else      ! meson

c initial quark configuration
       ivp0=1
       iap0=1
       idp0=0
       isp0=0
       
      endif
      
      if(idtg.gt.1000)then    ! baryon

c initial quark configuration
       ivt0=3
       iat0=0
       idt0=1
       ist0=0
       
      elseif(idtg.lt.-1000)then   ! antibaryon

c initial quark configuration
       ivt0=0
       iat0=3
       idt0=1
       ist0=0
       
      else       ! meson
      
c initial quark configuration
       ivt0=1
       iat0=1
       idt0=0
       ist0=0
       
      endif


c eikonal parameters

       dcel=dble(chad(iclpro)*chad(icltar)) 

c counters

       antot=0.
       ansh=0.
       ansf=0.
       antotf=0.
       anshf=0.
       ansff=0.
       pp4max=0.
       pp4ini=0.
       andropl=0.
       anstrg0=0.
       anstrg1=0.
       anreso0=0.
       anreso1=0.
       anghadr=0.
       antotre=0.
       anintdiff=0.
       anintsdif=0.
       anintine=0.

      return
      end

c-----------------------------------------------------------------------
      subroutine emsigr
c-----------------------------------------------------------------------
c initialize grid 
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)

      call utpri('emsigr',ish,ishini,5)   
      
      do k=1,koll  !----k-loop---->

c determine length of k-th line of grid

       o=max(1.e-5,min(sngl(om1intc(k)),float(npommx)))!if GFF used for propo
        if(ish.ge.7)write(ifch,*)'emsigr:k,o',k,o
       n=0
       if(o.le.50)then
         p=1./(exp(o)-1)
       else
         p=0.
       endif
10     n=n+1
       p=p*o/n
        if(ish.ge.7)write(ifch,*)'emsigr:n,p',n,p
       if((p.gt.1e-4.or.n.lt.int(o)).and.n.lt.npommx
     *.and.n.lt.nprmax)goto 10

       if(ish.ge.5)write(ifch,*)'emsigr:nmax,b',n,bk(k)

       npr(0,k)=n
       nprmx(k)=n               
       nprt(k)=0
       do i=1,3
        npr(i,k)=0
       enddo
     

c initial value for interaction type

       itpr(k)=0

c initialize grid


       do n=1,nprmx(k)         
        idpr(n,k)=0
        idfpr(n,k)=0
        ivpr(n,k)=1
        nppr(n,k)=0
        nbkpr(n,k)=0
        nvpr(n,k)=0
        idsppr(n,k)=0
        idstpr(n,k)=0
        idrpr(n,k)=0
        idhpr(n,k)=0
        bhpr(n,k)=0.
        xpr(n,k)=0d0
        ypr(n,k)=0d0
        xppr(n,k)=0d0
        xmpr(n,k)=0d0
        xp1pr(n,k)=0d0
        xp2pr(n,k)=0d0
        xm1pr(n,k)=0d0
        xm2pr(n,k)=0d0
        xp1pr(n,k)=0d0
        xp2pr(n,k)=0d0
        xm1pr(n,k)=0d0
        xm2pr(n,k)=0d0
        idp1pr(n,k)=0
        idp2pr(n,k)=0
        idm1pr(n,k)=0
        idm2pr(n,k)=0
        xxp1pr(n,k)=0d0
        xyp1pr(n,k)=0d0
        xxp2pr(n,k)=0d0
        xyp2pr(n,k)=0d0
        xxm1pr(n,k)=0d0
        xym1pr(n,k)=0d0
        xxm2pr(n,k)=0d0
        xym2pr(n,k)=0d0
       enddo
       
      enddo !  <----k-loop-----

      call utprix('emsigr',ish,ishini,5)   
      return
      end

c-----------------------------------------------------------------------
      subroutine emsipt
c-----------------------------------------------------------------------
c initialize projectile and target 
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      common/cems5/plc,s
      common/ems3/dcel,ad
      common/ems6/ivp0,iap0,idp0,isp0,ivt0,iat0,idt0,ist0
      common /cncl/xproj(mamx),yproj(mamx),zproj(mamx)
     *            ,xtarg(mamx),ytarg(mamx),ztarg(mamx)

      double precision dcel,s,plc
      
c initialize projectile

      do i=1,maproj
       idp(i)=idp0
       ivp(i)=ivp0
       iap(i)=iap0
       isp(i)=isp0
       iep(i)=-1
       ifp(i)=0
       kolp(i)=0
       npp(i)=0
       npproj(i)=0
       xxp(i)=0d0
       xyp(i)=0d0
       xpmn(i)=(amemn(idp(i),0)+ampmn(idp(i)))**2/s
       xpmx(i)=dmin1(1d0,(amemx(0)+ampmn(idp(i)))**2/s)
       xpos(i)=0.9d0*(amemx(0)+ampmn(idp(i)))**2/s
       xppmx(i)=0.5d0/(1d0+1d0/dble(maproj)**0.3d0)!1d0-dsqrt(xpmn(i))/maproj
       xmpmx(i)=0.5d0/(1d0+1d0/dble(matarg)**0.3d0)!1d0-dsqrt(xpmn(i))/matarg
       xmpmn(i)=xpmn(i)/xppmx(i)
       xppmn(i)=xpmn(i)/xmpmx(i)
       xpp(i)=1d0
       xmp(i)=0d0
       xppst(i)=0.d0
       xmpst(i)=0.d0
       xposst(i)=0.d0
      enddo

c initialize target

      do j=1,matarg
       idt(j)=idt0
       ivt(j)=ivt0
       iat(j)=iat0
       ist(j)=ist0
       iet(j)=-1
       ift(j)=0
       kolt(j)=0
       npt(j)=0
       nptarg(j)=0
       xxt(j)=0d0
       xyt(j)=0d0
       xtmn(j)=(amemn(idt(j),0)+amtmn(idt(j)))**2/s
       xtmx(j)=dmin1(1d0,(amemx(0)+amtmn(idt(j)))**2/s)
       xtos(j)=0.9d0*(amemx(0)+amtmn(idt(j)))**2/s
       xmtmx(j)=0.5d0/(1d0+1d0/dble(matarg)**0.3d0)!1d0-dsqrt(xtmn(j))/matarg
       xptmx(j)=0.5d0/(1d0+1d0/dble(maproj)**0.3d0)!1d0-dsqrt(xtmn(j))/maproj
       xptmn(j)=xtmn(j)/xmtmx(j)
       xmtmn(j)=xtmn(j)/xptmx(j)
       xmt(j)=1d0
       xpt(j)=0d0
       xmtst(j)=0.d0
       xptst(j)=0.d0
       xtosst(j)=0.d0
      enddo

      return
      end

      
c-----------------------------------------------------------------------
      subroutine emszz
c-----------------------------------------------------------------------      
c     completes /cptl/ for nucleons, checks for no interaction
c     writes   /cevt/
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      common/nucl3/phi,bimp
      common/col3/ncol,kolpt
      integer kolpz(mamx),koltz(mamx)

      call utpri('emszz ',ish,ishini,6)   
 
c     write /cptl/
c     ------------

      if(iokoll.eq.1)then   ! precisely matarg collisions

c nothing to do

      else
 
c determine ncol

       ncolx=ncol
       ncol=0
       ncoli=0
       do 8 k=1,koll
       if(ish.ge.7)write(ifch,*)'k,itpr,ncol,ncolx',k,itpr(k),ncol,ncolx  
        if(itpr(k).eq.0)goto 8
        if(itpr(k).eq.1.or.itpr(k).eq.3)ncoli=ncoli+1
        ncol=ncol+1
        i=iproj(k)
        j=itarg(k)
        istptl(i)=1
        iorptl(i)=-1
        tivptl(2,i)=coord(4,k)
        istptl(maproj+j)=1
        iorptl(maproj+j)=-1
        tivptl(2,maproj+j)=coord(4,k)
8      continue
       if(ncolx.ne.ncol)write(6,*)'ncolx,ncol:', ncolx,ncol 
       if(ncolx.ne.ncol)call utstop('********ncolx.ne.ncol********&')
       if(ncol.eq.0)goto1001
      
c determine npj, ntg

       do ip=1,maproj
        kolpz(ip)=0
       enddo 
       do it=1,matarg
        koltz(it)=0
       enddo 
      do k=1,koll
       if(itpr(k).gt.0)then
        ip=iproj(k)
        it=itarg(k)
        kolpz(ip)=kolpz(ip)+1
        koltz(it)=koltz(it)+1
       endif
      enddo
      npj=0
      do ip=1,maproj
       if(kolpz(ip).gt.0)npj=npj+1
      enddo 
      ntg=0
      do it=1,matarg
       if(koltz(it).gt.0)ntg=ntg+1
      enddo 
c     write(6,*)'npj,ntg,npj+ntg:',npj,ntg,npj+ntg
      
       endif
           
c     write /cevt/
c     ------------

      nevt=1
      bimevt=bimp
      phievt=phi
      kolevt=ncol
      koievt=ncoli
      npjevt=npj
      ntgevt=ntg
      pmxevt=pnll
      egyevt=engy
      !print*,' ===== ',kolevt,koievt,' ====='

c     exit
c     ---- 

      if(ish.ge.7)then
      do n=1,nptl
      write(ifch,115)iorptl(n),jorptl(n),n,istptl(n)
     *,tivptl(1,n),tivptl(2,n)
      enddo
  115 format(1x,'/cptl/',2i6,2i10,2(e10.3,1x))  
      endif

1000  continue
      call utprix('emszz ',ish,ishini,6)
      return

1001  continue 
      if(ish.ge.3)then
      write(ifch,*)
      write(ifch,*)'   ***** no interaction!!!'
      write(ifch,*)'   ***** ncol=0 detected in emszz'
      write(ifch,*)
      endif
      goto1000

      end
      
c-----------------------------------------------------------------------
      subroutine ProCop(i,ii)
c-----------------------------------------------------------------------
c Propose Coordinates of remnants from active projectile nucleons
c-----------------------------------------------------------------------

      include 'epos.incems'
      include 'epos.inc'

      double precision xmptmp,aproj
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icrmn(2)
      double precision s,plc

      nptl=nptl+1
      npproj(i)=nptl
      idptl(nptl)=idptl(ii)*100+99  !100*10**idp(i)+iep(i)
      istptl(nptl)=40
      ityptl(nptl)=40
      iorptl(nptl)=ii
      jorptl(nptl)=0 
      ifrptl(1,nptl)=0 
      ifrptl(2,nptl)=0 

      istptl(ii)=1

c     determine kolz
 
      if(lproj(i).gt.1)then
        zmax=-ainfin
        kolz=0
        do l=1,lproj(i)
          k=kproj(i,l)
          z=coord(3,k)
          if(itpr(k).gt.0.and.z.gt.zmax)then
            zmax=z
            kolz=k
          endif
        enddo
      else
        kolz=1
      endif
c      if(kolz.eq.0)call utstop(' kolz=0 (proj)&')
      if(kolz.eq.0)then
        t=0.
      else
        t=coord(4,kolz)
      endif

      xorptl(1,nptl)=xorptl(1,ii)
      xorptl(2,nptl)=xorptl(2,ii)
      xorptl(3,nptl)=xorptl(3,ii)
      xorptl(4,nptl)=t
      tivptl(1,nptl)=t
      tivptl(2,nptl)=t

      icrmn(1)=icproj(1,i)
      icrmn(2)=icproj(2,i)
      aproj=dble(max(amproj,fremnux2(icrmn)))
c      aprojex=max(ampmn(idp(i))+amemn(idp(i),iep(i))
c     &           ,dble(fremnux(icrmn)))     
      xmptmp=(aproj**2+xxp(i)*xxp(i)+xyp(i)*xyp(i))
     &       /(xpp(i)*s)
      if(iep(i).eq.6)xmptmp=max(xmptmp,xmp(i))
      xpos(i)=1.1d0*xpp(i)*xmptmp
      if(xmptmp.gt.1.d0)then
        xmptmp=0.d0
      if(ish.ge.1)write(ifmt,*)'Warning in ProCop, Remnant mass too low'
      endif

      pptl(1,nptl)=sngl(xxp(i))
      pptl(2,nptl)=sngl(xyp(i))
      pptl(3,nptl)=sngl((xpp(i)-xmptmp)*plc/2d0)
      pptl(4,nptl)=sngl((xpp(i)+xmptmp)*plc/2d0)
      pptl(5,nptl)=amproj

c      write(ifmt,*)'ProCop',i,nptl

      return

      end

c-----------------------------------------------------------------------
      subroutine ProCot(j,jj)
c-----------------------------------------------------------------------
c Propose Coordinates of remnants from active targets nucleons
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      double precision xpttmp,atarg
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      integer icrmn(2)
      double precision s,plc

      nptl=nptl+1
      nptarg(j)=nptl

      idptl(nptl)=idptl(jj)*100+99    !100*10**idt(j)+iet(j)
      istptl(nptl)=40
      ityptl(nptl)=50
      iorptl(nptl)=jj
      jorptl(nptl)=0 
      ifrptl(1,nptl)=0
      ifrptl(2,nptl)=0

      istptl(jj)=1

c     determine kolz
 
      if(ltarg(j).gt.1)then
        zmin=ainfin
        kolz=0
        do l=1,ltarg(j)
          k=ktarg(j,l)
          z=coord(3,k)
          if(itpr(k).gt.0.and.z.lt.zmin)then
            zmin=z
            kolz=k
          endif
        enddo
      else
        kolz=1
      endif
c      if(kolz.eq.0)call utstop(' kolz=0 (targ)&')
      if(kolz.eq.0)then
        t=0.
      else
        t=coord(4,kolz)
      endif

      xorptl(1,nptl)=xorptl(1,jj)
      xorptl(2,nptl)=xorptl(2,jj)
      xorptl(3,nptl)=xorptl(3,jj)
      xorptl(4,nptl)=t
      tivptl(1,nptl)=t
      tivptl(2,nptl)=t

      icrmn(1)=ictarg(1,j)
      icrmn(2)=ictarg(2,j)
      atarg=dble(max(amtarg,fremnux2(icrmn)))
c      atargex=max(amtmn(idt(j))+amemn(idt(j),iet(j))
c     &           ,dble(fremnux(icrmn)))    
      xpttmp=(atarg**2+xxt(j)*xxt(j)+xyt(j)*xyt(j))
     &       /(xmt(j)*s)
      if(iet(j).eq.6)xpttmp=max(xpttmp,xpt(j))
      xtos(j)=1.1d0*xpttmp*xmt(j)
      if(xpttmp.gt.1.d0)then
        xpttmp=0.d0
      if(ish.ge.1)write(ifch,*)'Warning in ProCot, Remnant mass too low'
      endif

      pptl(1,nptl)=sngl(xxt(j))
      pptl(2,nptl)=sngl(xyt(j))
      pptl(3,nptl)=sngl((xpttmp-xmt(j))*plc/2d0)
      pptl(4,nptl)=sngl((xpttmp+xmt(j))*plc/2d0)
      pptl(5,nptl)=amtarg
 
c      write(ifmt,*)'ProCot',j,nptl

      return
      end

c-----------------------------------------------------------------------
      subroutine emswrp(i,ii)
c-----------------------------------------------------------------------

      include 'epos.incems'
      include 'epos.inc'

      double precision p5sq
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,plc
      parameter(eps=1.e-5)

      if(npproj(i).eq.0)then
        write(*,*)'emswrp i ii',i,ii
        call utstop('emswrp with npproj=0 should never happen !&')

c        t=xorptl(4,kolp(i))
c        istptl(ii)=1
c        iorptl(ii)=-1
c        tivptl(2,ii)=t
c        nptl=nptl+1
c        npproj(i)=nptl
c        idptl(nptl)=idptl(ii)*100+99 !100*10**idp(i)+iep(i)
c        istptl(nptl)=40
c        ityptl(nptl)=40
c        iorptl(nptl)=ii
c        jorptl(nptl)=kolp(i)
c        ifrptl(1,nptl)=0 
c        ifrptl(2,nptl)=0 
c        xorptl(1,nptl)=xorptl(1,ii)
c        xorptl(2,nptl)=xorptl(2,ii)
c        xorptl(3,nptl)=xorptl(3,ii)
c        xorptl(4,nptl)=t
c        tivptl(1,nptl)=t
c        tivptl(2,nptl)=t
c        mm=nptl
c        kolp(i)=1
      else
        mm=npproj(i)
      endif
      pptl(1,mm)=sngl(xxp(i))
      pptl(2,mm)=sngl(xyp(i))
      pptl(3,mm)=sngl((xpp(i)-xmp(i))*plc/2d0)
      pptl(4,mm)=sngl((xpp(i)+xmp(i))*plc/2d0)
      if(pptl(4,mm).lt.-eps)call utstop('E pro<0 !&')
      p5sq=xpp(i)*plc*xmp(i)*plc-xxp(i)*xxp(i)-xyp(i)*xyp(i)
      if(p5sq.gt.1.d-10)then
        pptl(5,mm)=sngl(dsqrt(p5sq)) 
      else
        if(ish.ge.2)then
          write(ifch,*)'problem with mass for projectile, '
     &         ,'continue with zero mass'
          write(ifch,*)i,mm,xxp(i),xyp(i),xpp(i),xmp(i),p5sq
        endif
        pptl(5,mm)=0.
      endif 

      do l=1,4
       ibptl(l,mm)=0
      enddo
           
      return

      end

c-----------------------------------------------------------------------
      subroutine emswrt(j,jj)
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      double precision p5sq
      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,plc
      parameter(eps=1.e-5)

      if(nptarg(j).eq.0)then

        write(*,*)'emswrt j jj',j,jj
        call utstop('emswrt with nptarg=0 should never happen !&')

c        t=xorptl(4,kolt(j))
c        istptl(jj)=1
c        iorptl(jj)=-1
c        tivptl(2,jj)=t
c        nptl=nptl+1
c        nptarg(j)=nptl
c        idptl(nptl)=idptl(jj)*100+99 !100*10**idp(i)+iep(i)
c        istptl(nptl)=40
c        ityptl(nptl)=50
c        iorptl(nptl)=jj
c        jorptl(nptl)=kolt(j)
c        ifrptl(1,nptl)=0 
c        ifrptl(2,nptl)=0 
c        xorptl(1,nptl)=xorptl(1,jj)
c        xorptl(2,nptl)=xorptl(2,jj)
c        xorptl(3,nptl)=xorptl(3,jj)
c        xorptl(4,nptl)=t
c        tivptl(1,nptl)=t
c        tivptl(2,nptl)=t
c        mm=nptl
c        kolt(j)=1
      else
        mm=nptarg(j)
      endif
      pptl(1,mm)=sngl(xxt(j))
      pptl(2,mm)=sngl(xyt(j))
      pptl(3,mm)=sngl((xpt(j)-xmt(j))*plc/2d0)
      pptl(4,mm)=sngl((xpt(j)+xmt(j))*plc/2d0)
      if(pptl(4,mm).lt.-eps)call utstop('E targ<0 !&')
      p5sq=xpt(j)*plc*xmt(j)*plc-xxt(j)*xxt(j)-xyt(j)*xyt(j)
      if(p5sq.gt.1.d-10)then
        pptl(5,mm)=sngl(dsqrt(p5sq)) 
      else
        if(ish.ge.2)then
          write(ifch,*)'problem with mass for target, '
     &            ,'continue with zero mass'
          write(ifch,*)j,mm,xxt(j),xyt(j),xpt(j),xmt(j),p5sq
        endif
        pptl(5,mm)=0.
      endif 

      do l=1,4
       ibptl(l,mm)=0
      enddo
 
      return
      end

c-----------------------------------------------------------------------
      subroutine emswrpom(k,i,j)
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      common/cems5/plc,s
      common/emsptl/nppr(npommx,kollmx),npproj(mamx),nptarg(mamx)
      double precision s,px,py,plc

      do 30 n=1,nprmx(k)
       if(idpr(n,k).eq.0.or.ivpr(n,k).eq.0)goto30
       nptl=nptl+1
       nppr(n,k)=nptl
       px=xxp1pr(n,k)+xxp2pr(n,k)+xxm1pr(n,k)+xxm2pr(n,k)
       py=xyp1pr(n,k)+xyp2pr(n,k)+xym1pr(n,k)+xym2pr(n,k)
       pptl(1,nptl)=sngl(px)
       pptl(2,nptl)=sngl(py)
       pptl(3,nptl)=sngl(dsqrt(xpr(n,k))*dsinh(ypr(n,k))*plc)
       pptl(4,nptl)=sngl(dsqrt(xpr(n,k))*dcosh(ypr(n,k))*plc)
       pptl(5,nptl)=sngl(dsqrt(xpr(n,k)*plc*plc-px*px-py*py))
   !    print*,pptl(5,nptl)/plc
       idptl(nptl)=idpr(n,k)*10000   
     &     +idp1pr(n,k)*1000
     &     +idp2pr(n,k)*100
     &     +idm1pr(n,k)*10
     &     +idm2pr(n,k)
       idptl(nptl)=idptl(nptl)*100+99
       istptl(nptl)=30
       iorptl(nptl)=i
       jorptl(nptl)=j
       ifrptl(1,nptl)=0
       ifrptl(2,nptl)=0
       xorptl(1,nptl)=coord(1,k)
       xorptl(2,nptl)=coord(2,k)
       xorptl(3,nptl)=coord(3,k)
       xorptl(4,nptl)=coord(4,k)
       tivptl(1,nptl)=coord(4,k)
       tivptl(2,nptl)=coord(4,k)
       if(idpr(n,k).eq.1)then
        ityptl(nptl)=20
       elseif(idpr(n,k).eq.2)then
        ityptl(nptl)=25
       elseif(idpr(n,k).eq.3)then
        ityptl(nptl)=30
       else
        call utstop('emswrpom: unknown id&')
       endif
       do l = 1,4
        ibptl(l,nptl)=0
       enddo
30    continue

      return
      end

cc--------------------------------------------------------------------------
c      subroutine reaction(idpj,idtg,ireac)   
cc--------------------------------------------------------------------------
cc returns reaction code ireac 
cc--------------------------------------------------------------------------
c      iap=iabs(idpj/10)
c      iat=iabs(idtg/10)
c      isp=idpj/10/iap
c      ist=idtg/10/iat
c      call idchrg(idpj,cp)
c      call idchrg(idtg,ct)
c      ac=abs(cp+ct)
c      if(iap.gt.100)then
c       if(iat.gt.100)then
c        if(isp.eq.1)then
c         if(ist.eq.1)then
c          ireac=1
c         else
c          ireac=6
c         endif
c        else
c         if(ist.eq.1)then
c          ireac=6
c         else
c          ireac=1
c         endif
c        endif
c       elseif(iat.eq.11.or.iat.eq.12.or.iat.eq.22)then
c        if(ac.ge.2.)then
c         ireac=2
c        else 
c         ireac=3
c        endif
c       else
c        if(ac.ge.2.)then
c         ireac=4
c        else 
c         ireac=5
c        endif
c       endif
c      elseif(iap.eq.11.or.iap.eq.12.or.iap.eq.22)then
c       if(iat.gt.100)then
c        if(ac.ge.2.)then
c         ireac=2
c        else 
c         ireac=3
c        endif
c       elseif(iat.eq.11.or.iat.eq.12.or.iat.eq.22)then
c        ireac=7
c       else
c        ireac=8
c       endif
c      else
c       if(iat.gt.100)then
c        if(ac.ge.2.)then
c         ireac=4
c        else 
c         ireac=5
c        endif
c       elseif(iat.eq.11.or.iat.eq.12.or.iat.eq.22)then
c        ireac=8
c       else
c        ireac=9
c       endif
c      endif
c
c      end
c
c-----------------------------------------------------------------------
      subroutine xEmsI1(iii,kc,omlog)
c-----------------------------------------------------------------------
c plot omlog vs iter
c plot  nr of pomerons vs iter  
c plot number of collisions vs iter  
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'

      parameter(nbin=100)
      common/cmc/ot(0:nbin),zz(0:nbin),i(0:nbin)
     *,yt1,yt2,kx(0:nbin)
      parameter(nbim=100)
      common/cmc1/xp(0:nbim),xt(0:nbim),x(0:nbim),o(0:nbim)
     *,y1,y2,car
      character car*5
      double precision xp,xt,x,omlog,om1intbc
      character ce*8
      double precision plc,s,seedp
      common/cems5/plc,s
      
c      if(iemsi2.eq.0)call utstop('ERROR in XemsI1: iemsi2 = 0&')

       if(iii.eq.1)then

      o(kc)=sngl(omlog)
      nptk=0
      kollx=0
      do ko=1,koll
      nptk=nptk+nprt(ko)
c      if(itpr(ko).gt.0)then
      if(nprt(ko).gt.0)then
       kollx=kollx+1
      endif
      enddo
      zz(kc)=nptk
      kx(kc)=kollx      

        elseif(iii.eq.2)then
      
      call ranfgt(seedp)
      sum=0
      kollx=0
      sumg=0
      kollg=0
      do ko=1,koll
ctp060829       ip=iproj(ko)
ctp060829       it=itarg(ko)
       om1i=sngl(om1intbc(bk(ko)))
ctp060829         wk=1.
ctp060829         wp=0.
ctp060829         wt=0.
       om1g=sngl(om1intbc(bk(ko)))
       sum=sum+om1i
       sumg=sumg+om1g
       if(rangen().lt.1.-exp(-om1i))then
        kollx=kollx+1
       endif
       if(rangen().lt.1.-exp(-om1g))then
        kollg=kollg+1
       endif
      enddo
      call ranfst(seedp)

      x1=0
      x2=nbin
      write(ce,'(f8.2)')sngl(plc)
      
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   log omega       for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name omega-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')       'yrange auto auto '
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis ln[W]"'
      write(ifhi,'(a,a)')  'text 0.5 0.90 "E ='//ce//'"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),o(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i3)')'! nr of coll`s  for event ',nrevt+1
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name coll-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis nr of collisions"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kx(k))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kollx)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),float(kollg)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i3)')'! nr of pom`s  for event ',nrevt+1
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name pom-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis nr of Pomerons"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),zz(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      if(sum.lt.4*zz(nbin))then
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),sum
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbin
      write(ifhi,'(2e11.3)')float(k),sumg
         enddo
      write(ifhi,'(a)')    '  endarray'
      endif
      write(ifhi,'(a)')    'closehisto plot 0'

        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsI2(iii,kc)
c-----------------------------------------------------------------------
c plot quanities vs iter
c   plot 1: <x> for Pomeron vs iter  
c   plot 2: <x> for projectile vs iter  
c   plot 3: <x> for target vs iter 
c arguments: 
c   iii:   modus (1,2)
c   kc:    iteration step
c   omega: config probability
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'

      parameter(nbim=100)
      common/cmc1/xp(0:nbim),xt(0:nbim),x(0:nbim),o(0:nbim)
     *,y1,y2,car
      character car*5
      double precision xp,xt,x,xpo,xpj,xtg
      common/cemsi2/xpo,xpj,xtg

        if(iii.eq.1)then
        
      npom=0
      xpo=0
      do k=1,koll
c       ip=iproj(k)
c       it=itarg(k)
       if(nprmx(k).gt.0)then
        do n=1,nprmx(k)
         if(idpr(n,k).gt.0.and.ivpr(n,k).gt.0)then
          xpo=xpo+xpr(n,k)
          npom=npom+1
         endif
        enddo
       endif
      enddo 
      if(npom.gt.0)xpo=xpo/npom

      npk=0
      xpj=0d0
      do i=1,maproj
       if(xpp(i).lt.0.999)then
        xpj=xpj+xpp(i)!*xmp(i)
        npk=npk+1
       endif
      enddo
      if(npk.gt.0)xpj=xpj/dble(npk)
      
      ntk=0
      xtg=0d0
      do j=1,matarg
       if(xmt(j).lt.0.999)then
        xtg=xtg+xmt(j)!*xpt(j)
        ntk=ntk+1
       endif
      enddo
      if(ntk.gt.0)xtg=xtg/dble(ntk)
     
      x(kc)=xpo
      xp(kc)=xpj
      xt(kc)=xtg

        elseif(iii.eq.2)then
      
      x1=0
      x2=nbim
 
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x  Pom   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxPom-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x Pomeron"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),x(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x proj   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxProj-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x proj"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),xp(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
 
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i3)')    '!   average x targ   for event ',nrevt+1
      write(ifhi,'(a)')       '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name avxTarg-',nrevt+1
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',x1,x2
      write(ifhi,'(a)')    'text 0 0 "xaxis iteration"'
      write(ifhi,'(a)')    'text 0 0 "yaxis average x targ"'
      write(ifhi,'(a)')       'array 2'
         do k=0,nbim
      write(ifhi,'(2e11.3)')float(k),xt(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
        endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsRx(iii,id,xp,xm)
c-----------------------------------------------------------------------
c plot  x+, x-, x, y distribution of remnants
c-----------------------------------------------------------------------

      include 'epos.inc'

      parameter(nbix=50,nbiy=50,nid=2)
      common/cxp/nxp(nid),nxm(nid),nx(nid),ny(nid)
     *,wxp(nbix,nid),wxm(nbix,nid),wx(nbix,nid),wy(nbiy,nid)
     *,xpu,xpo,xmu,xmo,xu,xo,yu,yo,dy
      
      if(iemsrx.eq.0)call utstop('ERROR in XemsRx: iemsrx = 0&')

        if(iii.eq.0)then

      xpu=10/engy**2
      xpo=1
      xmu=10/engy**2
      xmo=1
      xu=10/engy**2
      xo=1
      yu=-alog(engy**2)       
      yo=alog(engy**2)       
      dy=(yo-yu)/nbiy
      do j=1,nid
       nxp(j)=0
       nxm(j)=0
       nx(j)=0
       do i=1,nbix
        wxp(i,j)=0
        wxm(i,j)=0
        wx(i,j)=0
       enddo
       ny(j)=0
       do i=1,nbiy
        wy(i,j)=0
       enddo
      enddo
      
        elseif(iii.eq.1)then
      
      i=0  
      if(xp.lt.xpu)goto1
      i=1+int(alog(xp/xpu)/alog(xpo/xpu)*nbix)
      if(i.gt.nbix)goto1
      if(i.lt.1)goto1
      wxp(i,id)=wxp(i,id)+1
      nxp(id)=nxp(id)+1
1     continue

      if(xm.lt.xmu)goto2
      i=1+int(alog(xm/xmu)/alog(xmo/xmu)*nbix)
      if(i.gt.nbix)goto2
      if(i.lt.1)goto2
      wxm(i,id)=wxm(i,id)+1
      nxm(id)=nxm(id)+1
2     continue

      x=xp*xm
      if(x.lt.xu)goto3
      i=1+int(alog(x/xu)/alog(xo/xu)*nbix)
      if(i.gt.nbix)goto3
      if(i.lt.1)goto3
      wx(i,id)=wx(i,id)+1
      nx(id)=nx(id)+1
3     continue

      if(xm.le.0.)goto4
      if(xp.le.0.)goto4
      y=0.5*alog(xp/xm)
      if(y.lt.yu)goto4
      i=int((y-yu)/dy)+1
      if(i.gt.nbiy)goto4
      if(i.lt.1)goto4
      wy(i,id)=wy(i,id)+1
      ny(id)=ny(id)+1
4     continue
       
        elseif(iii.eq.2)then
      
      do j=1,nid
      if(j.eq.1)iclrem=iclpro
      if(j.eq.2)iclrem=icltar
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant xp distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xpRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xpu,xpo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x+"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x+)"'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      x=xpu*(xpo/xpu)**((i-0.5)/nbix)
      dx=xpu*(xpo/xpu)**(1.*i/nbix)*(1.-(xpo/xpu)**(-1./nbix))
      if(nxp(j).ne.0)write(ifhi,'(2e11.3)')x,wxp(i,j)/dx/nxp(j)
      if(nxp(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')       'openhisto'
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      write(ifhi,'(2e11.3)')x,x**alplea(iclrem)*(1+alplea(iclrem))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
      
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant xm distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xmRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xmu,xmo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x-"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x-)"'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      x=xmu*(xmo/xmu)**((i-0.5)/nbix)
      dx=xmu*(xmo/xmu)**(1.*i/nbix)*(1.-(xmo/xmu)**(-1./nbix))
      if(nxm(j).ne.0)write(ifhi,'(2e11.3)')x,wxm(i,j)/dx/nxm(j)
      if(nxm(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
 
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant x distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name xRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod log ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',xu,xo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant x"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      x=xu*(xo/xu)**((i-0.5)/nbix)
      dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
      if(nx(j).ne.0)write(ifhi,'(2e11.3)')x,wx(i,j)/dx/nx(j)
      if(nx(j).eq.0)write(ifhi,'(2e11.3)')x,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
 
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a)')      '!   remnant y distribution      '
      write(ifhi,'(a)')      '!----------------------------------'
      write(ifhi,'(a,i1)')    'openhisto name yRemnant-',j
      write(ifhi,'(a)')       'htyp lin'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',yu,yo
      write(ifhi,'(a)')    'text 0 0 "xaxis remnant y"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
      write(ifhi,'(a)')       'array 2'
         do i=1,nbix
      y=yu+dy/2.+(i-1)*dy
      if(ny(j).ne.0)write(ifhi,'(2e11.3)')y,wy(i,j)/dy/ny(j)
      if(ny(j).eq.0)write(ifhi,'(2e11.3)')y,0.
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      enddo

        endif
              
      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsPm(iii,ko,nmc)
c-----------------------------------------------------------------------
c m (pomeron number) distribution for different b-bins.
c arguments:
c   iii:  modus (0,1,2)
c   ko:   pair number (1 - AB)
c   nmc:  number of pomerons
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbin=npommx)
      parameter(nbib=32)
      common/cn/wn(0:nbin,nbib),wnmc(0:nbin,nbib),npmx(nbib),nn(nbib)
     &         ,nn2(nbib)
      common/cb1/db,b1,b2,bb(nbib),nbibx
      double precision plc,s,om1intbc
      character ce*8,cb*4
      common/cems5/plc,s
      common/cemspm/sumb(nbib)
      
      if(iemspm.eq.0)call utstop('ERROR in XemsPm: iemspm = 0&')
      
        if(iii.eq.0)then

      do k=1,nbib
       nn(k)=0
       nn2(k)=0
       sumb(k)=0
       do i=0,nbin
        wnmc(i,k)=0
       enddo 
      enddo
      nbibx=6
      b1=0
      b2=2
      db=(b2-b1)/nbibx


        elseif(iii.eq.1)then
        
      k=int((bk(ko)-b1)/db)+1
      if(k.gt.nbibx)k=nbibx
      if(k.lt.1)k=1
      if(nmc.gt.nbin)return
      if(nmc.lt.0)return
      nn(k)=nn(k)+1
      wnmc(nmc,k)=wnmc(nmc,k)+1
      sumb(k)=sumb(k)+bk(ko)


        elseif(iii.eq.2)then
      
      do 1 k=1,nbibx
      
       bb(k)=b1+(k-0.5)*db
       if(maproj.eq.1.and.matarg.eq.1.and.bmaxim.eq.0.)bb(k)=b1
       om1i=sngl(om1intbc(bb(k)))
       do i=0,nbin
        if(i.eq.0)then
         wn(i,k)=exp(-om1i)
        else
         wn(i,k)=wn(i-1,k)*om1i/i
        endif
       if(wn(i,k).gt.0.000001*(1.-exp(-om1i)))npmx(k)=i
       enddo
      
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! distr of Pomeron number vs b'
      write(ifhi,'(a)')   '!##################################'
      write(ce,'(f8.2)')sngl(plc)
      write(cb,'(f4.2)')bb(k)
      if(nn(k).gt.0)then
      write(ifhi,'(a,i1)')    'openhisto name mPom-',k
      write(ifhi,'(a)')       'htyp lru'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(npmx(k))
      write(ifhi,'(a)')    'text 0 0 "xaxis number m of Pomerons"'
      write(ifhi,'(a)')    'text 0 0 "yaxis prob(m)"'
      if(k.eq.1)
     *write(ifhi,'(a,a)')     'text 0.5 0.90 "E ='//ce//'"'
      write(ifhi,'(a,a)')     'text 0.5 0.80 "b ='//cb//'"'
      write(ifhi,'(a)')       'array 2'
         do i=0,nbin
      write(ifhi,'(2e11.3)')float(i),wnmc(i,k)/max(1,nn(k))
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      endif

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! distr of Pomeron number vs b'
      write(ifhi,'(a)')   '!   traditional approach'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1)')    'openhisto name mPomTradi-',k
      write(ifhi,'(a)')       'htyp lba'
      write(ifhi,'(a)')       'xmod lin ymod log'
      write(ifhi,'(a,2e11.3)')'xrange',0.,float(npmx(k))
      write(ifhi,'(a)')       'array 2'
         do i=0,nbin
      write(ifhi,'(2e11.3)')float(i),wn(i,k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

 1    continue


      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsB(iii,jjj,ko)
c-----------------------------------------------------------------------
c b distribution at different stages
c arguments:
c   iii:  modus (0,1,2)
c   jjj:  stage or type of interaction 
c     just after Metropolis:
c           1 ... all 
c           2 ... interaction 
c     after defining diffraction:
c           3 ... nothing
c           4 ... cut
c           5 ... diffr
c           6 ... cut + diffr
c   ko:   pair number (1 - AB)
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      include 'epos.incsem'
      parameter(njjj=6)
      parameter(nbib=32)
      common/cxemsb1/w(0:njjj,nbib),nn(njjj)
      common/cxemsb2/db,b1,b2
      common/cxemsb3/njjj1
      double precision PhiExact,om1intbi,PhiExpo!,PhiUnit
      common/geom/rmproj,rmtarg,bmax,bkmx
      dimension uua2(nbib),uuo2(nbib),uu3(nbib)
      
      if(iemsb.eq.0)call utstop('ERROR in XemsB: iemsB = 0&')
      
        if(iii.eq.0)then

      do k=1,nbib
       do j=0,njjj
        w(j,k)=0
       enddo 
      enddo
      do j=1,njjj
       nn(j)=0
      enddo 
      njjj1=0

        elseif(iii.eq.1)then
        
      b1=0
      b2=bkmx*1.2
      db=(b2-b1)/nbib
      k=int((bk(ko)-b1)/db)+1
      if(k.gt.nbib)return
      if(k.lt.1)return
      w(jjj,k)=w(jjj,k)+1
      nn(jjj)=nn(jjj)+1
      if(jjj.eq.1)njjj1=1

        elseif(iii.eq.2)then
        
      if(njjj1.ne.1)call utstop
     &('xEmsB must be called also with jjj=1&')  
      ymax=0
      kollini=koll
      koll=1
      do k=1,nbib 
       x=b1+(k-0.5)*db  
       y=w(1,k)/nn(1)/(pi*((x+0.5*db)**2-(x-0.5*db)**2)) 
       ymax=max(ymax,y)
      enddo
      fk=bkmx**2*pi
      ymax=1.4
      
      do 1 j=1,njjj
       if(nn(j).eq.0)goto1

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr exact theory '
      write(ifhi,'(a)')   '!##################################'
         if(j.ge.2.and.j.le.6)then
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Exact'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db 
      if(j.eq.2)then
        uuo2(k)=sngl(PhiExpo(1.,1.d0,1.d0,engy**2,b))
        uua2(k)=min(uuo2(k),max(0.,
     &          sngl(PhiExact(1.,1.d0,1.d0,engy**2,b))))
        uu3(k)=sngl(min(50d0,exp(om1intbi(b,2)/dble(r2hads(iclpro)
     &                                             +r2hads(icltar)))))
      endif
      if(j.eq.2)y=(1.-uua2(k))
      if(j.eq.3)y=uua2(k)
      if(j.eq.4)y=(1.-uua2(k)*uu3(k))
      if(j.eq.5)y=uua2(k)*(uu3(k)-1.)
      if(j.eq.6)y=(1.-uua2(k))
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
         endif
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr unitarized theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Unit'
      write(ifhi,'(a)')       'htyp lbf xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db   
      if(j.eq.1)y=1
      if(j.eq.2)y=(1.-uuo2(k))
      if(j.eq.3)y=uuo2(k)
      if(j.eq.4)y=(1.-uuo2(k)*uu3(k))
      if(j.eq.5)y=uuo2(k)*(uu3(k)-1.)
      if(j.eq.6)y=(1.-uuo2(k))
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr for cross section '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Unit'
      write(ifhi,'(a)')       'htyp lge xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1+(k-0.5)*db   
      if(j.eq.1)y=1
      if(j.eq.2)y=(1.-(uuo2(k)+uua2(k))*0.5)
      if(j.eq.3)y=(uuo2(k)+uua2(k))*0.5
      if(j.eq.4)y=(1.-(uuo2(k)+uua2(k))*0.5*uu3(k))
      if(j.eq.5)y=(uuo2(k)+uua2(k))*0.5*(uu3(k)-1.)
      if(j.eq.6)y=(1.-(uuo2(k)+uua2(k))*0.5)
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name b',j,'Simu'
      write(ifhi,'(a)')       'htyp lrf xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.0,b2
      write(ifhi,'(a,2e11.3)')'yrange',0.,ymax
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      if(j.eq.1)write(ifhi,'(a)')'text 0.1 0.35 "after Metropolis"'
      if(j.eq.1)write(ifhi,'(a)')'text 0.2 0.20 "all "'
      if(j.eq.2)write(ifhi,'(a)')'text 0.3 0.85 "after Metropolis"'
      if(j.eq.2)write(ifhi,'(a)')'text 0.5 0.70 "interaction "'
      if(j.eq.3)write(ifhi,'(a)')'text 0.3 0.85 "nothing"'
      if(j.eq.4)write(ifhi,'(a)')'text 0.3 0.85 "cut"'
      if(j.eq.5)write(ifhi,'(a)')'text 0.3 0.85 "diffr"'
      if(j.eq.6)write(ifhi,'(a)')'text 0.3 0.85 "cut + diffr"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      x=b1+(k-0.5)*db  
      if(j.eq.1)y=fk*w(j,k)/nn(1)/(pi*((x+0.5*db)**2-(x-0.5*db)**2)) 
      if(j.ne.1)y=0.
      if(j.ne.1.and.w(1,k).ne.0.)y=w(j,k)/w(1,k)
      if(nn(j).gt.0)write(ifhi,'(2e11.3)')x,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'
        
   1  continue

      koll=kollini

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsBg(iii,jjj,ko)
c-----------------------------------------------------------------------
c b distribution at different stages for different group
c arguments:
c   iii:  modus (0,1,2,3)
c   jjj:  group of interaction (1,2 ... ,7)
c   ko:   pair number (1 - AB)
c-----------------------------------------------------------------------
      include 'epos.inc'
      include 'epos.incems'
      parameter(njjj=7)
      parameter(nbib=16)
      common/cxemsb4/wg(-1:njjj,nbib),nng(nbib),uug(nbib),kollx
      common/cxemsb5/dbg,b1g,b2g
      common/cxemsb6/njjj0
      double precision seedp,PhiExpo!,PhiExact
      common/geom/rmproj,rmtarg,bmax,bkmx
      
      if(iemsbg.eq.0)call utstop('ERROR in XemsBg: iemsbg = 0&')
      
        if(iii.eq.0)then

      do k=1,nbib
       nng(k)=0
       do j=-1,njjj
        wg(j,k)=0
       enddo 
      enddo
      njjj0=0
      kollx=0

        elseif(iii.eq.1)then
        
      b1g=0
      b2g=bkmx*1.2
      dbg=(b2g-b1g)/nbib
      k=int((bk(ko)-b1g)/dbg)+1
      if(k.gt.nbib)return
      if(k.lt.1)return
      if(jjj.eq.-1.or.jjj.eq.0)then
        wg(jjj,k)=wg(jjj,k)+1
      else
        wg(jjj,k)=wg(jjj,k)+1
        nng(k)=nng(k)+1
      endif
      if(jjj.eq.0)njjj0=1

        elseif(iii.eq.3)then

          call ranfgt(seedp)
          do k=1,koll
            om1i=sngl(om1intc(k))
            if(rangen().lt.1.-exp(-om1i))then
c            om1i=sngl(PhiExpo(1.,1.d0,1.d0,engy*engy,bk(k)))
c            if(rangen().lt.1.-om1i)then
              kollx=kollx+1
            endif
          enddo
          call ranfst(seedp)
        
        elseif(iii.eq.2)then
        
      if(njjj0.ne.1)call utstop
     &('xEmsBg must be called also with jjj=0&')  
      ymax=1.4
      kollini=koll
      koll=1

      wtot=1.
      if(matarg+maproj.gt.2)then
      wtot=0.
      do k=1,nbib
       wtot=wtot+wg(-1,k)
      enddo
      wtot=wtot/float(kollx)
      endif

      do 1 j=1,njjj

      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distribution simulation'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name bg',j,'Simu'
      write(ifhi,'(a)')       'htyp lin xmod lin ymod lin'
      write(ifhi,'(a,2e11.3)')'xrange',0.,b2g
      write(ifhi,'(a,2e11.3)')'yrange',0.,ymax
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      if(wtot.gt.0.d0)
     &write(ifhi,'(a,f7.4,a)')    'text 0.5 0.8 "alpha=',1./wtot,'"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg
      y=0.
      if(nng(k).ne.0.and.wg(0,k).ne.0)
     &              y=wg(j,k)/float(nng(k))*wg(-1,k)/wg(0,k)!/wtot
c      if(wg(0,k).ne.0..and.nng(k).ne.0)y=wg(j,k)/nng(k)*wg(-1,k)/wg(0,k) 
c!???????????? better normalization ? probability to have an interaction 
c in epos compared to eikonal probability, instead of normalized by the 
c probability of a collision for a pair (the number collision/number 
c active pair).  
      uug(k)=uug(k)+y
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
   1  continue
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr tot simul theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')  'openhisto name btotSimu'
      write(ifhi,'(a)')       'htyp pfc xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg  
      write(ifhi,'(2e11.3)')b,uug(k)
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0-'
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a)')   '! b distr unitarized theory '
      write(ifhi,'(a)')   '!##################################'
      write(ifhi,'(a,i1,a)')  'openhisto name bg',j,'Unit'
      write(ifhi,'(a)')       'htyp lba xmod lin ymod lin'
      write(ifhi,'(a)')    'text 0 0 "xaxis impact parameter b"'
      write(ifhi,'(a)')    'text 0 0 "yaxis P(b)"'
      write(ifhi,'(a)')       'array 2'
         do k=1,nbib
      b=b1g+(k-0.5)*dbg   
c      a1=PhiExact(1.,1.d0,1.d0,engy**2,b)
       a1=sngl(PhiExpo(1.,1.d0,1.d0,engy**2,b))
      y=(1.-a1)
      write(ifhi,'(2e11.3)')b,y
         enddo
      write(ifhi,'(a)')    '  endarray'
      write(ifhi,'(a)')    'closehisto plot 0'

      koll=kollini

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsPx(iii,xmc,ymc,npos)
c-----------------------------------------------------------------------
c plot  x-distribution and y-distribution of Pomerons 
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incems'
      common/geom/rmproj,rmtarg,bmax,bkmx
      
      parameter(nbix=30,nbib=51)
      common/cx/x(2,nbix),dx(2,nbix),wxmc(2,nbix),wxmcI(2,nbix)
     * ,xl(2,nbix),dxl(2,nbix),wxp(2,nbix),wxm(2,nbix),wxpI(2,nbix)
     *,wxmI(2,nbix),wxpY(2,nbix),wxmY(2,nbix),wxmcY(2,nbix)
      parameter(nbiy=50)
      common/cy/y(nbiy),wymc(nbiy),wymcY(nbiy),wymcI(nbiy),nyp,nym
      double precision PomIncXExact,PomIncPExact,PomIncMExact,dcel
      double precision PomIncXIExact,PomIncPIExact,PomIncMIExact
      common/ems3/dcel,ad
      common/cemspx/xu,xo,yu,yo,dy,xlu,xlo,bb,nn,db,mm,nm,nt
      character mod*5, imod*5, txtxm*6
      
      nposi=5
      
      if(iemspx.eq.0)call utstop('ERROR in XemsPx: iemspx = 0&')
      
      if(iii.eq.0)then
      
       xu=0.1/engy**2
       xo=1.
       xlu=0.01/engy
       xlo=1.
       yu=-alog(engy**2)       
       yo=alog(engy**2)       
       dy=(yo-yu)/nbiy
        do i=1,nbix
        x(1,i)=xu*(xo/xu)**((i-0.5)/nbix)
        x(2,i)=xu+(xo-xu)*((i-0.5)/nbix)
        dx(1,i)=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
        dx(2,i)=(xo-xu)/nbix
        wxmc(1,i)=0.
        wxmc(2,i)=0.
        wxmcI(1,i)=0.
        wxmcI(2,i)=0.
        wxmcY(1,i)=0.
        wxmcY(2,i)=0.
       enddo
       do i=1,nbix
        xl(1,i)=xlu*(xlo/xlu)**((i-0.5)/nbix)
        xl(2,i)=xlu+(xlo-xlu)*((i-0.5)/nbix)
        dxl(1,i)=xlu*(xlo/xlu)**(1.*i/nbix)*(1.-(xlo/xlu)**(-1./nbix))
        dxl(2,i)=(xlo-xlu)/nbix
        wxp(1,i)=0.
        wxp(2,i)=0.
        wxm(1,i)=0.
        wxm(2,i)=0.
        wxpI(1,i)=0.
        wxpI(2,i)=0.
        wxmI(1,i)=0.
        wxmI(2,i)=0.
        wxpY(1,i)=0.
        wxpY(2,i)=0.
        wxmY(1,i)=0.
        wxmY(2,i)=0.
       enddo
       do i=1,nbiy
        y(i)=yu+dy/2.+float(i-1)*dy
        wymc(i)=0.
        wymcI(i)=0.
        wymcY(i)=0.
       enddo
       mm=0
       nt=0
       nyp=0
       nym=0
       db=bkmx*2./float(nbib-1)
       
      elseif(iii.eq.1)then
      
       xp=sqrt(xmc)*exp(ymc)
       xm=sqrt(xmc)*exp(-ymc)
       mm=mm+1
       
       if(xmc.lt.xu)goto11
       i=1+int(alog(xmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxmc(1,i)=wxmc(1,i)+1.
       if(npos.eq.1)    wxmcI(1,i)=wxmcI(1,i)+1.
       if(npos.eq.nposi)wxmcY(1,i)=wxmcY(1,i)+1.
1      continue
       i=1+int((xmc-xu)/(xo-xu)*nbix)
       if(i.gt.nbix)goto11
       if(i.lt.1)goto11
       wxmc(2,i)=wxmc(2,i)+1.
       if(npos.eq.1)    wxmcI(2,i)=wxmcI(2,i)+1.
       if(npos.eq.nposi)wxmcY(2,i)=wxmcY(2,i)+1.
11     continue

       if(xp.lt.xlu)goto12
       i=1+int(alog(xp/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto2
       if(i.lt.1)goto2
       wxp(1,i)=wxp(1,i)+1.
       if(npos.eq.1)    wxpI(1,i)=wxpI(1,i)+1.
       if(npos.eq.nposi)wxpY(1,i)=wxpY(1,i)+1.
2      continue
       i=1+int((xp-xlu)/(xlo-xlu)*nbix)
       if(i.gt.nbix)goto12
       if(i.lt.1)goto12
       wxp(2,i)=wxp(2,i)+1.
       if(npos.eq.1)    wxpI(2,i)=wxpI(2,i)+1.
       if(npos.eq.nposi)wxpY(2,i)=wxpY(2,i)+1.
12     continue

       if(xm.lt.xlu)goto13
       i=1+int(alog(xm/xlu)/alog(xlo/xlu)*nbix)
       if(i.gt.nbix)goto3
       if(i.lt.1)goto3
       wxm(1,i)=wxm(1,i)+1.
       if(npos.eq.1)    wxmI(1,i)=wxmI(1,i)+1.
       if(npos.eq.nposi)wxmY(1,i)=wxmY(1,i)+1.
3      continue
       i=1+int((xm-xlu)/(xlo-xlu)*nbix)
       if(i.gt.nbix)goto13
       if(i.lt.1)goto13
       wxm(2,i)=wxm(2,i)+1.
       if(npos.eq.1)    wxmI(2,i)=wxmI(2,i)+1.
       if(npos.eq.nposi)wxmY(2,i)=wxmY(2,i)+1.
13     continue

       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)return
       if(i.lt.1)return
       wymc(i)=wymc(i)+1
       if(npos.eq.1)    wymcI(i)=wymcI(i)+1
       if(npos.eq.nposi)wymcY(i)=wymcY(i)+1
       if(ymc.gt.0)nyp=nyp+1
       if(ymc.lt.0)nym=nym+1
       
      elseif(iii.eq.2)then
       
       if(maproj.eq.1.and.matarg.eq.1.and.bminim.eq.bmaxim)then
        mmmm=1  
        bb=bmaxim
        ff=float(nrevt)/float(ntevt)
        imod='   dn'
       elseif(maproj.eq.1.and.matarg.eq.1)then
        mmmm=3    
        ff=1.
        imod='   dn'
       elseif(bminim.lt.0.001.and.bmaxim.gt.20)then
        mmmm=2   
        area=pi*(rmproj+rmtarg)**2
        ff=area*float(nrevt)/float(ntevt)/(maproj*matarg)/sigine*10 
        imod='   dn'
       else
        write(ifmt,*)'xEmsPx ignored' 
        return 
       endif 
       
       kk1=nint(xpar1)
       kk2=nint(xpar2)

       do kk=kk1,kk2
      
       if(kk.eq.1)mod=' log '
       if(kk.eq.2)mod=' lin '
      
       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Pomeron x distribution    '//mod
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a)')  'openhisto name xPomSimuL'//mod(3:4)
       write(ifhi,'(a)')  'htyp lru xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')    'text 0 0 "xaxis x?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx?PE!"'
       if(kk.eq.1)write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "f=',ff,'"'
       if(kk.eq.2)write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "f=',ff,'"'
       write(ifhi,'(a)')       'array 2'
       s1=0
       do i=1,nbix
       u=x(kk,i)
       z=ff*wxmc(kk,i)/dx(kk,i)/nrevt
       s1=s1+z*dx(kk,i)
        write(ifhi,'(2e11.3)')u,z  
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0-'
              
       write(ifhi,'(a)')       'openhisto name xPomUnitL'//mod(3:4)
       write(ifhi,'(a)')  'htyp lba xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xu,xo
       write(ifhi,'(a)')    'text 0 0 "xaxis x?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx?PE!"'
       write(ifhi,'(a)')       'array 2'
       s2=0
       do i=1,nbix
        u=x(kk,i)
        if(mmmm.eq.1)z=PomIncXExact(dble(u),bb)   
        if(mmmm.eq.2)z=PomIncXIExact(dble(u))/sigine*10     
        if(mmmm.eq.3)z=PomIncXIExact(dble(u))/sigine*10     
        s2=s2+dx(kk,i)*z
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f5.3,a,f5.3,a)')
     *                       'text .1 .85 "I= ',s1,' (',s2,')"'
       write(ifhi,'(a)')    'closehisto plot 0'
       
       write(ifhi,'(a)')           '!--------------------------------'
       write(ifhi,'(a)')           '!   Pomeron y distribution   '//mod
       write(ifhi,'(a)')           '!--------------------------------'
       
       write(ifhi,'(a)')       'openhisto name yPomSimuL'//mod(3:4)
       write(ifhi,'(a)')       'htyp lru xmod lin ymod'//mod
       write(ifhi,'(a,2e11.3)')'xrange',yu,yo
       write(ifhi,'(a)')    'text 0 0 "xaxis y?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom!/dy?PE!"'
       write(ifhi,'(a,f5.2,a)')'text 0.1 0.7 "f=',ff,'"'
       write(ifhi,'(a)')       'array 2'
       s1=0
       do i=1,nbiy
       u=y(i)
       z=ff*wymc(i)/dy/nrevt
       s1=s1+z*dy
        write(ifhi,'(2e11.3)')u,z       
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0'
       
       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   Pomeron x+ distribution    '//mod
       write(ifhi,'(a)')       '!----------------------------------'
       
       write(ifhi,'(a)')   'openhisto name xpPomSimuL'//mod(3:4)
       write(ifhi,'(a)')   'htyp lru xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       if(kk.eq.1)write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "f=',ff,'"'
       if(kk.eq.2)write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "f=',ff,'"'
       write(ifhi,'(a)')       'array 2'
       s1=0
       do i=1,nbix
       u=xl(kk,i)
       z=ff*wxp(kk,i)/dxl(kk,i)/nrevt 
       s1=s1+z*dxl(kk,i)
        write(ifhi,'(2e11.3)')u,z       
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0-'
             
       write(ifhi,'(a)')       'openhisto name xpPomUnitL'//mod(3:4)
       write(ifhi,'(a)')   'htyp lba xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x+?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx+?PE!"'
       write(ifhi,'(a)')       'array 2'
       s2=0
       do i=1,nbix
        u=xl(kk,i)
        if(mmmm.eq.1)z=PomIncPExact(dble(u),bb) 
        if(mmmm.eq.2)z=PomIncPIExact(dble(u))/sigine*10    
        if(mmmm.eq.3)z=PomIncPIExact(dble(u))/sigine*10   
        s2=s2+dxl(kk,i)*z
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f5.3,a,f5.3,a)')
     *                       'text .1 .85 "I= ',s1,' (',s2,')"'
       write(ifhi,'(a)')    'closehisto plot 0'

       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a)')       '!   x-?PE! distribution    '//mod
       write(ifhi,'(a)')       '!----------------------------------'
       
       write(ifhi,'(a)')   'openhisto name xmPomSimuL'//mod(3:4)
       write(ifhi,'(a)')   'htyp lru xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x-?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx-?PE!"'
       if(kk.eq.1)write(ifhi,'(a,f5.2,a)')'text 0.1 0.3 "f=',ff,'"'
       if(kk.eq.2)write(ifhi,'(a,f5.2,a)')'text 0.1 0.1 "f=',ff,'"'
       write(ifhi,'(a)')       'array 2'
       s1=0
       do i=1,nbix
       u=xl(kk,i)
       z=ff*wxm(kk,i)/dxl(kk,i)/nrevt
       s1=s1+z*dxl(kk,i)
        write(ifhi,'(2e11.3)')u,z       
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a)')    'closehisto plot 0-'
              
       write(ifhi,'(a)')       'openhisto name xmPomUnitL'//mod(3:4)
       write(ifhi,'(a)')   'htyp lba xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange',xlu,xlo
       write(ifhi,'(a)')    'text 0 0 "xaxis x-?PE!"'
       write(ifhi,'(a)') 'text 0 0 "yaxis'//imod//'?Pom! / dx-"'
       write(ifhi,'(a)')       'array 2'
       s2=0
       do i=1,nbix
        u=xl(kk,i)
        if(mmmm.eq.1)z=PomIncMExact(dble(u),bb) 
        if(mmmm.eq.2)z=PomIncMIExact(dble(u))/sigine*10   
        if(mmmm.eq.3)z=PomIncMIExact(dble(u))/sigine*10 
        s2=s2+dxl(kk,i)*z
        write(ifhi,'(2e11.3)')u,z
       enddo
       write(ifhi,'(a)')    '  endarray'
       write(ifhi,'(a,f5.3,a,f5.3,a)')
     *                       'text .1 .85 "I= ',s1,' (',s2,')"'
       write(ifhi,'(a)')    'closehisto plot 0'

  !................................................................
  
       xm=-1. !xm integration
       txtxm='xm int'
       do jjb=0,3
       b=jjb*0.5
       do jj=0,2
       
       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a,3i1)')   '!   ffom11    '//mod,jjb,jj
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,2i1)')'openhisto name ffom11L'//mod(3:4),jjb,jj+8
       write(ifhi,'(a)')    'htyp lin xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange ',xlu,xlo
       write(ifhi,'(a)')'txt "xaxis  x+?PE!"'
       write(ifhi,'(a)')'txt "yaxis dn?Pom! / dx+?PE! "'
       write(ifhi,'(a)')'text 0.05 0.1  "fit and exact, all contrib."'
       if(jjb.lt.3)write(ifhi,'(a,f4.1,3a)')
     *             'txt "title ffom11   b =',b,'   ',txtxm,'"'
       if(jjb.ge.3)write(ifhi,'(3a)')
     *             'txt "title ffom11   b aver   ',txtxm,'"'
       write(ifhi,'(a)')       'array 2'
       do i=1,nbix
       u=xl(kk,i)
       if(jjb.lt.3.and.jj.eq.0)z= ffom11(u,xm,b,-1,-1) 
       if(jjb.lt.3.and.jj.eq.1)z= ffom11(u,xm,b,0,5) 
       if(jjb.lt.3.and.jj.eq.2)z= ffom11(u,xm,b,0,4) 
       if(jjb.eq.3.and.jj.eq.0)z=ffom11a(u,xm,-1,-1) 
       if(jjb.eq.3.and.jj.eq.1)z=ffom11a(u,xm,0,5) 
       if(jjb.eq.3.and.jj.eq.2)z=ffom11a(u,xm,0,4) 
        write(ifhi,'(2e11.3)')u,z  
       enddo
       write(ifhi,'(a)')    '  endarray'
       if(jj.le.1)write(ifhi,'(a)')    'closehisto plot 0-'
       if(jj.eq.2)write(ifhi,'(a)')    'closehisto plot 0'

       enddo
       enddo

       do jjb=0,3
       b=jjb*0.5
       do jjj=1,6
       jj=jjj
       if(jjj.eq.6)jj=0
       
       write(ifhi,'(a)')       '!----------------------------------'
       write(ifhi,'(a,3i1)')   '!   ffom11    '//mod,jjb,jj
       write(ifhi,'(a)')       '!----------------------------------'

       write(ifhi,'(a,3i1)')'openhisto name om1ffL'//mod(3:4),jjb,jj
       if(jj.ne.0)write(ifhi,'(a)')    'htyp lin xmod'//mod//'ymod log'
       if(jj.eq.0)write(ifhi,'(a)')    'htyp lro xmod'//mod//'ymod log'
       write(ifhi,'(a,2e11.3)')'xrange ',xlu,xlo
       if(jj.eq.1)then
       write(ifhi,'(a)') 'txt "xaxis  x+?PE!"'
       write(ifhi,'(a)') 'txt "yaxis  dn?Pom! / dx+?PE!  "'
       if(kk.eq.2)then
        write(ifhi,'(a)') 'text 0.1 0.2  "soft sea-sea"'
        write(ifhi,'(a)') 'text 0.1 0.1  "val-sea sea-val val-val"'
       else
        write(ifhi,'(a)') 'text 0.05 0.8  "soft"'
        write(ifhi,'(a)') 'text 0.05 0.7  "diff"'
        write(ifhi,'(a)') 'text 0.05 0.6  "sea-sea"'
        write(ifhi,'(a)') 'text 0.05 0.5  "val-sea"'
        write(ifhi,'(a)') 'text 0.05 0.4  "sea-val"'
        write(ifhi,'(a)') 'text 0.05 0.3  "val-val"'
      endif
       if(jjb.lt.3)write(ifhi,'(a,f4.1,3a)')
     *             'txt "title ffom11   b =',b,'  ',txtxm,'"'
       if(jjb.ge.3)write(ifhi,'(3a)')
     *             'txt "title ffom11   b aver  ',txtxm,'"'
       endif
       write(ifhi,'(a)')       'array 2'
       do i=1,nbix
       u=xl(kk,i)
       if(jjb.lt.3)z= ffom11(u,xm,b,jj,jj)
       if(jjb.eq.3)z=ffom11a(u,xm,jj,jj)
       write(ifhi,'(2e11.3)')u,z  
       enddo
       write(ifhi,'(a)')    '  endarray'
       if(jjj.ne.6)write(ifhi,'(a)')    'closehisto plot 0-'
       if(jjj.eq.6)write(ifhi,'(a)')    'closehisto plot 0'

       enddo
       enddo
       
      enddo

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsP2(iii,jaa,jex,xpd,xmd,xpb,xmb,pt1,pt2)
c-----------------------------------------------------------------------
c plot  x+ distributions of Pomeron ends (PE) (xpd)
c          and Pomeron's in Born (IB) partons (xpb), 
c     and pt dist of Pomeron's out Born (OB) partons 
c       integrated over x- bins (xmd,xmb)
c  iii=0: initialize
c  ii=1: fill arrays
c  iii>=2: make histogram 
c           (2 - Pomeron end PE, 3 - in Born IB, 4 - out Born OB)
c  jaa: type of semihard Pomeron 
c         1= sea-sea, 2= val=sea, 3= sea-val, 4= val-val
c         5= all  for iii=2
c  jex: emission type 
c         1= no emission, 2= proj emis, 3= targ emis, 4= both sides
c         5= all  for iii=2
c-----------------------------------------------------------------------

      include 'epos.inc'
      include 'epos.incsem'
      include 'epos.incems'
      common/geom/rmproj,rmtarg,bmax,bkmx
      parameter(nbixp=25,nbixm=5,nbipt=20)
      common/cxb/xlp(2,nbixp),dxlp(2,nbixp)
     *          ,xlm(2,nbixm),dxlm(2,nbixm) 
     *          ,wxb(2,4,4,nbixp,nbixm)
     *          ,wxe(2,4,4,nbixp,nbixm)
      common/cptb/ptu,pto,ptob(nbipt),wptob(4,4,nbipt)
      common/cemspbx/xlub1,xlub2,xlob
ctp060829      character imod*5

      if(iemspbx.eq.0)call utstop('ERROR in xEmsP2: iemspbx = 0&')
      
      if(iii.eq.0)then
      
       xlub1=0.01/engy
       xlub2=0.
       xlob=1.
       do i=1,nbixp
        xlp(1,i)=xlub1*(xlob/xlub1)**((i-0.5)/nbixp)
        xlp(2,i)=xlub2+(xlob-xlub2)*((i-0.5)/nbixp)
        dxlp(1,i)=xlub1*(xlob/xlub1)**(1.*i/nbixp)
     *             *(1.-(xlob/xlub1)**(-1./nbixp))
        dxlp(2,i)=(xlob-xlub2)/nbixp
       enddo        
       do i=1,nbixm
        xlm(1,i)=xlub1*(xlob/xlub1)**((i-0.5)/nbixm)
        xlm(2,i)=xlub2+(xlob-xlub2)*((i-0.5)/nbixm)
        dxlm(1,i)=xlub1*(xlob/xlub1)**(1.*i/nbixm)
     *             *(1.-(xlob/xlub1)**(-1./nbixm))
        dxlm(2,i)=(xlob-xlub2)/nbixm
       enddo        
       do i=1,nbixp
       do j=1,nbixm
       do jaai=1,4
       do jexi=1,4
        wxb(1,jaai,jexi,i,j)=0.
        wxb(2,jaai,jexi,i,j)=0.
        wxe(1,jaai,jexi,i,j)=0.
        wxe(2,jaai,jexi,i,j)=0.
       enddo
       enddo
       enddo
       enddo
       ptu=2
       pto=20
       do i=1,nbipt
       ptob(i)=ptu+(pto-ptu)*(i-0.5)/nbipt
       do jaai=1,4
       do jexi=1,4
       wptob(jaai,jexi,i)=0
       enddo
       enddo
       enddo
       
      elseif(iii.eq.1)then
      
       xp=xpb
       xm=xmb
       if(xp.lt.xlub1)goto2
       if(xm.lt.xlub1)goto2
       i=1+int(alog(xp/xlub1)/alog(xlob/xlub1)*nbixp)
       if(i.gt.nbixp)goto2
       if(i.lt.1)goto2
       j=1+int(alog(xm/xlub1)/alog(xlob/xlub1)*nbixm)
       if(j.gt.nbixm)goto2
       if(j.lt.1)goto2
       wxb(1,jaa,jex,i,j)=wxb(1,jaa,jex,i,j)+1.
2      continue

       if(xp.lt.xlub2)goto12
       if(xm.lt.xlub2)goto12
       i=1+int((xp-xlub2)/(xlob-xlub2)*nbixp)
       if(i.gt.nbixp)goto12
       if(i.lt.1)goto12
       j=1+int((xm-xlub2)/(xlob-xlub2)*nbixm)
       if(j.gt.nbixm)goto12
       if(j.lt.1)goto12
       wxb(2,jaa,jex,i,j)=wxb(2,jaa,jex,i,j)+1.
12     continue

       xp=xpd
       xm=xmd
       if(xp.lt.xlub1)goto22
       if(xm.lt.xlub1)goto22
       i=1+int(alog(xp/xlub1)/alog(xlob/xlub1)*nbixp)
       if(i.gt.nbixp)goto22
       if(i.lt.1)goto22
       j=1+int(alog(xm/xlub1)/alog(xlob/xlub1)*nbixm)
       if(j.gt.nbixm)goto22
       if(j.lt.1)goto22
       wxe(1,jaa,jex,i,j)=wxe(1,jaa,jex,i,j)+1.
  22   continue

       if(xp.lt.xlub2)goto32
       if(xm.lt.xlub2)goto32
       i=1+int((xp-xlub2)/(xlob-xlub2)*nbixp)
       if(i.gt.nbixp)goto32
       if(i.lt.1)goto32
       j=1+int((xm-xlub2)/(xlob-xlub2)*nbixm)
       if(j.gt.nbixm)goto32
       if(j.lt.1)goto32
       wxe(2,jaa,jex,i,j)=wxe(2,jaa,jex,i,j)+1.
  32   continue

       do m=1,2
       if(m.eq.1)pt=pt1
       if(m.eq.2)pt=pt2
       i=1+int((pt-ptu)/(pto-ptu)*nbipt)
       if(i.lt.1)goto42
       if(i.gt.nbipt)goto42
       wptob(jaa,jex,i)=wptob(jaa,jex,i)+1
   42  continue   
       enddo
         
      elseif(iii.ge.2)then

       if(maproj.eq.1.and.matarg.eq.1.and.bminim.eq.bmaxim)then
ctp060829        mmmm=1  
ctp060829        bb=bmaxim
        ff=float(nrevt)/float(ntevt)
ctp060829        imod='   dn'
       elseif(maproj.eq.1.and.matarg.eq.1)then
ctp060829        mmmm=3    
        ff=1.
ctp060829        imod='   dn'
       elseif(bminim.lt.0.001.and.bmaxim.gt.20)then
ctp060829        mmmm=2   
        area=pi*(rmproj+rmtarg)**2
        ff=area*float(nrevt)/float(ntevt)/(maproj*matarg)/sigine*10 
ctp060829        imod='   dn'
       else
        write(ifmt,*)'xEmsP2 ignored' 
        return 
       endif 
       
       j1=1  !nint(xpar1)   !first xminus bin
       j2=5  !nint(xpar2)   !last xminus bin
       if(iii.eq.4)j2=1
       kkk=2 !nint(xpar3)   !1 (log binning) 2 (lin binning)
       if(kkk.eq.1)then
ctp060829         xmi1=xlub1*(xlob/xlub1)**((j1-1.)/nbixm)
ctp060829         xmi2=xlub1*(xlob/xlub1)**((j2-0.)/nbixm)
         xlub=xlub1
       elseif(kkk.eq.2)then
ctp060829         xmi1=xlub2+(xlob-xlub2)*((j1-1.)/nbixm)
ctp060829         xmi2=xlub2+(xlob-xlub2)*((j2-0.)/nbixm)
         xlub=xlub2
       endif

       jaa1=jaa
       jaa2=jaa
       jex1=jex
       jex2=jex
       if(jaa.eq.5)then
       jaa1=1
       jaa2=4
       endif
       if(jex.eq.5)then
       jex1=1
       jex2=4
       endif

       if(jex.eq.1)then
        je1=0
        je2=0   
       elseif(jex.eq.2)then
        je1=1
        je2=0     
       elseif(jex.eq.3)then
        je1=0
        je2=1     
       elseif(jex.eq.4)then
        je1=1
        je2=1     
       elseif(jex.eq.5)then
        je1=2
        je2=2   
       endif   

       if(iii.eq.2)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   PE    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

        sum=ffom12aii(jaa,je1,je2)
        write(ifhi,'(a,2i1)')'openhisto name ffom12a',jaa,jex
        write(ifhi,'(a)')'htyp lin xmod lin ymod log'
        write(ifhi,'(a,2e11.3)')'xrange ',xlub,xlob
        write(ifhi,'(a)')    'txt "xaxis  x+?PE!"'
        write(ifhi,'(a)')    'txt "yaxis dn?semi! / dx+?PE!    "'
       write(ifhi,'(a,2i1,a)')'txt "title ffom12a + MC   (',jaa,jex,')"'
        write(ifhi,'(a)')    'array 2'
        do i=1,nbixp
         u=xlp(kkk,i)
         z=ffom12ai(u,jaa1,jaa2,je1,je2)
         write(ifhi,'(2e11.3)')u,z  
        enddo
        write(ifhi,'(a)')    '  endarray'
        if(jex.eq.5)then
          write(ifhi,'(a)')    'closehisto plot 0-'
          write(ifhi,'(a,2i1)')'openhisto name ffom11',jaa,jex
          write(ifhi,'(a)')'htyp lba'
          write(ifhi,'(a)')'text 0.05 0.5 "+ ffom11a "'
          write(ifhi,'(a)')'array 2'
          do i=1,nbixp
           u=xlp(kkk,i)
           z=ffom11a(u,-1.,jaa1,jaa2)
           write(ifhi,'(2e11.3)')u,z  
          enddo
          write(ifhi,'(a)')    '  endarray'
        endif

       elseif(iii.eq.3)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   IB    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

    !.......total integral
        s2min=4*q2min
        zmin=s2min/engy**2
        zmax=1
        xpmin0 = 0.01/engy
        xpmax=1
        ig1=3
        ig2=3
        r1=0         
        do i1=1,ig1
        do m1=1,2
          z=zmin*(zmax/zmin)**(.5+tgss(ig1,i1)*(m1-1.5))
          xpmin=max(z,xpmin0)        
          r2=0  
          if(xpmin.lt.xpmax)then
          do i2=1,ig2
          do m2=1,2
            xp=xpmin*(xpmax/xpmin)**(.5+tgss(ig2,i2)*(m2-1.5))
            xm=z/xp
            r2=r2+wgss(ig2,i2)*ffsigiut(xp,xm,jaa,je1,je2)
          enddo
          enddo
          endif
          r2=r2*0.5*log(xpmax/xpmin)
          r1=r1+wgss(ig1,i1)*r2*z
        enddo
        enddo
        r1=r1*0.5*log(zmax/zmin)
        res=  r1 * factk * .0390  /sigine*10 
        sum=res
   !.......plot
        xx2min = 0.01/engy     !max(xpar1,0.01/engy)
        xx2max = 1             !xpar2
        xx1min = 0.01/engy     !max(xpar3,0.01/engy)
        xx1max = 1             !xpar4
        nbins  = 10            !nint(xpar5)

        write(ifhi,'(a,2i1)') 'openhisto xrange 0 1 name ffsig',jaa,jex           
        write(ifhi,'(a)') 'yrange auto auto htyp lin xmod lin ymod log'
        write(ifhi,'(a)') 'txt "xaxis x+?IB!         "              '
        write(ifhi,'(a)') 'txt "yaxis dn?semi! / dx+?IB!  "'
        write(ifhi,'(a,2i1,a)')'txt "title ffsig + MC   (',jaa,jex,')"'
        write(ifhi,'(a)') 'array 2'
        del=(xx1max-xx1min)/nbins
        do ii=1,nbins
          xx1=xx1min+(ii-0.5)*del
          ig2=3
          r2=0
          do i2=1,ig2
          do m2=1,2
            xx2=xx2min*(xx2max/xx2min)**(.5+tgss(ig2,i2)*(m2-1.5))
            r2=r2+wgss(ig2,i2)*ffsigiut(xx1,xx2,jaa,je1,je2)*xx2
          enddo
          enddo
          sig=r2*0.5*log(xx2max/xx2min)
          sig   = sig * factk * .0390   /sigine*10 
          write(ifhi,'(2e12.4)')xx1,sig
        enddo
        write(ifhi,'(a)')  '  endarray'
          
       elseif(iii.eq.4)then

        write(ifhi,'(a)')       '!----------------------------------'
        write(ifhi,'(a,3i1)')   '!   OB    ',jaa,jex
        write(ifhi,'(a)')       '!----------------------------------'

      !...... integral 
        y2     = 10
        ptmin  = 2
        ptmax  = 6
        sum=0   
        ig=2
        do i=1,ig
        do m=1,2
              pt=ptmin*(ptmax/ptmin)**(.5+tgss(ig,i)*(m-1.5))
          sig=ffsigi(pt**2,y2)     
          sig   =sig    * factk * .0390 /sigine*10  * 2   ! 2 partons!
              sum=sum+wgss(ig,i)*sig*pt 
        enddo
        enddo
        sum=sum*0.5*log(ptmax/ptmin)
      !...... pt distr      
        y2     = 10
        ptmin  = 2
        ptmax  = 20
        nbins  = 18
        sx=engy**2
        do jj=3,1,-1
        write(ifhi,'(a,i1)')'openhisto name jet',jj
        write(ifhi,'(a)')'xrange 0 20 xmod lin ymod log '
        write(ifhi,'(a)') 'txt "xaxis pt?OB!         "           '
        write(ifhi,'(a)') 'txt "yaxis dn?ptn! / dpt?OB!  "'
        if(jj.eq.1)write(ifhi,'(a)')'htyp lro'
        if(jj.eq.2)write(ifhi,'(a)')'htyp lgo'
        if(jj.eq.3)write(ifhi,'(a)')'htyp lyo'
        write(ifhi,'(a,f7.2,a)')  'text 0.05 0.1 "1/f=',1./ff,'"'
        write(ifhi,'(a)')'array 2'
        delpt=(ptmax-ptmin)/nbins
        do i=1,nbins
          pt=ptmin+(i-0.5)*delpt
          sig=1 
          if(jj.eq.1)then
            sig=ffsigi(pt**2,y2)      ! our stuff
          elseif(jj.eq.2)then
            if(engy.ge.10.)sig=psjvrg1(pt**2,sx,y2) ! grv
          elseif(jj.eq.3)then
            if(engy.ge.10.)sig=psjwo1(pt**2,sx,y2)   !duke-owens
          endif
          sig   =sig    * factk * .0390 /sigine*10 * 2 
          write(ifhi,'(2e12.4)')pt,sig
        enddo
        write(ifhi,'(a)')       '  endarray'
        if(jj.ne.1)write(ifhi,'(a)')       'closehisto'
        if(jj.ne.1)write(ifhi,'(a)')  'plot 0-'
        enddo

       endif
       
       x=0.1+(min(3,iii)-2)*0.30
       y=0.2+(min(3,iii)-2)*0.55
       if(engy.gt.100.)then
       write(ifhi,'(a,2f5.2,a,f6.3,a)')'text',x,y,' "   form ',sum,'"'
       else
       write(ifhi,'(a,2f5.2,a,f6.5,a)')'text',x,y,' "   form ',sum,'"'
       endif
       write(ifhi,'(a)')  'closehisto plot 0-'
       
       write(ifhi,'(a)') "!-----------------------------"
       write(ifhi,'(a)') "! MC   "
       write(ifhi,'(a)') "!-----------------------------"
       
       if(iii.eq.2)
     *  write(ifhi,'(a,i1,i1)')'openhisto name dndxPE',jaa,jex
       if(iii.eq.3)
     *  write(ifhi,'(a,i1,i1)')'openhisto name dndxIB',jaa,jex
       if(iii.eq.4)
     *  write(ifhi,'(a,i1,i1)')'openhisto name dndptOB',jaa,jex
       write(ifhi,'(a)')     'htyp prs'
       write(ifhi,'(a)')     'array 2'
       sum=0
       imax=nbixp
       if(iii.eq.4)imax=nbipt
       do i=1,imax
        u=xlp(kkk,i)
        if(iii.eq.4)u=ptob(i)
        z=0
        do j=j1,j2
        do jaai=jaa1,jaa2
        do jexi=jex1,jex2
         if(iii.eq.2)z=z+wxe(kkk,jaai,jexi,i,j)
         if(iii.eq.3)z=z+wxb(kkk,jaai,jexi,i,j)
         if(iii.eq.4)z=z+wptob(jaai,jexi,i)
        enddo
        enddo
        enddo
        del=dxlp(kkk,i)
        if(iii.eq.4)del=(pto-ptu)/nbipt
        z=z/del*ff/nrevt 
        write(ifhi,'(2e11.3)')u,z    
        sum=sum+z*del  
       enddo
       write(ifhi,'(a)')    '  endarray'
       x=0.1+(min(3,iii)-2)*0.30
       y=0.1+(min(3,iii)-2)*0.55
       if(engy.gt.100)then
       write(ifhi,'(a,2f5.2,a,f6.3,a)')'text',x,y,' "   simu ',sum,'"'
       else
       write(ifhi,'(a,2f5.2,a,f6.5,a)')'text',x,y,' "   simu ',sum,'"'
       endif
       write(ifhi,'(a)')    'closehisto'
       
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsSe(iii,xmc,ptmc,ih,iqq)
c-----------------------------------------------------------------------
c     iqq = 1 : String End mass and rapidity
c     iqq = 2 : String mass and rapidity
c-----------------------------------------------------------------------

      include 'epos.inc'

      parameter(nbix=50)
      common/cxpar/nx(2),x(nbix),wxmc(nbix,2),xmn,xmx,xu,xo
      parameter(nbiy=40)
      common/cypar/ny(2),y(nbiy),wymc(nbiy,2),ymin,ymax,dy,yu,yo 

      s=engy**2

      if(iii.eq.0)then
     
       nx(iqq)=0
       xu=0.1/engy**2
       xo=1.
       do i=1,nbix
         x(i)=xu*(xo/xu)**((i-0.5)/nbix)
         wxmc(i,iqq)=0
       enddo
       yo=alog(s)
       yu=-yo
       dy=(yo-yu)/nbiy
       ny(iqq)=0
       do i=1,nbiy
         y(i)=yu+dy/2.+(i-1)*dy
         wymc(i,iqq)=0
       enddo
       
      elseif(iii.eq.1)then
      
       if(xmc.lt.xu)return
       if(ptmc.eq.0.)return
       if(iqq.eq.1)ymc=0.5*alog(xmc*s/ptmc)*ih
       if(iqq.eq.2)ymc=0.5*alog(xmc/ptmc)
       i=1+int(alog(xmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxmc(i,iqq)=wxmc(i,iqq)+1
       nx(iqq)=nx(iqq)+1
1      continue
       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)return
       if(i.lt.1)return
       wymc(i,iqq)=wymc(i,iqq)+1
       ny(iqq)=ny(iqq)+1
       
      elseif(iii.eq.2)then
      
       write(ifhi,'(a)')        '!--------------------------------'
       write(ifhi,'(a)')        '!   string end x distr       '
       write(ifhi,'(a)')        '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(iqq.eq.1)write(ifhi,'(a)')    'text 0 0 "xaxis string end x"'
        if(iqq.eq.2)write(ifhi,'(a)')    'text 0 0 "xaxis string x"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nx(iqq).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxmc(i,iqq)/dx/nx(iqq)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod lin ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange',yu,yo
        if(iqq.eq.1)write(ifhi,'(a)')    'text 0 0 "xaxis string end y"'
        if(iqq.eq.2)write(ifhi,'(a)')    'text 0 0 "xaxis string y"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbiy
         if(ny(iqq).gt.0)
     *   write(ifhi,'(2e11.3)')y(i),wymc(i,iqq)/dy/ny(iqq)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine xEmsDr(iii,xpmc,xmmc,ie)
c-----------------------------------------------------------------------

      include 'epos.inc'

      parameter(nbix=50,nie=4)
      common/cxpardr/nxp(nie),nxm(nie),x(nbix),wxpmc(nbix,nie)
     &      ,wxmmc(nbix,nie),xmn,xmx,xu,xo,wxmc(nbix,nie),nx(nie)
      parameter(nbiy=40)
      common/cypardr/ny(nie),y(nbiy),wymc(nbiy,nie),ymin,ymax,dy,yu,yo 

      s=engy**2

      if(iii.eq.0)then
      
       do ni=1,nie
         nxp(ni)=0
         nxm(ni)=0
         nx(ni)=0
       enddo
       xu=0.1/engy**2
       xo=1.
       do i=1,nbix
         x(i)=xu*(xo/xu)**((i-0.5)/nbix)
         do ni=1,nie
           wxpmc(i,ni)=0
           wxmmc(i,ni)=0
           wxmc(i,ni)=0
         enddo
       enddo
       yo=alog(s)
       yu=-yo
       dy=(yo-yu)/nbiy
       do ni=1,nie
         ny(ni)=0
       enddo
       do i=1,nbiy
         y(i)=yu+dy/2.+(i-1)*dy
         do ni=1,nie
           wymc(i,ni)=0
         enddo
       enddo
       
      elseif(iii.eq.1)then

       if(ie.lt.1.or.ie.gt.nie)return
      
       if(xpmc.lt.xu)return
       i=1+int(alog(xpmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxpmc(i,ie)=wxpmc(i,ie)+1
       nxp(ie)=nxp(ie)+1
       if(xmmc.lt.xu)return
       i=1+int(alog(xmmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)goto1
       if(i.lt.1)goto1
       wxmmc(i,ie)=wxmmc(i,ie)+1
       nxm(ie)=nxm(ie)+1
1      continue
       if(xmmc.ge.xu)then
         ymc=0.5*alog(xpmc/xmmc)
       else
         return
       endif
       if(ymc.lt.yu)return
       i=int((ymc-yu)/dy)+1
       if(i.gt.nbiy)return
       if(i.lt.1)return
       wymc(i,ie)=wymc(i,ie)+1
       ny(ie)=ny(ie)+1

       xmc=xpmc*xmmc
       if(xmc.lt.xu)return
       i=1+int(alog(xmc/xu)/alog(xo/xu)*nbix)
       if(i.gt.nbix)return
       if(i.lt.1)return
       wxmc(i,ie)=wxmc(i,ie)+1
       nx(ie)=nx(ie)+1
       
      elseif(iii.eq.2)then
     
        do ii=1,nie

       if(ii.eq.1)write(ifhi,'(a)')'!-----  projectile droplet  ----'
       if(ii.eq.2)write(ifhi,'(a)')'!-----    target droplet    ----'
       if(ii.eq.3)write(ifhi,'(a)')'!-----  projectile string end  ----'
       if(ii.eq.4)write(ifhi,'(a)')'!-----    target string end    ----'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string x+ distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lru'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet x+"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end x+"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nxp(ii).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxpmc(i,ii)/dx/nxp(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0-'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string x- distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lba'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',xu,xo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet x-"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end x-"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nxm(ii).gt.0)
     *   write(ifhi,'(2e11.3)')x(i),wxmmc(i,ii)/dx/nxm(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string y distr       '
        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       'openhisto'
        write(ifhi,'(a)')       'htyp lin'
        write(ifhi,'(a)')       'xmod lin ymod lin'
        write(ifhi,'(a,2e11.3)')'xrange',yu,yo
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet y"'
        if(ii.eq.3.or.ii.eq.4)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end y"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(y)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbiy
         if(ny(ii).gt.0)
     *   write(ifhi,'(2e11.3)')y(i),wymc(i,ii)/dy/ny(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    'closehisto plot 0'
        
      enddo

        write(ifhi,'(a)')       '!--------------------------------'
        write(ifhi,'(a)')       '!   droplet/string mass distr       '
        write(ifhi,'(a)')       '!--------------------------------'
      do ii=1,nie


        if(ii.eq.2.or.ii.eq.4)write(ifhi,'(a)')    'closehisto plot 0-'
        if(ii.eq.3)write(ifhi,'(a)')    'closehisto plot 0'
        write(ifhi,'(a)')       'openhisto'
        if(ii.eq.1.or.ii.eq.3)write(ifhi,'(a)')       'htyp lru'
        if(ii.eq.2.or.ii.eq.4)write(ifhi,'(a)')       'htyp lba'
        write(ifhi,'(a)')       'xmod log ymod log'
        write(ifhi,'(a,2e11.3)')'xrange',sqrt(xu*s),sqrt(s*xo)
        if(ii.eq.1.or.ii.eq.2)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis droplet mass (GeV)"'
        if(ii.eq.4.or.ii.eq.3)
     *  write(ifhi,'(a)')    'text 0 0 "xaxis string end mass (GeV)"'
        write(ifhi,'(a)')    'text 0 0 "yaxis P(x)"'
        write(ifhi,'(a)')       'array 2'
        do i=1,nbix
         dx=xu*(xo/xu)**(1.*i/nbix)*(1.-(xo/xu)**(-1./nbix))
         if(nx(ii).gt.0)
     *   write(ifhi,'(2e11.3)')sqrt(x(i)*s),wxmc(i,ii)/dx/nx(ii)
        enddo
        write(ifhi,'(a)')    '  endarray'
      enddo
       write(ifhi,'(a)')    'closehisto plot 0'

      endif

      return
      end

cc--------------------------------------------------------------------------
c      subroutine xtype(k,n,i1,i2,text)
cc--------------------------------------------------------------------------
c      
c      include 'epos.inc'
c      include 'epos.incems'
c      parameter(itext=40)
c      character  text*40
c
c      imax=itext+1
c      do i=itext,1,-1
c      if(text(i:i).eq.'&')imax=i
c      enddo
c      
c      ip=iproj(k)
c      it=itarg(k)
c       
c      if(i1.eq.1)then
c         write(ifch,*)
c         write(ifch,*)('-',ll=1,27)
c         write(ifch,*)'  '//text(1:imax-1)
c         write(ifch,*)('-',ll=1,27)
c      endif
c      
c      if(i2.eq.1)then
c         write(ifch,*)
c         write(ifch,*)'k:',k,'   n:',n,'   ip:',ip,'   it:',it
c         write(ifch,*)'bk:',bk(k)
c         if(n.ne.0)write(ifch,*)'idpr:',idpr(n,k)
c         write(ifch,*)'iep:',iep(ip),'   iet:',iet(it)
c         write(ifch,*)'idp:',idp(ip),'   idt:',idt(it)
c      endif
c      
c      end
c
c------------------------------------------------------------------------       
      subroutine XPrint(text) 
c------------------------------------------------------------------------  
      include 'epos.inc'
      include 'epos.incems'
      double precision xpptot,xmptot,xpttot,xmttot
      parameter(itext=15)
      character  text*15
      imax=itext+1
      do i=itext,1,-1
      if(text(i:i).eq.'&')imax=i
      enddo
      write(ifch,'(1x,a)')text(1:imax-1)
      
      write(ifch,'(a)')' npr0:   npr1:   nprmx:   Pomeron id lattice:' 
      do k=1,koll
       write(ifch,'(1x,i6,1x,i2,6x,i2,6x,i2,7x,$)')
     *              k,npr(0,k),npr(1,k),nprmx(k)
       do n=1,nprmx(k)
        write(ifch,'(i2,$)')idpr(n,k)
       enddo
       write(ifch,*)' '
      enddo
      
      xpptot=0d0
      xmptot=0d0
      xpttot=0d0
      xmttot=0d0
      write(ifch,'(a)')' Pomeron xy lattice:' 
      do k=1,koll
       do n=1,nprmx(k)
       xpptot=xpptot+xppr(n,k)
       xmttot=xmttot+xmpr(n,k)
        write(ifch,'(i6,1x,i2,1x,d10.3,1x,d10.3,3x,$)')
     *                  k,n,xpr(n,k),ypr(n,k)
       enddo
       write(ifch,*)' '
      enddo
      
      write(ifch,'(a)')' projectile remnants x+,x-,px,py,x,iep:' 
      do ip=1,maproj
       xpptot=xpptot+xpp(ip)
       xmptot=xmptot+xmp(ip)
       write(ifch,'(i3,2x,5d12.3,i3)')ip,xpp(ip),xmp(ip),xxp(ip),xyp(ip)
     *                             ,xpos(ip),iep(ip)
      enddo
      
      write(ifch,'(a)')' target remnants x-,x+,px,py,x,iet:' 
      do it=1,matarg
       xpttot=xpttot+xpt(it)
       xmttot=xmttot+xmt(it)
       write(ifch,'(i3,2x,5d12.3,i3)')it,xmt(it),xpt(it),xxt(it),xyt(it)
     *                             ,xtos(it),iet(it)
      enddo
      
      write(ifch,*)' remnant balance x+,x-:'
     &,(xpptot+xpttot)/dble(maproj)
     &,(xmptot+xmttot)/dble(matarg)
      end


c-------------------------------------------------------------------------
      subroutine xfom
c-------------------------------------------------------------------------
      include 'epos.inc'
      double precision fom,x
      write(ifhi,'(a)')     '!##################################'
      write(ifhi,'(a,i3)')  '!   fom     '
      write(ifhi,'(a)')     '!##################################'
      b=0.
      do i=1,6
        z=0.2*exp(0.8*i)
        xi=0.01+0.16*float(i-1)
        write(ifhi,'(a,i1)') 'openhisto name fom',i
        write(ifhi,'(a)')    'htyp lin xmod lin ymod log'
        write(ifhi,'(a)')    'xrange 0 1'
        write(ifhi,'(a)')    'yrange 0.1 1000 '
        write(ifhi,'(a)')    'text 0 0 "xaxis x "'
        write(ifhi,'(a)')    'text 0 0 "yaxis fom"'
        if(z.lt.10.)
     &   write(ifhi,'(a,f4.2,a,f4.1,a)')'text ',xi,' 0.9 "',z,'"'
        if(z.ge.10.)
     &   write(ifhi,'(a,f4.2,a,f4.0,a)')'text ',xi,' 0.9 "',z,'"'
        write(ifhi,'(a)')    'array 2'
        do n=1,99
          x=dble(n)*0.01d0    
          write(ifhi,'(2e11.3)')x,fom(z,x,b)
        enddo
        write(ifhi,'(a)')    '  endarray'
        write(ifhi,'(a)')    '  closehisto '
        if(i.lt.6)write(ifhi,'(a)')    'plot 0-'
        if(i.eq.6)write(ifhi,'(a)')    'plot 0'
      enddo
      end            


