! -*- F90 -*-


      subroutine evolvePDF(x,Q,f) 
      implicit none 
      include 'parmsetup.inc'
      real*8 gridx(nmxgridx),gridq(nmxgridq)
      integer nxgrid,nqgrid
      integer nset,imem,Eorder,IP2 
      real*8 x,Q,P2,Q2fit,f(-6:6),alfas,a,photon,gluino 
      nset = 1 
      call evolvePDFM(nset,x,Q,f) 
      return 
!                                                                       
      entry evolvePDFp(x,Q,P2,IP2,f) 
      nset = 1 
      call evolvePDFpM(nset,x,Q,P2,IP2,f) 
      return 
!                                                                       
      entry evolvePDFa(x,Q,a,f) 
      nset = 1 
      call evolvePDFaM(nset,x,Q,a,f) 
      return 
!                                                                       
      entry evolvePDFphoton(x,Q,f,photon) 
      nset = 1 
      call evolvePDFphotonM(nset,x,Q,f,photon) 
      return 
!                                                                       
      entry evolvePDFgluino(x,Q,f,gluino) 
      nset = 1 
      call evolvePDFgluinoM(nset,x,Q,f,gluino) 
      return 
!                                                                       
      entry initPDF(imem) 
      nset = 1 
      call initPDFM(nset,imem) 
      return

      entry getgrid(nxgrid,nqgrid,gridx,gridq)
      nset = 1
      call getgridM(nset,nxgrid,nqgrid,gridx,gridq)

      return 
      END                                           
!                                                                       
      subroutine evolvePDFaM(nset,xin,Qin,a,f) 
      implicit none 
      real*8 x,Q,a,f(-6:6) 
      real*8 ruv,rdv,ru,rd,rs,rc,rb,rt,rg 
      integer nset,iimem,j 
      real*8 xin,qin,q2in 
      character*20 lparm 
      real*8 xmin,xmax,q2min,q2max 
      integer iorder,ipset,ia
     
                                                                        
      call getlhaparm(18,lparm) 
      if(lparm.ne.'EXTRAPOLATE') then 
        call getnmem(nset,iimem) 
        call getminmaxm(nset,iimem,xmin,xmax,q2min,q2max) 
        x = max(xmin,min(xmax,xin)) 
        q2in = qin**2 
        q = sqrt(max(0d0,q2min,min(q2max,q2in))) 
      else 
        x = xin 
        q = qin 
      endif 
                                                                        
      call getlhaparm(15,lparm) 
      if(lparm.eq.'EPS08') then 
        call eps08(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
      else if(lparm(1:5).eq.'EPS09') then 
        if(lparm.eq.'EPS09LO') then 
          iorder=1
          ipset=1 
        else if(lparm.eq.'EPS09NLO') then 
          iorder=2
          ipset=1 
        else if(lparm(1:8).eq.'EPS09LO,') then 
          iorder=1
          read(lparm(9:),*)ipset 
        else if(lparm(1:9).eq.'EPS09NLO,') then 
          iorder=2
          read(lparm(10:),*)ipset 
        else
          iorder=2
          ipset=1 
        endif
        ia = a
        call eps09(iorder,ipset,ia,x,q,ruv,rdv,ru,rd,rs,rc,rb,rg)
        rt = 1.0d0
      else
        call eks98(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
      endif
      
      call evolvePDFM(nset,x,Q,f)
      f(0) = f(0)*rg
      f(1) = f(1)*rdv-f(-1)*(rdv-rd)
      f(2) = f(2)*ruv-f(-2)*(ruv-ru)
      f(3) = f(3)*rs      
      f(4) = f(4)*rc      
      f(5) = f(5)*rb      
      f(6) = f(6)*rt      
      f(-6) = f(-6)*rt
      f(-5) = f(-5)*rb
      f(-4) = f(-4)*rc
      f(-3) = f(-3)*rs
      f(-2) = f(-2)*ru
      f(-1) = f(-1)*rd

      return
      end
!      
      subroutine evolvePDFM(nset,xin,Qin,f)
      implicit none
      include 'parmsetup.inc'
      integer Eorder,index,imem
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem,ip2
      common/NAME/name,nmem,ndef,mem
      integer iset,iimem
      common/SET/iset,iimem
      integer nset,j
      real*8 x,xin,Q,Qin,Q2fit,alfas,p2,q2in
      real*8 f(-6:6),photon,gluino,xphoton
      character*20 lparm
      real*8 xmin,xmax,q2min,q2max
      real*8 gridx(nmxgridx),gridq(nmxgridq)
      integer nxgrid,nqgrid
      character*512 setpath
      integer nnpdf,nnpdf100,nnpdf1000
       data nnpdf,nnpdf100,nnpdf1000/0,0,0/
      save
!
      call setnset(nset)
!      
!      print *,'this is evolvePDFM, name=',nset,name(nset)
!   set all f's to 0.0d0 at start
!      do j = -6,6
!        f(j) = 0.0d0
!      enddo 

      call getlhaparm(18,lparm)
      if(lparm.ne.'EXTRAPOLATE') then
        call getnmem(nset,iimem)
        call getminmaxm(nset,iimem,xmin,xmax,q2min,q2max)
        x = max(xmin,min(xmax,xin))
        q2in = qin**2
        q = sqrt(max(0d0,q2min,min(q2max,q2in)))
      else
        x = xin
        q = qin
      endif

#ifdef QCDNUM 
      if (name(nset).eq.'QCDNUM') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3evolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4evolve(x,Q,f)
#endif
#ifdef GRV
      if (name(nset).eq.'GRV') call GRVevolve(x,Q,f)
#endif
#ifdef ZEUS     
      if (name(nset)(1:12).eq.'QCDNUM_ZEUS_') call ZEUSevolve(x,Q,f)
#endif
#ifdef HERA
      if (name(nset)(1:12).eq.'QCDNUM_HERA_') call HERAevolve(x,Q,f)
      if (name(nset)(1:8).eq.'HERAGRID') call HERAGRIDevolve(x,Q,f)
#endif
#ifdef CTEQ
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQevolve(x,Q,f)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5evolve(x,Q,f)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6evolve(x,Q,f)
      if (name(nset).eq.'CTEQ65grid') call CTEQ65evolve(x,Q,f)
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65cevolve(x,Q,f)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65sevolve(x,Q,f)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6evolve(x,Q,f)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65evolve(x,Q,f)
      if (name(nset).eq.'CT10grid') call CT12evolve(x,Q,f)
      if (name(nset).eq.'CT12grid') call CT12evolve(x,Q,f)
#endif
#ifdef MRST
      if (name(nset).eq.'MRST') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'MRSTpdf') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'MRSTgrid') call MRSTevolve(x,Q,f)
      if (name(nset).eq.'MRST3grid') call MRSTevolve(x,Q,f)
      if (name(nset).eq.'MRST4grid') call MRSTevolve(x,Q,f)
#endif
#ifdef MRSTQED
      if (name(nset).eq.'MRST4qed') call MRSTqedevolve(x,Q,f,xphoton)
#endif
#ifdef MRST98
      if (name(nset).eq.'MRST98grid') call MRST98evolve(x,Q,f)
#endif
#ifdef MRST06
      if (name(nset).eq.'MRST2006grid') call MRST2006evolve(x,Q,f)
#endif
#ifdef MSTW
      if (name(nset).eq.'MSTWgrid') call MSTWevolve(x,Q,f,xphoton)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Mevolve(x,Q,f)
      if (name(nset).eq.'ABKM09') call ABKM09evolve(x,Q,f)
      if (name(nset).eq.'ABM11') call ABM11evolve(x,Q,f)
#endif
#ifdef H1
      if (name(nset).eq.'H12000') call H1evolve(x,Q,f)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRevolve(x,Q,f)
#endif
#ifdef NNPDF
      if (name(nset).eq.'NNPDF') call NNPDFevolve(x,Q,f)
      if (name(nset).eq.'NNPDFint') call NNPDFINTevolve(x,Q,f)
      if (name(nset).eq.'NNPDF20int') call NNPDFINT20evolve(x,Q,f)
      if (name(nset).eq.'NNPDF20intqed') call NNPDFINT20qedevolve(x,Q,f,xphoton)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hknevolve(x,Q,f)
#endif
#ifdef PIONS
      if (name(nset).eq.'OWP') call OWPevolve(x,Q,f)
      if (name(nset).eq.'SMRSP') call SMRSPevolve(x,Q,f)
      if (name(nset).eq.'GRVP0') call GRVP0evolve(x,Q,f)
      if (name(nset).eq.'GRVP1') call GRVP1evolve(x,Q,f)
      if (name(nset).eq.'ABFKWP') call ABFKWPevolve(x,Q,f)
#endif
#ifdef USER
      if (name(nset).eq.'USER') call USERevolve(x,Q,f)
      if (name(nset)(1:8).eq.'USERGRID') call USERGRIDevolve(x,Q,f)
#endif
      return
!
      entry evolvePDFpM(nset,xin,Qin,P2,IP2,f)
!
      call setnset(nset)
!      
      call getlhaparm(18,lparm)
      if(lparm.ne.'EXTRAPOLATE') then
        call getnmem(nset,iimem)
        call getminmaxm(nset,iimem,xmin,xmax,q2min,q2max)
        x = max(xmin,min(xmax,xin))
        q2in = qin**2
        q = sqrt(max(0d0,q2min,min(q2max,q2in)))
      else
        x = xin
        q = qin
      endif
!
#ifdef PHOTONS
      if(name(nset).eq.'SASG') call SASGevolvep(x,Q,P2,IP2,f)
      if(name(nset).eq.'GRVG0') call GRVGevolvep0(x,Q,P2,IP2,f)
      if(name(nset).eq.'GRVG1') call GRVGevolvep1(x,Q,P2,IP2,f)
      if (name(nset).eq.'DOG0') call DOGevolvep0(x,Q,P2,IP2,f)
      if (name(nset).eq.'DOG1') call DOGevolvep1(x,Q,P2,IP2,f)
      if (name(nset).eq.'DGG') call DGGevolvep(x,Q,P2,IP2,f)
      if (name(nset).eq.'LACG') call LACGevolvep(x,Q,P2,IP2,f)
      if (name(nset).eq.'GSG0') call GSGevolvep0(x,Q,P2,IP2,f)
      if (name(nset).eq.'GSG1') call GSGevolvep1(x,Q,P2,IP2,f)
      if (name(nset).eq.'GSG960') call GSG96evolvep0(x,Q,P2,IP2,f)
      if (name(nset).eq.'GSG961') call GSG96evolvep1(x,Q,P2,IP2,f)
      if (name(nset).eq.'ACFGP') call ACFGPevolvep(x,Q,P2,IP2,f)
      if (name(nset).eq.'WHITG') call WHITGevolvep(x,Q,P2,IP2,f)
#endif
      return
!
      entry evolvePDFphotonM(nset,xin,qin,f,photon)      

!
      call setnset(nset)
!
      call getlhaparm(18,lparm)
      if(lparm.ne.'EXTRAPOLATE') then
        call getnmem(nset,iimem)
        call getminmaxm(nset,iimem,xmin,xmax,q2min,q2max)
        x = max(xmin,min(xmax,xin))
        q2in = qin**2
        q = sqrt(max(0d0,q2min,min(q2max,q2in)))
      else
        x = xin
        q = qin
      endif
!
#ifdef MRSTQED      
      if(name(nset).eq.'MRST4qed') then
        call MRSTqedevolve(x,Q,f,photon)
      else if (name(nset).ne.'NNPDF20intqed') then
        photon = 0.0d0
      endif
#endif
!
#ifdef NNPDF      
      if(name(nset).eq.'NNPDF20intqed') then
        call NNPDFINT20qedevolve(x,Q,f,photon)
      else if (name(nset).ne.'MRST4qed') then
        photon = 0.0d0
      endif
#endif
!
      return

      entry evolvePDFgluinoM(nset,xin,qin,f,gluino)      

!
      call setnset(nset)
!
      call getlhaparm(18,lparm)
      if(lparm.ne.'EXTRAPOLATE') then
        call getnmem(nset,iimem)
        call getminmaxm(nset,iimem,xmin,xmax,q2min,q2max)
        x = max(xmin,min(xmax,xin))
        q2in = qin**2
        q = sqrt(max(0d0,q2min,min(q2max,q2in)))
      else
        x = xin
        q = qin
      endif
!
#ifdef CTEQ
      if(name(nset).eq.'CTEQ6LGgrid') then
        call CTEQ6LGevolve(x,Q,f,gluino)
      else
        gluino = 0.0d0
      endif
#endif
      return


      entry readevolve(nset)
!
      call getsetpath(setpath)
!
      if(index(setpath,'NNPDF').gt.0) then
          if(nset.gt.1) then
              if((index(setpath,'1000.').gt.0.and.nnpdf.gt.0).or.(index(setpath,'100.').gt.0.and.nnpdf1000.gt.0)) then
                  print *,'LHAPDF ERROR: MULTISET-INITIALIZATION with NNPDF 1000 sets IS NOT AVAIABLE (AT THE MOMENT)!'
                  STOP
              endif
          endif
          nnpdf=nnpdf+1
          if(index(setpath,'1000').gt.0) then
              nnpdf1000=nnpdf1000+1
          else
              nnpdf100=nnpdf100+1
          endif
      endif

!
      read(1,*) name(nset)   
!      print *, 'this is readevolve', name(nset)   
!
      call setnset(nset)
!      
#ifdef QCDNUM
      if (name(nset).eq.'QCDNUM') call QCDNUMread(nset)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMread(nset)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3read(nset)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4read(nset)
#endif
#ifdef GRV
      if (name(nset).eq.'GRV') call GRVread(nset)
#endif
#ifdef ZEUS
      if (name(nset)(1:12).eq.'QCDNUM_ZEUS_') call ZEUSread(nset)
#endif
#ifdef HERA
      if (name(nset)(1:12).eq.'QCDNUM_HERA_') call HERAread(nset)
      if (name(nset)(1:8).eq.'HERAGRID') call HERAGRIDread(nset)
#endif
#ifdef CTEQ
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQread(nset)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5read(nset)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6read(nset)
      if (name(nset).eq.'CTEQ65grid') call CTEQ65read(nset)
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65cread(nset)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65read(nset)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6read(nset)
      if (name(nset).eq.'CTEQ66grid') call CTEQ66read(nset)
      if (name(nset).eq.'CTEQ6LGgrid') call CTEQ6LGread(nset)
      if (name(nset).eq.'CT10grid') call CT12read(nset)
      if (name(nset).eq.'CT12grid') call CT12read(nset)
#endif
#ifdef MRST
      if (name(nset).eq.'MRST') call QCDNUMread(nset)
      if (name(nset).eq.'MRSTpdf') call QCDNUMread(nset)
      if (name(nset).eq.'MRSTgrid') call MRSTread(nset)
      if (name(nset).eq.'MRST3grid') call MRSTread(nset)
      if (name(nset).eq.'MRST4grid') call MRSTread(nset)
#endif
#ifdef MRSTQED
      if (name(nset).eq.'MRST4qed') call MRSTqedread(nset)
#endif
#ifdef MRST98
      if (name(nset).eq.'MRST98grid') call MRST98read(nset)
#endif
#ifdef MRST06
      if (name(nset).eq.'MRST2006grid') call MRST2006read(nset)
#endif
#ifdef MSTW
      if (name(nset).eq.'MSTWgrid') call MSTWread(nset)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Mread(nset)
      if (name(nset).eq.'ABKM09') call ABKM09read(nset)
      if (name(nset).eq.'ABM11') call ABM11read(nset)
#endif
#ifdef H1      
      if (name(nset).eq.'H12000') call H1read(nset)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRread(nset)
#endif
#ifdef NNPDF
      if (name(nset).eq.'NNPDF') call NNPDFread(nset)
      if (name(nset).eq.'NNPDFint') call NNPDFINTread(nset)
      if (name(nset).eq.'NNPDF20int') call NNPDFINT20read(nset)
      if (name(nset).eq.'NNPDF20intqed') call NNPDFINT20qedread(nset)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hknread(nset)
#endif
#ifdef PHOTONS
      if (name(nset).eq.'SASG') call SASGread(nset)
      if (name(nset).eq.'GRVG0' .OR. &
     &    name(nset).eq.'GRVG1') call GRVGread(nset)
      if (name(nset).eq.'DOG0' .OR. &
     &    name(nset).eq.'DOG1') call DOGread(nset)
      if (name(nset).eq.'DGG') call DGGread(nset)
      if (name(nset).eq.'LACG') call LACGread(nset)
      if (name(nset).eq.'GSG0' .OR. & 
     &    name(nset).eq.'GSG1') call GSGread(nset)
      if (name(nset).eq.'GSG960' .OR. &
     &    name(nset).eq.'GSG961') call GSG96read(nset)
      if (name(nset).eq.'ACFGP') call ACFGPread(nset)
      if (name(nset).eq.'WHITG') call WHITGread(nset)
#endif
#ifdef PIONS
      if (name(nset).eq.'OWP') call OWPread(nset)
      if (name(nset).eq.'SMRSP') call SMRSPread(nset)
      if (name(nset).eq.'GRVP0' .OR. & 
     &    name(nset).eq.'GRVP1') call GRVPread(nset)
      if (name(nset).eq.'ABFKWP') call ABFKWPread(nset)
#endif
#ifdef USER
      if (name(nset).eq.'USER') call USERread(nset)
      if (name(nset)(1:8).eq.'USERGRID') call USERGRIDread(nset)
#endif
      return
!
      entry alfasevolve(nset,alfas,Qin)
!
      call setnset(nset)
        q = qin
#ifdef QCDNUM
      if (name(nset).eq.'QCDNUM') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3alfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4alfa(alfas,Q)
#endif
#ifdef GRV
      if (name(nset).eq.'GRV') call GRValfa(alfas,Q)
#endif
#ifdef ZEUS
      if (name(nset)(1:12).eq.'QCDNUM_ZEUS_') call ZEUSalfa(alfas,Q)
#endif
#ifdef HERA
      if (name(nset)(1:12).eq.'QCDNUM_HERA_') call HERAalfa(alfas,Q)
      if (name(nset)(1:8).eq.'HERAGRID') call HERAGRIDalfa(alfas,Q)
#endif
#ifdef CTEQ
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQalfa(alfas,Q)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5alfa(alfas,Q)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6alfa(alfas,Q)
      if (name(nset).eq.'CTEQ65grid') call CTEQ65alfa(alfas,Q)
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65alfa(alfas,Q)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65alfa(alfas,Q)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65alfa(alfas,Q)
      if (name(nset).eq.'CTEQ6LGgrid') call CTEQ6LGalfa(alfas,Q)
      if (name(nset).eq.'CT10grid') call CT12alfa(alfas,Q)
      if (name(nset).eq.'CT12grid') call CT12alfa(alfas,Q)
#endif
#ifdef MRST
      if (name(nset).eq.'MRST') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'MRSTpdf') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'MRSTgrid') call MRSTalfa(5,alfas,Q)
      if (name(nset).eq.'MRST3grid') call MRSTalfa(3,alfas,Q)
      if (name(nset).eq.'MRST4grid') call MRSTalfa(4,alfas,Q)
#endif
#ifdef MRSTQED
      if (name(nset).eq.'MRST4qed') call MRSTqedalfa(4,alfas,Q)
#endif
#ifdef MRST98
      if (name(nset).eq.'MRST98grid') call MRST98alfa(alfas,Q)
#endif
#ifdef MRST06
      if (name(nset).eq.'MRST2006grid') call MRST2006alfa(5,alfas,Q)
#endif
#ifdef MSTW
      if (name(nset).eq.'MSTWgrid') call MSTWalfa(alfas,Q)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Malfa(alfas,Q)
      if (name(nset).eq.'ABKM09') call ABKM09alfa(alfas,Q)
      if (name(nset).eq.'ABM11') call ABM11alfa(alfas,Q)
#endif
#ifdef H1
      if (name(nset).eq.'H12000') call H1alfa(alfas,Q)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRalfa(alfas,Q)
#endif
#ifdef NNPDF
      if (name(nset).eq.'NNPDF') call NNPDFalfa(alfas,Q)
      if (name(nset).eq.'NNPDFint') call NNPDFINTalfa(alfas,Q)
      if (name(nset).eq.'NNPDF20int') call NNPDFINT20alfa(alfas,Q)
      if (name(nset).eq.'NNPDF20intqed') call NNPDFINT20qedalfa(alfas,Q)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hknalfa(alfas,Q)
#endif
#ifdef PHOTONS
      if (name(nset).eq.'SASG') call SASGalfa(alfas,Q)
      if (name(nset).eq.'GRVG0' .OR. &
     &    name(nset).eq.'GRVG1') call GRVGalfa(alfas,Q)
      if (name(nset).eq.'DOG0' .OR. &
     &    name(nset).eq.'DOG1') call DOGalfa(alfas,Q)
      if (name(nset).eq.'DGG') call DGGalfa(alfas,Q)
      if (name(nset).eq.'LACG') call LACGalfa(alfas,Q)
      if (name(nset).eq.'GSG0' .OR. &
     &    name(nset).eq.'GSG1') call GSGalfa(alfas,Q)
      if (name(nset).eq.'GSG960' .OR. &
     &    name(nset).eq.'GSG961') call GSG96alfa(alfas,Q)
      if (name(nset).eq.'ACFGP') call ACFGPalfa(alfas,Q)
      if (name(nset).eq.'WHITG') call WHITGalfa(alfas,Q)
#endif
#ifdef PIONS
      if (name(nset).eq.'OWP') call OWPalfa(alfas,Q)
      if (name(nset).eq.'SMRSP') call SMRSPalfa(alfas,Q)
      if (name(nset).eq.'GRVP0' .OR. & 
     &    name(nset).eq.'GRVP1') call GRVPalfa(alfas,Q)
      if (name(nset).eq.'ABFKWP') call ABFKWPalfa(alfas,Q)
#endif
#ifdef USER
      if (name(nset).eq.'USER') call USERalfa(alfas,Q)
      if (name(nset)(1:8).eq.'USERGRID') call USERGRIDalfa(alfas,Q)
#endif
      return
!
      entry initevolution(nset,Eorder,Q2fit)
!
      call setnset(nset)
!            
#ifdef QCDNUM
      if (name(nset).eq.'QCDNUM') call QCDNUMinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'QCDNUM_MRST') then 
        call QCDNUMinit(nset,Eorder,Q2fit) 
        call QNLSET('BMARK',.TRUE.) 
      endif 
      if (name(nset).eq.'QCDNUM_MRST3') then 
        call QCDNUM3init(nset,Eorder,Q2fit) 
        call QNLSET('BMARK',.TRUE.) 
      endif 
      if (name(nset).eq.'QCDNUM_MRST4') then 
        call QCDNUM4init(nset,Eorder,Q2fit) 
        call QNLSET('BMARK',.TRUE.) 
      endif 
#endif
#ifdef GRV
      if (name(nset).eq.'GRV') call GRVinit(Eorder,Q2fit)
#endif
#ifdef ZEUS
      if (name(nset)(1:12).eq.'QCDNUM_ZEUS_') &
     & call ZEUSinit(nset,Eorder,Q2fit)
#endif
#ifdef HERA
      if (name(nset)(1:12).eq.'QCDNUM_HERA_') & 
     & call HERAinit(nset,Eorder,Q2fit)
      if (name(nset)(1:8).eq.'HERAGRID')call HERAGRIDinit(nset,Eorder,Q2fit)
#endif
#ifdef CTEQ
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ65grid') call CTEQ65init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ6LGgrid') call CTEQ6LGinit(Eorder,Q2fit)
      if (name(nset).eq.'CT10grid') call CT12init(Eorder,Q2fit)
      if (name(nset).eq.'CT12grid') call CT12init(Eorder,Q2fit)
#endif
#ifdef MRST
      if (name(nset).eq.'MRSTgrid') call MRSTinit(Eorder,Q2fit)
      if (name(nset).eq.'MRST3grid') call MRSTinit(Eorder,Q2fit)
      if (name(nset).eq.'MRST4grid') call MRSTinit(Eorder,Q2fit)
#endif
#ifdef MRSTQED
      if (name(nset).eq.'MRST4qed') call MRSTqedinit(Eorder,Q2fit)
#endif
#ifdef MRST98
      if (name(nset).eq.'MRST98grid') call MRST98init(Eorder,Q2fit)
#endif
#ifdef MRST06
      if (name(nset).eq.'MRST2006grid') call MRST2006init(Eorder,Q2fit)
#endif
#ifdef MSTW
      if (name(nset).eq.'MSTWgrid') call MSTWinit(Eorder,Q2fit)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Minit
      if (name(nset).eq.'ABKM09') call ABKM09init
      if (name(nset).eq.'ABM11') call ABM11init
#endif
#ifdef H1
      if (name(nset).eq.'H12000') call H1init(Eorder,Q2fit)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRinit(Eorder,Q2fit)
#endif
#ifdef NNPDF
      if (name(nset).eq.'NNPDF') call NNPDFinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'NNPDFint') call NNPDFINTinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'NNPDF20int') call NNPDFINT20init(nset,Eorder,Q2fit)
      if (name(nset).eq.'NNPDF20intqed') call NNPDFINT20qedinit(nset,Eorder,Q2fit)
#endif
#ifdef PHOTONS
      if (name(nset).eq.'SASG') call SASGinit(Eorder,Q2fit)
      if (name(nset).eq.'GRVG0' .OR. &
     &    name(nset).eq.'GRVG1') call GRVGinit(Eorder,Q2fit)
      if (name(nset).eq.'DOG0' .OR. &
     &    name(nset).eq.'DOG1') call DOGinit(Eorder,Q2fit)
      if (name(nset).eq.'DGG') call DGGinit(Eorder,Q2fit)
      if (name(nset).eq.'LACG') call LACGinit(Eorder,Q2fit)
      if (name(nset).eq.'GSG0' .OR. &
     &    name(nset).eq.'GSG1') call GSGinit(Eorder,Q2fit)
      if (name(nset).eq.'GSG960' .OR. &
     &    name(nset).eq.'GSG961') call GSG96init(Eorder,Q2fit)
      if (name(nset).eq.'ACFGP') call ACFGPinit(Eorder,Q2fit)
      if (name(nset).eq.'WHITG') call WHITGinit(Eorder,Q2fit)
#endif
#ifdef PIONS
      if (name(nset).eq.'OWP') call OWPinit(Eorder,Q2fit)
      if (name(nset).eq.'SMRSP') call SMRSPinit(Eorder,Q2fit)
      if (name(nset).eq.'GRVP0' .OR. & 
     &    name(nset).eq.'GRVP1') call GRVPinit(Eorder,Q2fit)
      if (name(nset).eq.'ABFKWP') call ABFKWPinit(Eorder,Q2fit)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hkninit(nset,Eorder,Q2fit)
#endif
#ifdef USER
      if (name(nset).eq.'USER') call USERinit(nset,Eorder,Q2fit)
      if (name(nset)(1:8).eq.'USERGRID') call USERGRIDinit(nset,Eorder,Q2fit)
#endif
      return
!
      entry initPDFM(nset,imem)
!
      call setnset(nset)
      call setnmem(nset,imem)
!            
      iimem = imem
#ifdef QCDNUM
      if (name(nset).eq.'QCDNUM') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'QCDNUM_MRST') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'QCDNUM_MRST3') then
         call InitEvolvePDF(nset,imem)
         call QCDNUM3pdf(nset)
      endif
      if (name(nset).eq.'QCDNUM_MRST4') then
         call InitEvolvePDF(nset,imem)
         call QCDNUM4pdf(nset)
      endif
#endif
#ifdef ZEUS
      if (name(nset)(1:12).eq.'QCDNUM_ZEUS_') then
         call InitEvolvePDF(nset,imem)
         call ZEUSpdf(nset)
      endif
#endif
#ifdef HERA
      if (name(nset)(1:12).eq.'QCDNUM_HERA_') then
         call InitEvolvePDF(nset,imem)
         call HERApdf(nset)
      endif
      if (name(nset)(1:8).eq.'HERAGRID') call HERAGRIDpdf(imem)
#endif
#ifdef MRST
      if (name(nset).eq.'MRST') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'MRSTpdf') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'MRSTgrid') call MRSTpdf(imem)
      if (name(nset).eq.'MRST3grid') call MRSTpdf(imem)
      if (name(nset).eq.'MRST4grid') call MRSTpdf(imem)
#endif
#ifdef MRSTQED
      if (name(nset).eq.'MRST4qed') call MRSTqedpdf(imem)
#endif
#ifdef MRST98
      if (name(nset).eq.'MRST98grid') call MRST98pdf(imem)
#endif
#ifdef MRST06
      if (name(nset).eq.'MRST2006grid') call MRST2006pdf(imem)
#endif
#ifdef CTEQ
      if (name(nset).eq.'EVLCTEQ') then
         call InitEvolvePDF(nset,imem)
         call EVLCTEQpdf(nset)
!          call EVLCTEQpdf(nset,imem)
      endif
      if (name(nset).eq.'CTEQ65grid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CTEQ65cgrid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CTEQ65sgrid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CTEQ6ABgrid') then
         call CTEQ6NewAlpha(nset,imem)
!         call CTEQ6pdf(nset)
      endif
      if (name(nset).eq.'CTEQ66grid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CT10grid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CT12grid') then
         call CTEQ6NewAlpha(nset,imem)
      endif
      if (name(nset).eq.'CTEQ5grid') call CTEQ5pdf(imem)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6pdf(imem)
      if (name(nset).eq.'CTEQ65grid') call CTEQ65pdf(imem)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65pdf(imem)
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65pdf(imem)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65pdf(imem)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6pdf(imem)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65pdf(imem)
      if (name(nset).eq.'CTEQ6LGgrid') call CTEQ6LGpdf(imem)
      if (name(nset).eq.'CT10grid') call CT12pdf(imem)
      if (name(nset).eq.'CT12grid') call CT12pdf(imem)
#endif
#ifdef H1
      if (name(nset).eq.'H12000') then
         call InitEvolvePDF(nset,imem)
         call H1pdf(imem)
      endif         
#endif
#ifdef NNPDF
      if (name(nset).eq.'NNPDF') then
         call InitEvolvePDF(nset,imem)
         call NNPDFpdf(nset)
      endif
      if (name(nset).eq.'NNPDFint') then
         call InitEvolvePDF(nset,imem)
         call NNPDFINTpdf(imem)
      endif
      if (name(nset).eq.'NNPDF20int') then
         call InitEvolvePDF(nset,imem)
         call NNPDFINT20pdf(imem)
      endif
      if (name(nset).eq.'NNPDF20intqed') then
         call InitEvolvePDF(nset,imem)
         call NNPDFINT20qedpdf(imem)
      endif
#endif
#ifdef MSTW
      if (name(nset).eq.'MSTWgrid') call MSTWpdf(imem)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRpdf(imem)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Mpdf(imem)
      if (name(nset).eq.'ABKM09') call ABKM09pdf(imem)
      if (name(nset).eq.'ABM11') call ABM11pdf(imem)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hknpdf(imem)
#endif
#ifdef PHOTONS
!      if (name(nset).eq.'GRV0' .OR. &
!     &    name(nset).eq.'GRV1') call GRVpdf(imem)
      if (name(nset).eq.'SASG') call SASGpdf(imem)
      if (name(nset).eq.'GRVG') call GRVGpdf(imem)
      if (name(nset).eq.'DOG0' .OR. &
     &    name(nset).eq.'DOG1') call DOGpdf(imem)
      if (name(nset).eq.'DGG') call DGGpdf(imem)
      if (name(nset).eq.'LACG') call LACGpdf(imem)
      if (name(nset).eq.'GSG0' .OR. &
     &    name(nset).eq.'GSG1') call GSGpdf(imem)
      if (name(nset).eq.'GSG960' .OR. &
     &    name(nset).eq.'GSG961') call GSG96pdf(imem)
      if (name(nset).eq.'ACFGP') call ACFGPpdf(imem)
      if (name(nset).eq.'WHITG') call WHITGpdf(imem)
#endif
#ifdef PIONS
      if (name(nset).eq.'OWP') call OWPpdf(imem)
      if (name(nset).eq.'SMRSP') call SMRSPpdf(imem)
      if (name(nset).eq.'GRVP0' .OR. &
     &    name(nset).eq.'GRVP1') call GRVPpdf(imem) 
      if (name(nset).eq.'ABFKWP') call ABFKWPpdf(imem)
#endif
#ifdef USER
      if (name(nset).eq.'USER') then
         call InitEvolvePDF(nset,imem)
         call USERpdf(imem)
      endif
      if (name(nset)(1:8).eq.'USERGRID') then
         call InitEvolvePDF(nset,imem)
         call USERGRIDpdf(imem)
      endif
#endif
      return
!
      entry getGridM(nset,nxgrid,nqgrid,gridx,gridq)
#ifdef MRST
      if(name(nset).eq.'MRSTgrid') call MRSTgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if(name(nset).eq.'MRST3grid') call MRSTgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if(name(nset).eq.'MRST4grid') call MRSTgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif 
#ifdef MRSTQED
      if(name(nset).eq.'MRST4qed') call MRSTqedgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif 
#ifdef MRST98
      if(name(nset).eq.'MRST98grid') call MRST98getgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif 
#ifdef MRST06
      if(name(nset).eq.'MRST2006grid') call MRST2006getgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif 
#ifdef MSTW
      if(name(nset).eq.'MSTWgrid') call MSTWgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif 
#ifdef CTEQ
      if (name(nset).eq.'CTEQ6grid') call CTEQ6getgrid(nset,nxgrid,nqgrid,gridx,gridq)  
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6getgrid(nset,nxgrid,nqgrid,gridx,gridq) 
      if (name(nset).eq.'CTEQ5grid') call CTEQ5getgrid(nset,nxgrid,nqgrid,gridx,gridq) 
      if (name(nset).eq.'CTEQ65grid') call CTEQ65getgrid(nset,nxgrid,nqgrid,gridx,gridq) 
      if (name(nset).eq.'CTEQ65cgrid') call CTEQ65getgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'CTEQ65sgrid') call CTEQ65getgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'CTEQ66grid') call CTEQ65getgrid(nset,nxgrid,nqgrid,gridx,gridq) 
      if (name(nset).eq.'CT10grid') call CT12getgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'CT12grid') call CT12getgrid(nset,nxgrid,nqgrid,gridx,gridq) 
#endif
#ifdef NNPDF
!      if (name(nset).eq.'NNPDFint') call NNPDFINTgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'NNPDF20int') call NNPDFINT20getgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'NNPDF20qedint') call NNPDFINT20qedgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
#ifdef HERA
      if (name(nset)(1:8).eq.'HERAGRID') call HERAGRIDgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
#ifdef ALEKHIN
      if (name(nset).eq.'A02M') call A02Mgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'ABKM09') call ABKM09getgrid(nset,nxgrid,nqgrid,gridx,gridq)
      if (name(nset).eq.'ABM11') call ABM11getgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
#ifdef GJR
      if (name(nset)(1:5).eq.'GJR08') call GJRgetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
#ifdef H1
      if (name(nset).eq.'H12000') call H1getgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
#ifdef HKN
      if (name(nset).eq.'HKNgrid') call hkngetgrid(nset,nxgrid,nqgrid,gridx,gridq)
#endif
      return
           
      end
