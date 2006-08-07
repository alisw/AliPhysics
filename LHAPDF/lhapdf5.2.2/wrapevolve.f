      subroutine evolvePDF(x,Q,f)
      implicit none
      integer nset,imem,Eorder,IP2
      real*8 x,Q,P2,Q2fit,f(-6:6),alfas,a
      nset = 1
      call evolvePDFM(nset,x,Q,f)
      return
c
      entry evolvePDFp(x,Q,P2,IP2,f)
      nset = 1
      call evolvePDFpM(nset,x,Q,P2,IP2,f)
      return
c
      entry evolvePDFa(x,Q,a,f)
      nset = 1
      call evolvePDFaM(nset,x,Q,a,f)
      return
c
c      entry readevolve
c      nset = 1
c      call readevolveM(nset)
c      return
c
c      entry alfasevolve(alfas,Q)
c      nset = 1
c      call alfasevolveM(nset,alfas,Q)
c      return
c
c      entry initevolution(Eorder,Q2fit)
c      nset = 1
c      call initevolutionM(nset,Eorder,Q2fit)
c      return
c
      entry initPDF(imem)
      nset = 1
      call initPDFM(nset,imem)
      
      return
      end
c
      subroutine evolvePDFaM(nset,x,Q,a,f)
      implicit none
      real*8 x,Q,a,f(-6:6)
      real*8 ruv,rdv,ru,rd,rs,rc,rb,rt,rg
      integer nset
      call eks98(x,q,a,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
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
c      
      subroutine evolvePDFM(nset,x,Q,f)
      implicit none
      include 'parmsetup.inc'
      integer Eorder,index,imem
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem,ip2
      common/NAME/name,nmem,ndef,mem
      integer iset,iimem
      common/SET/iset,iimem
      integer nset,j
      real*8 x,Q,Q2fit,alfas,p2
      real*8 f(-6:6)

      save
*
      call setnset(nset)
c      
c      print *,'this is evolvePDFM, name=',nset,name(nset)
c   set all f's to 0.0d0 at start
c      do j = -6,6
c        f(j) = 0.0d0
c      enddo 
      if (name(nset).eq.'QCDNUM') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3evolve(x,Q,f)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4evolve(x,Q,f)
      if (name(nset)(1:11).eq.'QCDNUM_ZEUS') call ZEUSevolve(x,Q,f)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5evolve(x,Q,f)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6evolve(x,Q,f)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6evolve(x,Q,f)
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQevolve(x,Q,f)
      if (name(nset).eq.'MRSTgrid') call MRSTevolve(x,Q,f)
      if (name(nset).eq.'MRST3grid') call MRSTevolve(x,Q,f)
      if (name(nset).eq.'MRST4grid') call MRSTevolve(x,Q,f)
      if (name(nset).eq.'MRST98grid') call MRST98evolve(x,Q,f)
      if (name(nset).eq.'MRSTpdf') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'MRST') call QCDNUMevolve(x,Q,f)
      if (name(nset).eq.'A02') call A02evolve(x,Q,f)
      if (name(nset).eq.'A02M') call A02Mevolve(x,Q,f)
      if (name(nset).eq.'H12000') call H1evolve(x,Q,f)
      if (name(nset).eq.'GRV') call GRVevolve(x,Q,f)
      if (name(nset).eq.'OWP') call OWPevolve(x,Q,f)
      if (name(nset).eq.'SMRSP') call SMRSPevolve(x,Q,f)
      if (name(nset).eq.'GRVP0') call GRVP0evolve(x,Q,f)
      if (name(nset).eq.'GRVP1') call GRVP1evolve(x,Q,f)
      if (name(nset).eq.'ABFKWP') call ABFKWPevolve(x,Q,f)
      return
*
      entry evolvePDFpM(nset,x,Q,P2,IP2,f)
c
      call setnset(nset)
c      
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
      return
c
      entry readevolve(nset)
*
      read(1,*) name(nset)   
c      print *, 'this is readevolve', name(nset)   
*
      call setnset(nset)
c      
      if (name(nset).eq.'QCDNUM') call QCDNUMread(nset)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMread(nset)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3read(nset)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4read(nset)
      if (name(nset)(1:11).eq.'QCDNUM_ZEUS') call ZEUSread(nset)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5read(nset)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6read(nset)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6read(nset)
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQread(nset)
      if (name(nset).eq.'MRSTgrid') call MRSTread(nset)
      if (name(nset).eq.'MRST3grid') call MRSTread(nset)
      if (name(nset).eq.'MRST4grid') call MRSTread(nset)
      if (name(nset).eq.'MRST98grid') call MRST98read(nset)
      if (name(nset).eq.'MRSTpdf') call QCDNUMread(nset)
      if (name(nset).eq.'MRST') call QCDNUMread(nset)
      if (name(nset).eq.'A02') call A02read(nset)
      if (name(nset).eq.'A02M') call A02Mread(nset)
      if (name(nset).eq.'H12000') call H1read(nset)
      if (name(nset).eq.'GRV') call GRVread(nset)
      if (name(nset).eq.'SASG') call SASGread(nset)
      if (name(nset).eq.'GRVG0' .OR.
     +    name(nset).eq.'GRVG1') call GRVGread(nset)
      if (name(nset).eq.'DOG0' .OR.
     +    name(nset).eq.'DOG1') call DOGread(nset)
      if (name(nset).eq.'DGG') call DGGread(nset)
      if (name(nset).eq.'LACG') call LACGread(nset)
      if (name(nset).eq.'GSG0' .OR. 
     +    name(nset).eq.'GSG1') call GSGread(nset)
      if (name(nset).eq.'GSG960' .OR.
     +    name(nset).eq.'GSG961') call GSG96read(nset)
      if (name(nset).eq.'ACFGP') call ACFGPread(nset)
      if (name(nset).eq.'WHITG') call WHITGread(nset)
      if (name(nset).eq.'OWP') call OWPread(nset)
      if (name(nset).eq.'SMRSP') call SMRSPread(nset)
      if (name(nset).eq.'GRVP0' .OR. 
     +    name(nset).eq.'GRVP1') call GRVPread(nset)
      if (name(nset).eq.'ABFKWP') call ABFKWPread(nset)
      return
*
      entry alfasevolve(nset,alfas,Q)
c
      call setnset(nset)
c      print *,'from alpfasevolveM',nset,Q,name(nset)
c
      if (name(nset).eq.'QCDNUM') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST3') call QCDNUM3alfa(alfas,Q)
      if (name(nset).eq.'QCDNUM_MRST4') call QCDNUM4alfa(alfas,Q)
      if (name(nset)(1:11).eq.'QCDNUM_ZEUS') call ZEUSalfa(alfas,Q)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5alfa(alfas,Q)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6alfa(alfas,Q)
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQalfa(alfas,Q)
      if (name(nset).eq.'MRSTgrid') call MRSTalfa(5,alfas,Q)
      if (name(nset).eq.'MRST3grid') call MRSTalfa(3,alfas,Q)
      if (name(nset).eq.'MRST4grid') call MRSTalfa(4,alfas,Q)
      if (name(nset).eq.'MRST98grid') call MRST98alfa(alfas,Q)
      if (name(nset).eq.'MRSTpdf') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'MRST') call QCDNUMalfa(alfas,Q)
      if (name(nset).eq.'A02') call A02alfa(alfas,Q)
      if (name(nset).eq.'A02M') call A02Malfa(alfas,Q)
      if (name(nset).eq.'H12000') call H1alfa(alfas,Q)
      if (name(nset).eq.'GRV') call GRValfa(alfas,Q)
      if (name(nset).eq.'SASG') call SASGalfa(alfas,Q)
      if (name(nset).eq.'GRVG0' .OR.
     +    name(nset).eq.'GRVG1') call GRVGalfa(alfas,Q)
      if (name(nset).eq.'DOG0' .OR.
     +    name(nset).eq.'DOG1') call DOGalfa(alfas,Q)
      if (name(nset).eq.'DGG') call DGGalfa(alfas,Q)
      if (name(nset).eq.'LACG') call LACGalfa(alfas,Q)
      if (name(nset).eq.'GSG0' .OR.
     +    name(nset).eq.'GSG1') call GSGalfa(alfas,Q)
      if (name(nset).eq.'GSG960' .OR.
     +    name(nset).eq.'GSG961') call GSG96alfa(alfas,Q)
      if (name(nset).eq.'ACFGP') call ACFGPalfa(alfas,Q)
      if (name(nset).eq.'WHITG') call WHITGalfa(alfas,Q)
      if (name(nset).eq.'OWP') call OWPalfa(alfas,Q)
      if (name(nset).eq.'SMRSP') call SMRSPalfa(alfas,Q)
      if (name(nset).eq.'GRVP0' .OR. 
     +    name(nset).eq.'GRVP1') call GRVPalfa(alfas,Q)
      if (name(nset).eq.'ABFKWP') call ABFKWPalfa(alfas,Q)
      return
*
      entry initevolution(nset,Eorder,Q2fit)
c
      call setnset(nset)
c            
      if (name(nset).eq.'QCDNUM') call QCDNUMinit(nset,Eorder,Q2fit)
      if (name(nset)(1:11).eq.'QCDNUM_ZEUS') 
     + call ZEUSinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'CTEQ5grid') call CTEQ5init(Eorder,Q2fit)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6init(Eorder,Q2fit)
      if (name(nset).eq.'EVLCTEQ') call EVLCTEQinit(nset,Eorder,Q2fit)
      if (name(nset).eq.'MRSTgrid') call MRSTinit(Eorder,Q2fit)
      if (name(nset).eq.'MRST3grid') call MRSTinit(Eorder,Q2fit)
      if (name(nset).eq.'MRST4grid') call MRSTinit(Eorder,Q2fit)
      if (name(nset).eq.'MRST98grid') call MRST98init(Eorder,Q2fit)
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
      if (name(nset).eq.'A02') call A02init(Eorder,Q2fit)
      if (name(nset).eq.'A02M') call A02Minit(Eorder,Q2fit)
      if (name(nset).eq.'H12000') call H1init(Eorder,Q2fit)
      if (name(nset).eq.'GRV') call GRVinit(Eorder,Q2fit)
      if (name(nset).eq.'SASG') call SASGinit(Eorder,Q2fit)
      if (name(nset).eq.'GRVG0' .OR.
     +    name(nset).eq.'GRVG1') call GRVGinit(Eorder,Q2fit)
      if (name(nset).eq.'DOG0' .OR.
     +    name(nset).eq.'DOG1') call DOGinit(Eorder,Q2fit)
      if (name(nset).eq.'DGG') call DGGinit(Eorder,Q2fit)
      if (name(nset).eq.'LACG') call LACGinit(Eorder,Q2fit)
      if (name(nset).eq.'GSG0' .OR.
     +    name(nset).eq.'GSG1') call GSGinit(Eorder,Q2fit)
      if (name(nset).eq.'GSG960' .OR.
     +    name(nset).eq.'GSG961') call GSG96init(Eorder,Q2fit)
      if (name(nset).eq.'ACFGP') call ACFGPinit(Eorder,Q2fit)
      if (name(nset).eq.'WHITG') call WHITGinit(Eorder,Q2fit)
      if (name(nset).eq.'OWP') call OWPinit(Eorder,Q2fit)
      if (name(nset).eq.'SMRSP') call SMRSPinit(Eorder,Q2fit)
      if (name(nset).eq.'GRVP0' .OR. 
     +    name(nset).eq.'GRVP1') call GRVPinit(Eorder,Q2fit)
      if (name(nset).eq.'ABFKWP') call ABFKWPinit(Eorder,Q2fit)
      return
*
      entry initPDFM(nset,imem)
c
      call setnset(nset)
      call setnmem(nset,imem)
c            
c      print *,'entered initPDFM,',nset,imem,name(nset)
      iimem = imem
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
      if (name(nset)(1:11).eq.'QCDNUM_ZEUS') then
         call InitEvolvePDF(nset,imem)
         call ZEUSpdf(nset)
      endif
      if (name(nset).eq.'MRST') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'MRSTpdf') then
         call InitEvolvePDF(nset,imem)
         call QCDNUMpdf(nset)
      endif
      if (name(nset).eq.'EVLCTEQ') then
         call InitEvolvePDF(nset,imem)
         call EVLCTEQpdf(nset)
c          call EVLCTEQpdf(nset,imem)
      endif
      if (name(nset).eq.'CTEQ6ABgrid') then
         call CTEQ6NewAlpha(nset,imem)
c         call CTEQ6pdf(nset)
      endif
      if (name(nset).eq.'H12000') then
         call InitEvolvePDF(nset,imem)
         call H1pdf(imem)
      endif
      if (name(nset).eq.'CTEQ5grid') call CTEQ5pdf(imem)
      if (name(nset).eq.'CTEQ6grid') call CTEQ6pdf(imem)
      if (name(nset).eq.'CTEQ6ABgrid') call CTEQ6pdf(imem)
      if (name(nset).eq.'MRSTgrid') call MRSTpdf(imem)
      if (name(nset).eq.'MRST3grid') call MRSTpdf(imem)
      if (name(nset).eq.'MRST4grid') call MRSTpdf(imem)
      if (name(nset).eq.'MRST98grid') call MRST98pdf(imem)
      if (name(nset).eq.'A02') call A02pdf(imem)
      if (name(nset).eq.'A02M') call A02Mpdf(imem)
      if (name(nset).eq.'GRV0' .OR.
     +    name(nset).eq.'GRV1') call GRVpdf(imem)
      if (name(nset).eq.'SASG') call SASGpdf(imem)
      if (name(nset).eq.'GRVG') call GRVGpdf(imem)
      if (name(nset).eq.'DOG0' .OR.
     +    name(nset).eq.'DOG1') call DOGpdf(imem)
      if (name(nset).eq.'DGG') call DGGpdf(imem)
      if (name(nset).eq.'LACG') call LACGpdf(imem)
      if (name(nset).eq.'GSG0' .OR.
     +    name(nset).eq.'GSG1') call GSGpdf(imem)
      if (name(nset).eq.'GSG960' .OR.
     +    name(nset).eq.'GSG961') call GSG96pdf(imem)
      if (name(nset).eq.'ACFGP') call ACFGPpdf(imem)
      if (name(nset).eq.'WHITG') call WHITGpdf(imem)
      if (name(nset).eq.'OWP') call OWPpdf(imem)
      if (name(nset).eq.'SMRSP') call SMRSPpdf(imem)
      if (name(nset).eq.'GRVP0' .OR.
     +    name(nset).eq.'GRVP1') call GRVPpdf(imem) 
      if (name(nset).eq.'ABFKWP') call ABFKWPpdf(imem)
      return
*
      end
