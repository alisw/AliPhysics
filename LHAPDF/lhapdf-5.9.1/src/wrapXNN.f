!***********************************************************************
!     Initial parameterization: NN (simplified wrt our code)            
!***********************************************************************
                                                                        
      subroutine LH_PDFIN(X,PDFTMP) 
      implicit none 
!                                                                       
      INCLUDE 'parmsetup.inc' 
      CHARACTER*16 name(nmxset) 
      INTEGER nmem(nmxset),ndef(nmxset),mmem 
      COMMON/NAME/name,nmem,ndef,mmem 
!                                                                       
      INTEGER KREP 
      integer MXREP 
      parameter(MXREP=1e3) 
      integer MXPDF 
      parameter(MXPDF=13) 
                                !# pdfs parametrized with NN            
      integer NTOTPDF 
      parameter(NTOTPDF=5) 
      integer MXPAR 
      parameter(MXPAR=2e2) 
!                                                                       
      integer JPDF 
      REAL*8 X 
      REAL*8 PDFTMP(MXPDF),PDFFLV(NTOTPDF) 
!                                                                       
!     Warnings                                                          
                                                                        
      KREP = MMEM 
      if(X.gt.1d0) then 
         write(6,*) "X = ",X 
         write(6,*) 'X greater than 1 in pdfin routine !' 
         call exit(-10) 
      elseif(X.le.0d0) then 
         write(6,*) "X = ",X 
         write(6,*) 'X less or equal 0 in pdfin routine !' 
         call exit(-10) 
      endif 
!                                                                       
!     Neural network parton distributions                               
      do JPDF=1,NTOTPDF 
         call LH_NNPDF(X,JPDF,KREP,PDFFLV(JPDF)) 
      enddo 
      call LH_PDFINPAR2EVLN(PDFFLV,PDFTMP) 
!                                                                       
      return 
      END                                           
                                                                        
!***********************************************************************
!     those commons might be changed depending on those put in inputPDF.
!***********************************************************************
                                                                        
      subroutine LH_NNPDF(X,JPDF,KREP,PDFNET) 
      implicit none 
!                                                                       
      integer KREP 
      integer MXREP 
      parameter(MXREP=1e3) 
      integer MXPDF 
      parameter(MXPDF=13) 
                                !# pdfs parametrized with NN            
      integer NTOTPDF 
      parameter(NTOTPDF=5) 
      integer MXPAR 
      parameter(MXPAR=2e2) 
!                                                                       
                                !NN basic parameters                    
      integer MXL,MXN 
      parameter(MXL=5,MXN=10) 
!                                                                       
                                            ! NN architecture           
      integer TNL(NTOTPDF),NEU(MXL,NTOTPDF) 
      REAL*8 XNOR(2,NTOTPDF),PDFNOR(2,NTOTPDF) 
      REAL*8 PDFEXP(2,NTOTPDF) 
      REAL*8 PDFDELTA(NTOTPDF) 
      REAL*8 XDELTALIN(NTOTPDF),XDELTALOG(NTOTPDF) 
      common/nnpdf10CNNARC/PDFEXP,PDFNOR,PDFDELTA,XDELTALIN             &
     &     ,XDELTALOG,XNOR,NEU,TNL                                      
!                                                                       
      real*8 PARTR(MXPAR,0:MXREP,MXPDF) 
      real*8 ANORMTR(MXPDF,0:MXREP) 
      common/nnpdf10CPARTR/PARTR,ANORMTR 
!                                                                       
      INTEGER ipdf,JPDF,NTMP(MXL),L 
      REAL*8 PDFOUT,X 
      REAL*8 XNNIN(MXN),XNNOUT,XNN(MXN,MXL) 
      REAL*8 PDFNET 
!                                                                       
!     Temporary array of neuron number                                  
      do L=1,TNL(JPDF) 
         NTMP(L)=NEU(L,JPDF) 
      enddo 
!                                                                       
      do IPDF=1,NTOTPDF 
!     Evaluate fixed coefficients used in the normalizaion of NN        
         XDELTALIN(IPDF)=1.d0/(XNOR(2,IPDF)-XNOR(1,IPDF))*0.8d0 
         XDELTALOG(IPDF)=1.d0/dlog(XNOR(2,IPDF)/XNOR(1,IPDF))*0.8d0 
         PDFDELTA(IPDF)=(PDFNOR(2,IPDF)-PDFNOR(1,IPDF))/0.8d0 
      enddo 
!                                                                       
      XNNIN(1)=(X-XNOR(1,JPDF))*XDELTALIN(JPDF)+0.1d0 
      XNNIN(2)=dlog(X/XNOR(1,JPDF))*XDELTALOG(JPDF)+0.1d0 
!                                                                       
!     Call neural network routine                                       
!                                                                       
      call LH_NNOUT(XNNIN,TNL(JPDF),NTMP,PARTR(1,KREP,JPDF),XNN,XNNOUT) 
!                                                                       
!     Output normalization: the delta is defined in pdfcreate.f         
!                                                                       
      PDFOUT=PDFNOR(1,JPDF)+(XNNOUT - 0.1d0)*PDFDELTA(JPDF) 
!                                                                       
!     Preprocessing exponents                                           
!                                                                       
      PDFNET=((1.-x)**PDFEXP(1,JPDF))*PDFOUT*x**(-PDFEXP(2,JPDF)) 
!                                                                       
!     PDF normalization factor (fixed in some cases by sum rules)       
!                                                                       
      PDFNET= ANORMTR(JPDF,KREP) * PDFNET 
!                                                                       
      return 
      END                                           
!                                                                       
!     file: NNSUB.f --> It evaluates NN output                          
!                                                                       
      subroutine LH_NNOUT(xnnin,tnl,n,PAR,xnn,xnnout) 
      implicit none 
!                                                                       
      integer MXL,MXN 
      parameter(MXL=5,MXN=10) 
      integer MXPAR 
      parameter(MXPAR=2e2) 
!                                                                       
      integer tnl,n,l,i,ii,j 
      REAL*8 xnnin,xnnout,xnn,LH_G,h 
      dimension n(mxl) 
      dimension xnnin(mxn) 
      dimension xnn(mxn,mxl) 
!                                                                       
      REAL*8 PAR(mxpar) 
      integer IPAR 
!                                                                       
      IPAR=0 
      do ii=1,n(1) 
         xnn(ii,1)=xnnin(ii) 
      enddo 
      do l=1,tnl-1 
         xnn(n(l)+1,l)=1d0 
      enddo 
                                                                        
      do   l=2,tnl-1 
         do  i=1,n(l) 
            h=0. 
            do j=1,n(l-1)+1 
               IPAR=IPAR+1 
               h=h+PAR(IPAR)*xnn(j,l-1) 
            enddo 
            xnn(i,l)=LH_G(h) 
         enddo 
      enddo 
      do i=1,n(tnl) 
         h=0. 
         do j=1,n(l-1)+1 
            IPAR=IPAR+1 
            h=h+PAR(IPAR)*xnn(j,l-1) 
         enddo 
         xnn(i,tnl)=h 
      enddo 
                                                                        
      xnnout=xnn(1,tnl) 
                                                                        
      return 
      END                                           
!                                                                       
!     Neural network activation function                                
!                                                                       
      function LH_g(x) 
      implicit none 
      REAL*8 LH_g,x 
!                                                                       
      LH_g=1./(1.+DExp(-x)) 
!                                                                       
      return 
      END                                           
