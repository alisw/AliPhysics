! -*- F90 -*-


      subroutine QCDNUM4evolve(x,Q,pdf) 
      implicit none 
      character*64 gridname 
      character*16 s1,dummy 
      integer Eorder,index 
      real*8 x,Q,Qlam,Q2,Q2fit,alfas0,scale0,alfas 
      real*8 mc,mc2,mb,mb2,mt,mt2,tc,tc2,tb,tb2 
      real*8 singlet,dm,um,dp,up,sp,ub,db,sb,cb,bb,bp,cp,gl,xx 
      real*8 QALFAS,QPDFXQ,XFROMIX 
      real*8 f(-6:6),pdf(-6:6) 
      integer iq0,iqc,iqb,nf,nf0,qnerr,NFLGET,iflag,ix,i,IQFROMQ 
      integer nx,nq 
      integer nset,iset,isetlast 
      data isetlast/-1/ 
      real*8 xmin,xmax,qmin,qmax,S 
      save iq0,iqc,iqb,nf0,mc2,mb2,tc2,tb2,mt2 
      save nx,xmin,xmax,nq,qmin,qmax,gridname 
!                                                                       
!      print *,'QCDNUM4evolve'                                          
      call getnset(iset) 
      if (iset.ne.isetlast) then 
        call get_pdfqcd(iset) 
        isetlast = iset 
      endif 
!                                                                       
      Q2=Q*Q 
      nf=4 
      if (Q2.lt.tc2) nf=3 
      pdf(0)=QPDFXQ('GLUON'  ,x,Q2,iflag) 
      singlet=QPDFXQ('SINGLET',x,Q2,iflag) 
      dm= QPDFXQ('DM',x,Q2,IFLAG) 
      um= QPDFXQ('UM',x,Q2,IFLAG) 
      dp= QPDFXQ('DP',x,Q2,IFLAG) 
      up= QPDFXQ('UP',x,Q2,IFLAG) 
      sp= QPDFXQ('SP',x,Q2,IFLAG) 
      ub=0.5d0*(up-um+singlet/dble(nf)) 
      db=0.5d0*(dp-dm+singlet/dble(nf)) 
      sb=0.5d0*(sp+singlet/dble(nf)) 
      cb=0.d0 
      if (nf.ge.4) then 
         cp= QPDFXQ('CP',X,Q2,IFLAG) 
         cb=0.5d0*(cp+singlet/dble(nf)) 
      end if 
                                                                        
      pdf(1)=dm+db 
      pdf(2)=um+ub 
      pdf(3)=sb 
      pdf(4)=cb 
      pdf(5)=0d0 
      pdf(6)=0d0 
      pdf(-1)=db 
      pdf(-2)=ub 
      pdf(-3)=sb 
      pdf(-4)=cb 
      pdf(-5)=0d0 
      pdf(-6)=0d0 
      return 
!                                                                       
      entry QCDNUM4alfa(alfas,Q) 
      Q2=Q*Q 
      nf=4 
      if (Q2.lt.mc2) nf=3 
      alfas=QALFAS(Q2,Qlam,nf,iflag) 
      return 
!                                                                       
      entry QCDNUM4read(nset) 
!      read(1,*) s1                                                     
!      print *,s1                                                       
!         gridname='large.grid'                                         
!         nx=400                                                        
!         xmin=1d-6                                                     
!         xmax=1d0                                                      
!         nq=112                                                        
!         qmin=1d0                                                      
!         qmax=1d10                                                     
!         read(1,*) dummy                                               
!      else                                                             
         read(1,*) gridname,nx,xmin,xmax,nq,qmin,qmax 
!*      endif                                                           
!      iset = nset                                                      
      return 
!                                                                       
      entry QCDNUM4init(nset,Eorder,Q2fit) 
      call qninit 
      call QNISET('ORDER',Eorder+1) 
      call grxdef(nx,xmin) 
      call grqdef(nq,qmin,qmax) 
      call grqinp(Q2fit,1) 
      call getQMassM(nset,4,mc) 
      mc2=mc*mc 
!      print *,mc                                                       
      CALL QNRSET('CMASS',mc) 
      CALL QNRSET('MCALF',mc) 
      call getThresholdM(nset,4,tc) 
      tc2=tc*tc 
      call GRQINP(tc2,1) 
      CALL GRQINP(tc2-1.0d-4,1) 
!      print *,'setting iq0',Q2fit                                      
      call qthres(tc2,2d10) 
      iq0=IQFROMQ(Q2fit) 
      iqc=IQFROMQ(tc2) 
      nf0= NFLGET(iq0) 
!      print *,iq0,iqc,nf0,q2fit,tc2                                    
!      nf0=3                                                            
!      print *,nf0                                                      
      CALL QNBOOK(2,'dm') 
      CALL QNBOOK(3,'um') 
      CALL QNBOOK(4,'dp') 
      CALL QNBOOK(5,'up') 
      CALL QNBOOK(6,'sp') 
      CALL QNBOOK(7,'cp') 
      CALL QNLSET('W1ANA',.TRUE.) 
      CALL QNLSET('W2NUM',.TRUE.) 
      CALL QNLSET('W2STF',.FALSE.) 
      if (index(gridname,'none').eq.1) then 
         call qnfilw(0,0) 
      else 
         qnerr=-1 
         open(unit=2,status='old',file=gridname,                        &
     &        form='unformatted',err=1)                                 
         call QNREAD(2,1,qnerr) 
    1    close(2) 
         if (qnerr.ne.0) then 
            write(*,*) 'Grid file problem: ',gridname 
            if (qnerr.lt.0) then 
               write(*,*) 'Grid file does not exist' 
               write(*,*) 'Calculating and creating grid file' 
               call qnfilw(0,0) 
               open(unit=2,status='unknown',file=gridname,              &
     &              form='unformatted')                                 
               call QNDUMP(2) 
               close(2) 
            else 
               write(*,*) 'Existing grid file is inconsistent' 
               if (qnerr.eq.1)                                          &
     &              write(*,*) 'Defined grid different'                 
               if (qnerr.eq.2)                                          &
     &              write(*,*) 'Heavy quark weight table different'     
               if (qnerr.eq.3)                                          &
     &              write(*,*) 'Charm mass different'                   
               if (qnerr.eq.4)                                          &
     &              write(*,*) 'Bottom mass different'                  
               stop 
            endif 
         endif 
      endif 
      return 
!                                                                       
      entry QCDNUM4pdf(nset) 
!      print *,'entering QCDNUMpdf',nset                                
      call GetAlfas(nset,alfas0,scale0) 
!      print *,alfas0,scale0                                            
      Q2=scale0*scale0 
      CALL QNRSET('ALFAS',alfas0) 
      CALL QNRSET('ALFQ0',Q2) 
      DO ix = 1,nx 
         xx = XFROMIX(ix) 
!         print *,'calling parmPDF',ix                                  
         call parmPDF(nset,xx,f) 
!         if(ix.lt.6) print *,nset,xx,f                                 
         singlet=0.d0 
         do i=1,nf0 
            singlet=singlet+f(i)+f(-i) 
!            print  *,i,singlet                                         
         end do 
         gl=f(0) 
         dm=f(1)-f(-1) 
         um=f(2)-f(-2) 
         dp=f(1)+f(-1)-singlet/dble(nf0) 
         up=f(2)+f(-2)-singlet/dble(nf0) 
         sp=f(3)+f(-3)-singlet/dble(nf0) 
!         print *,ix,iq0                                                
         CALL QNPSET('SINGLET',ix,iq0,singlet) 
         CALL QNPSET('GLUON',ix,iq0,gl) 
         CALL QNPSET('DM',ix,iq0,DM) 
         CALL QNPSET('UM',ix,iq0,UM) 
         CALL QNPSET('DP',ix,iq0,DP) 
         CALL QNPSET('UP',ix,iq0,UP) 
         CALL QNPSET('SP',ix,iq0,SP) 
      ENDDO 
      CALL EVOLSG(iq0,1,nq) 
      CALL EVOLNM('DM',iq0,1,nq) 
      CALL EVOLNM('UM',iq0,1,nq) 
      CALL EVLSEA4('dp',iq0,iqc,nq) 
      CALL EVLSEA4('up',iq0,iqc,nq) 
      CALL EVLSEA4('sp',iq0,iqc,nq) 
                                                                        
!      print *,'calling evols - heavy...'                               
                                                                        
!--   Heavy quark evolution                                             
                                                                        
      CALL EVOLCP4('cp',iqc,nq) 
!                                                                       
      call getnset(iset) 
      call save_pdfqcd(iset) 
      return 
!                                                                       
      END                                           
      subroutine EVLSEA4(name,IQ0,IQC,NQGRI) 
      implicit none 
                                                                        
      CHARACTER*(*) name 
      integer iq0,iqc,nqgri 
      real*8 f34,f45,f43,f54,factor 
      parameter(f34=1.d0/12.d0,f43=-1.d0/12.d0) 
                                                                        
!      print *,iq0,iqc,nqgri                                            
      If(IQ0.le.IQC)then 
         CAll EVPLUS(name,IQ0,1,IQC) 
         CALL QADDSI(name,IQC,f34) 
         CAll EVPLUS(name,IQC,IQC,NQGRI) 
      else if(IQ0.gt.IQC)then 
         CAll EVPLUS(name,IQ0,IQC,NQGRI) 
         CALL QADDSI(name,IQC,f43) 
         CAll EVPLUS(name,IQC,1,IQC) 
      end if 
                                                                        
      END                                           
                                                                        
                                                                        
      subroutine EVOLCP4(name,IQC,NQGRI) 
      implicit none 
                                                                        
      CHARACTER*(*) name 
      integer iq0,iqc,nqgri 
      real*8 f4 
      parameter(f4=-1.d0/4.d0) 
                                                                        
!                                                                       
!     First set to zero to avoid adding -1/4Singl at each iteration     
!                                                                       
                                                                        
      CAll QNPNUL(name) 
      CALL QADDSI(name,IQC,f4) 
      CAll EVPLUS(name,IQC,IQC,NQGRI) 
                                                                        
      END                                           
