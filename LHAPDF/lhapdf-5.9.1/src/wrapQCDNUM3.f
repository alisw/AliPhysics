! -*- F90 -*-


      subroutine QCDNUM3evolve(x,Q,pdf) 
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
      call getnset(iset) 
      if (iset.ne.isetlast) then 
        call get_pdfqcd(iset) 
        isetlast = iset 
      endif 
!                                                                       
      Q2=Q*Q 
      nf=3 
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
                                                                        
      pdf(1)=dm+db 
      pdf(2)=um+ub 
      pdf(3)=sb 
      pdf(4)=0d0 
      pdf(5)=0d0 
      pdf(6)=0d0 
      pdf(-1)=db 
      pdf(-2)=ub 
      pdf(-3)=sb 
      pdf(-4)=0d0 
      pdf(-5)=0d0 
      pdf(-6)=0d0 
      return 
!                                                                       
      entry QCDNUM3alfa(alfas,Q) 
      Q2=Q*Q 
      nf=3 
      alfas=QALFAS(Q2,Qlam,nf,iflag) 
      return 
!                                                                       
      entry QCDNUM3read(nset) 
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
      entry QCDNUM3init(nset,Eorder,Q2fit) 
      call qninit 
      call QNISET('ORDER',Eorder+1) 
      call grxdef(nx,xmin) 
      call grqdef(nq,qmin,qmax) 
      call grqinp(Q2fit,1) 
      iq0=IQFROMQ(Q2fit) 
!      nf0= NFLGET(iq0)                                                 
      nf0=3 
      CALL QNBOOK(2,'dm') 
      CALL QNBOOK(3,'um') 
      CALL QNBOOK(4,'dp') 
      CALL QNBOOK(5,'up') 
      CALL QNBOOK(6,'sp') 
      CALL QNLSET('W1ANA',.TRUE.) 
      CALL QNLSET('W2NUM',.TRUE.) 
      CALL QNLSET('W2STF',.FALSE.) 
      call qthres(1d10,2d10) 
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
      entry QCDNUM3pdf(nset) 
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
!                                                                       
!      CALL EVPLUS('DP',iq0,1,nq)                                       
!      CALL EVPLUS('UP',iq0,1,nq)                                       
!      CALL EVPLUS('SP',iq0,1,nq)                                       
      CALL EVOLNP('DP',iq0,1,nq) 
      CALL EVOLNP('UP',iq0,1,nq) 
      CALL EVOLNP('SP',iq0,1,nq) 
!                                                                       
      call getnset(iset) 
      call save_pdfqcd(iset) 
      return 
!                                                                       
      END                                           
