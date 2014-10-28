!***********************************************                        
!                                                                       
!                                                                       
!     wrapNNPDFgrid-lite.f 
!     Special low-memory version with only single fit ....                                                  
!     Routine called by LHAPDF package for calculating                  
!     the value of all xpdfs at x and Q from replica KREP               
!     in a (x,Q) point as called by NNPDF.LHgrid file.                  
!                                                                       
!     In 'wrapevolve.f' the package calls:                              
!     IF(NAME(NSET).EQ.'NNPDFint') call NNPDFINTevolve(x,Q,f)           
!     IF(NAME(NSET).EQ.'NNPDFint') call NNPDFINTread(nset)              
!     IF(NAME(NSET).EQ.'NNPDFint') call NNPDFINTalfa(alfas,Q)           
!     IF(NAME(NSET).EQ.'NNPDFint') call NNPDFINTinit(nset,Eorder,Q2fit)*
!                                                                       
!***********************************************                        
                                                                        
                                                                        
      subroutine NNPDFINTevolve(X,Q,XPDF) 
      IMPLICIT none 
!                                                                       
      include 'parmsetup.inc' 
      character*16 name(nmxset) 
      character*80 line 
      character*512 setpath 
      integer nmem(nmxset),ndef(nmxset),mmem,mem 
      common/NAME/name,nmem,ndef,mmem 
      integer nset,iset,pdfmem 
      real*8 parm(nopmax) 
      real*8 tmp(13)
!                                                                       
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
!                                                                       
      INTEGER order 
      REAL*8 alfas,alphaNNPDF 
      REAL*8 Eorder,Q2fit,mass 
!                                                                       
      INTEGER MXREP 
      PARAMETER(MXREP=nmxset-1) 
      INTEGER NREP 
      common/nnpdf10CNREP/NREP 
!                                                                       
      INTEGER NX,NQ2,NPL 
      PARAMETER(NX=60,NQ2=50) 
      PARAMETER(NPL=3000) 
      INTEGER NXX,NQQ2 
      REAL*8 Q2MIN,Q2MAX,XPDFMIN,XPDFMAX 
      REAL*8 XG(NX,nmxset),Q2G(NQ2,nmxset),XPDFEV(NX,NQ2,-6:6,0:MXREP) 
      INTEGER IX,IQ2 
      common/nnpdf10CPDFGR/XPDFEV,XG,Q2G,IX,IQ2 
!                                                                       
      INTEGER ipt,imodev,ivfn,itmc 
      COMMON/NNPDF10EVFLAGS/ipt,imodev,ivfn,itmc 
      REAL*8 q0,alfas0 
      REAL*8 q20,qth(4:6) 
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
      COMMON/nnpdf10vfns/q2th,asref,q2ref 
!                                                                       
      INTEGER I,J,K 
      INTEGER IPDF,KREP,LH_JISEARCH,IINTERP 
      INTEGER IDUM,JDUM 
      REAL*8 X,Q,QQ,Q2,QQ2,XPDF(-6:6) 
      REAL*8 AXB(NX,NQ2,-6:6), BXB(NX,NQ2,-6:6) 
      REAL*8 CXB(NX,NQ2,-6:6),TQ,DX 
      REAL*8 XPDF1(-6:6),XPDF2(-6:6) 
      REAL*8 XCH 
      PARAMETER(XCH=1D-1) 
!                                                                       
      integer m,n,nmax,mmax 
                                ! order of pol. interpolation           
      parameter(m=4,n=4) 
      parameter(nmax=1e3,mmax=1e3) 
      double precision dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax) 
      integer ix1a(m),ix2a(n) 
!                                                                       
!
!      KREP = mmem 
! always use the 0 element in this lite version
      KREP = 0 
!     Set correct scale                                                 
      Q2=Q**2d0 
                                      
      call getnset(iset)
                                                                        
!     Check kinematic point is within allowed range                     
                                                                        
      call GetXminM(iset,KREP,XPDFMIN) 
      call GetXmaxM(iset,KREP,XPDFMAX) 
      call GetQ2maxM(iset,KREP,Q2MAX) 
      call GetQ2minM(iset,KREP,Q2MIN) 
!                                                                       
      IF ( X.LT.XPDFMIN .OR. X.GT.XPDFMAX ) THEN 
         WRITE(6,2000) 
 2000    FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE -- STOP') 
         write(6,*) "X= ",X," XMAX, XMIN = ",XPDFMAX,XPDFMIN 
      ENDIF 
!                                                                       
      IF ( Q2.LT.Q2MIN .OR. Q2.GT.Q2MAX ) THEN 
         WRITE(6,2001) 
 2001    FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE -- STOP') 
         write(6,*) "Q2 ,Q2MIN, Q2MAX = ",Q2,Q2MIN,Q2MAX 
      ENDIF 
                                                                        
!     Select higher-order polynomial interpolation                      
                                                                        
      IINTERP = 1 
!                                                                       
      if(IINTERP.eq.0) then 
!                                                                       
!     CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION         
!                                                                       
!     spline coefficients                                               
         DO IPDF = -6,6,1 
            DO IQ2 = 1,NQ2 
!               CALL LH_JSPLINE(AXB,BXB,CXB,IPDF,IQ2,KREP) 
               CALL LH_JSPLINE(AXB,BXB,CXB,IPDF,IQ2,KREP+iset-1) 
            ENDDO 
         ENDDO 
                                                                        
!     Binary search of points in grid                                   
         IQ2 = LH_JISEARCH(NQ2,Q2G(1,iset),Q2) 
         IF (IQ2 .EQ. NQ2) IQ2 = NQ2-1 
         IX = LH_JISEARCH(NX,XG(1,iset),X) 
         DX = X - XG(IX,iset) 
                                                                        
!     Compute the values of xpdfs for two neighbouring values of Q2     
!     using splines of order 3                                          
                                                                        
         DO IPDF = -6,6,1 
            XPDF1(IPDF) = XPDFEV(IX,IQ2,IPDF,KREP+iset-1)                      &
     &           + DX*(AXB(IX,IQ2,IPDF) + DX*(BXB(IX,IQ2,IPDF)          &
     &           + DX*CXB(IX,IQ2,IPDF)) )                               
            XPDF2(IPDF) = XPDFEV(IX,IQ2+1,IPDF,KREP+iset-1)                    &
     &           + DX*(AXB(IX,IQ2+1,IPDF) + DX*(BXB(IX,IQ2+1,IPDF)      &
     &        + DX*CXB(IX,IQ2+1,IPDF)) )                                
         ENDDO 
                                                                        
!     Linear interpolation in log Q2                                    
                                                                        
         TQ = (DLOG(Q2)-DLOG(Q2G(IQ2,iset)))                                 &
     &        / (DLOG(Q2G(IQ2+1,iset))-DLOG(Q2G(IQ2,iset)))                       
!                                                                       
         DO IPDF = -6,6,1 
            XPDF(IPDF)  = (1.0D0-TQ)*XPDF1(IPDF) + TQ*XPDF2(IPDF) 
         ENDDO 
                                                                        
                                                                        
      elseif(IINTERP.eq.1) then 
!                                                                       
         IQ2 = LH_JISEARCH(NQ2,Q2G(1,iset),Q2) 
         IF (IQ2 .EQ. NQ2) IQ2 = NQ2-1 
         IX = LH_JISEARCH(NX,XG(1,iset),X) 
                                                                        
!     Assign grid for interpolation. M, N -> order of polyN interpolatio
                                                                        
         do I=1,M 
            if(IX.ge.M/2.and.IX.le.(NX-M/2)) IX1A(I) = IX - M/2 + I 
            if(IX.lt.M/2) IX1A(I) = I 
            if(IX.gt.(NX-M/2)) IX1A(I) = (NX - M) + I 
                                                                        
!     Check grids                                                       
            if(IX1A(I).le.0.or.IX1A(I).gt.NX) then 
               write(6,*) "Error in grids" 
               write(6,*) "I, IXIA(I) = ",I, IX1A(I) 
               call exit(-10) 
            endif 
         enddo 
                                                                        
         do J=1,N 
            if(IQ2.ge.N/2.and.IQ2.le.(NQ2-N/2)) IX2A(J) = IQ2 - N/2 + J 
            if(IQ2.lt.N/2) IX2A(J) = J 
            if(IQ2.gt.(NQ2-N/2)) IX2A(J) = (NQ2 - N) + J 
!     Check grids                                                       
            if(IX2A(J).le.0.or.IX2A(J).gt.NQ2) then 
               write(6,*) "Error in grids" 
               write(6,*) "J, IXIA(J) = ",J,IX2A(J) 
               call exit(-10) 
            endif 
         enddo 
                                                                        
!     uncomment for 3rd order interp.                                   
                                                                        
!         IQ2 = LH_JISEARCH(NQ2,Q2G,Q2)                                 
!         IX  = LH_JISEARCH(NX,XG,X)                                    
!         IF((IX+(M-1)/2).GT.NX) IX = NX - (M-1)/2                      
!         IF((IX-(M-1)/2).LT.1) IX = (M+1)/2                            
!         IDUM = 0                                                      
!         DO I = -(M-1)/2,(M-1)/2,1                                     
!            IDUM = IDUM +1                                             
!            IX1A(IDUM) = IX + I                                        
!         ENDDO                                                         
                                                                        
!         IF((IQ2+(N-1)/2).GT.NQ2) IQ2 = NQ2 - (N-1)/2                  
!         IF((IQ2-(N-1)/2).LT.1) IQ2 = (N+1)/2                          
                                                                        
!         JDUM = 0                                                      
!         DO J = -(N-1)/2,(N-1)/2,1                                     
!            JDUM = JDUM +1                                             
!            IX2A(JDUM) = IQ2 + J                                       
!         ENDDO                                                         
!                                                                       
!     Define points where to evaluate interpolation                     
!     Choose between linear or logarithmic (x,Q2) interpolation         
                                                                        
         IF(X.LT.XCH)THEN 
            X1=dlog(X) 
         ELSE 
            X1=X 
         ENDIF 
         X2=dlog(Q2) 
                                                                        
         DO IPDF = -6,6,1 
                                                                        
!     Choose between linear or logarithmic (x,Q2) interpolation         
                                                                        
            DO I=1,M 
               IF(X.LT.XCH)THEN 
                  X1A(I)= dlog(XG(IX1A(I),iset)) 
               ELSE 
                  X1A(I)= XG(IX1A(I),iset) 
               ENDIF 
               DO J=1,N 
                  X2A(J) = dlog(Q2G(IX2A(J),iset)) 
                  YA(I,J) = XPDFEV(IX1A(I),IX2A(J),IPDF,KREP+iset-1) 
               enddo 
            enddo 
                                                                        
!     2D polynomial interpolation                                       
            call lh_polin2(x1a,x2a,ya,m,n,x1,x2,y,dy) 
            XPDF(IPDF) = y 
                                                                        
         enddo 
      endif 
!                                                                       
      RETURN 
                                                                        
!********************************************************               
      entry NNPDFINTgetgrid(nset,ngridx,ngridq,gridx,gridq)
      do jx=1,nx
          gridx(jx)=xg(jx,nset)
      enddo
      do jq=1,nq2
          gridq(jq)=q2g(jq,nset)
      enddo
      ngridx=nx
      ngridq=nq2        
      return
!********************************************************               
                                                                        
      ENTRY NNPDFINTread(nset) 
!                                                                       
      READ(1,*)nmem(nset),ndef(nset) 
                                                                        
!     Set number of members                                             
      call setnmem(nset,nmem) 
                                                                        
!     Read the grid in x                                                
      READ(1,*) nxx 
      IF(NXX.NE.NX)WRITE(*,*)"WARNING CHANGE NX ACCORDING TO .LHgrid"
      DO ix = 1,nxx 
         READ(1,*) xg(ix,nset) 
      ENDDO 
                                                                        
!     Read the grid in Q2                                               
      READ(1,*) nqq2 
      IF(NQQ2.NE.NQ2)WRITE(*,*)"WARNING CHANGE NQ2 ACCORDING TO .LHgrid"
      READ(1,*) q2min 
      DO iq2 = 1,nqq2 
         READ(1,*) q2g(iq2,nset) 
      ENDDO 
                                                                        
!     Read the number of replicas                                       
      READ(1,*) NREP 
                                                                        
! - dummy read in to get to End: (stream 1 is still open)               
      DO K=0,NREP 
         DO IX=1,NX 
            DO IQ2=1,NQ2 
               READ(1,'(a)') line 
            ENDDO 
         ENDDO 
      ENDDO 
!                                                                       
      RETURN 
                                                                        
!********************************************************               
                                                                        
      ENTRY NNPDFINTalfa(alfas,Q) 
!                                                                       
      QQ = Q 
      alfas =  alphaNNPDF(QQ) 
!                                                                       
      RETURN 
                                                                        
!********************************************************               
                                                                        
      ENTRY NNPDFINTinit(nset,Eorder,Q2fit) 
!                                                                       
      IMODEV = 0 
      IVFN = 1 
      ITMC = 1 
!                                                                       
      CALL GetOrderPDFM(nset,order) 
      IPT = order 
!                                                                       
      CALL GetQ2fitM(nset,QQ2) 
      Q2fit = QQ2 
      Q20   = QQ2 
!                                                                       
      call GetQmassM(nset,4,mass) 
      QTH(4) = mass 
      call GetQmassM(nset,5,mass) 
      QTH(5) = mass 
      call GetQmassM(nset,6,mass) 
      QTH(6) = mass 
!                                                                       
      DO i = 4,6 
         q2th(i) = qth(i)**2d0 
      ENDDO 
!                                                                       
                                    ! added for filling Fparm->asref    
      call initEVOLVEpdf(nset,mmem) 
      CALL GetAlfas(nset,alfas0,Q0) 
      asref = alfas0 
      q2ref = q0**2d0 
!                                                                       
      RETURN 
                                                                        
!********************************************************               
!                                                                       
      entry NNPDFINTpdf(mem)
! have to reopen stream 1
        call getnset(iset)                                               
        call setnmem(iset,mem) 
   
        call getsetpath(setpath) 
        open(1,file=setpath(1:len_trim(setpath)),action='READ') 
        line = '' 
	
        do while (line(3:12).ne.'Evolution:') 
           read(1,'(a)')line
        enddo 
        read(1,'(a)')line 
        read(1,'(a)')line 
!                                                                       
      READ(1,*)nmem(iset),ndef(iset) 
      pdfmem = mem 
      IF (pdfmem.LT.0) THEN 
         WRITE(*,*) 'NNPDF set:' 
         WRITE(*,*) 'PDF member out of range:' 
         WRITE(*,*) 'member = ',pdfmem 
         STOP 
      ENDIF 
                                                                        
!     Read the grid in x                                                
      READ(1,*) nxx 
      IF(NXX.NE.NX)WRITE(*,*)"WARNING CHANGE NX ACCORDING TO .LHgrid"
      DO ix = 1,nxx 
         READ(1,*) xg(ix,iset) 
      ENDDO 
                                                                        
!     Read the grid in Q2                                               
      READ(1,*) nqq2 
      IF(NQQ2.NE.NQ2)WRITE(*,*)"WARNING CHANGE NQ2 ACCORDING TO .LHgrid"
      READ(1,*) q2min 
      DO iq2 = 1,nqq2 
         READ(1,*) q2g(iq2,iset) 
      ENDDO 
                                                                        
!     Read the number of replicas                                       
      READ(1,*) NREP 
                                                                        
! - dummy read up rto the member requested               

      DO K=0,mem-1 
         DO IX=1,NX 
            DO IQ2=1,NQ2 
               READ(1,'(a)') line 
            ENDDO 
         ENDDO 
      ENDDO 


!- read in the data of the requested member 
      DO IX=1,NX 
         DO IQ2=1,NQ2 
            READ(1,*) tmp
            DO ipdf=-6,6,1
               xpdfev(ix,iq2,ipdf,iset-1) = tmp(ipdf+7)
            ENDDO
         ENDDO 
      ENDDO

!                                                                       
      close(1) 
      RETURN 
        
      END                                           
                                                                        
                                                                        
!****************************************************************       
!                                                                       
! CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.     
! INTERPOLATION SUBROUTINES ARE TAKEN FROM                              
! G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,                           
! COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977). 
!                                                                       
! SUBROUTINE TAKEN FROM AAC GROUP (KUMANO et al.)                       
!                                                                       
!*******************************************************************    
!                                                                       
      SUBROUTINE LH_JSPLINE(B,C,D,I,J,KREP) 
      IMPLICIT none 
      include 'parmsetup.inc' 
!                                                                       
      INTEGER MXREP 
      PARAMETER(MXREP=nmxset-1) 
      INTEGER NREP 
      common/nnpdf10CNREP/NREP 
!                                                                       
      INTEGER NX,NQ2,NPL,IX,IQ2 
      PARAMETER(NX=60,NQ2=50) 
      PARAMETER(NPL=3000) 
      REAL*8 XG(NX,nmxset),Q2G(NQ2,nmxset),XPDFEV(NX,NQ2,-6:6,0:MXREP) 
      common/nnpdf10CPDFGR/XPDFEV,XG,Q2G,IX,IQ2 
!                                                                       
      INTEGER I,J,NM1,K,IB,KREP 
      REAL*8 B,C,D,T 
      DIMENSION B(NX,NQ2,-6:6), C(NX,NQ2,-6:6), D(NX,NQ2,-6:6) 
      integer iset
!                                                                       
      call getnset(iset)
      NM1=NX-1 
      IF(NX.LT.2) RETURN 
      IF(NX.LT.3) GOTO 250 
      D(1,J,I)=XG(2,iset)-XG(1,iset) 
      C(2,J,I)=(XPDFEV(2,J,I,KREP)-XPDFEV(1,J,I,KREP))/D(1,J,I) 
      DO 210 K=2,NM1 
         D(K,J,I)=XG(K+1,iset)-XG(K,iset) 
         B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I)) 
         C(K+1,J,I)=(XPDFEV(K+1,J,I,KREP)-XPDFEV(K,J,I,KREP))/D(K,J,I) 
         C(K,J,I)=C(K+1,J,I)-C(K,J,I) 
  210 END DO 
      B(1,J,I)=-D(1,J,I) 
      B(NX,J,I)=-D(NX-1,J,I) 
      C(1,J,I)=0.0D0 
      C(NX,J,I)=0.0D0 
      IF(NX.EQ.3) GOTO 215 
      C(1,J,I)=C(3,J,I)/(XG(4,iset)-XG(2,iset))-C(2,J,I)/(XG(3,iset)-XG(1,iset)) 
      C(NX,J,I)=C(NX-1,J,I)/(XG(NX,iset)-XG(NX-2,iset))                           &
     &     -C(NX-2,J,I)/(XG(NX-1,iset)-XG(NX-3,iset))                             
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(XG(4,iset)-XG(1,iset)) 
      C(NX,J,I)=-C(NX,J,I)*D(NX-1,J,I)**2.0D0/(XG(NX,iset)-XG(NX-3,iset)) 
  215 CONTINUE 
      DO 220 K=2,NX 
         T=D(K-1,J,I)/B(K-1,J,I) 
         B(K,J,I)=B(K,J,I)-T*D(K-1,J,I) 
         C(K,J,I)=C(K,J,I)-T*C(K-1,J,I) 
  220 END DO 
      C(NX,J,I)=C(NX,J,I)/B(NX,J,I) 
      DO 230 IB=1,NM1 
         K=NX-IB 
         C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I) 
  230 END DO 
      B(NX,J,I)=(XPDFEV(NX,J,I,KREP)-XPDFEV(NM1,J,I,KREP))/D(NM1,J,I)   &
     &     +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(NX,J,I))                     
      DO 240 K=1,NM1 
         B(K,J,I)=(XPDFEV(K+1,J,I,KREP)-XPDFEV(K,J,I,KREP))/D(K,J,I)    &
     &        -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))                     
         D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I) 
         C(K,J,I)=3.0D0*C(K,J,I) 
  240 END DO 
      C(NX,J,I)=3.0D0*C(NX,J,I) 
      D(NX,J,I)=D(NX-1,J,I) 
      RETURN 
  250 CONTINUE 
      B(1,J,I)=(XPDFEV(2,J,I,KREP)-XPDFEV(1,J,I,KREP))/(XG(2,iset)-XG(1,iset)) 
      C(1,J,I)=0.0D0 
      D(1,J,I)=0.0D0 
      B(2,J,I)=B(1,J,I) 
      C(2,J,I)=0.0D0 
      D(2,J,I)=0.0D0 
      RETURN 
      END                                           
!                                                                       
!***********************************************************************
!     THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION           
!     X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.                      
!                                                                       
!     FUNCTION TAKEN FROM AAC GROUP (KUMANO et al.)                     
!***********************************************************************
                                                                        
      INTEGER FUNCTION LH_JISEARCH(N,X,Y) 
!                                                                       
      IMPLICIT REAL*8(A-H,O-Z) 
!     Dynamical memory allocation                                       
      REAL*8 X(*) 
!                                                                       
      MIN=1 
      MAX=N+1 
!                                                                       
   10 CONTINUE 
      MID=(MIN+MAX)/2 
      IF(Y.LT.X(MID)) THEN 
        MAX=MID 
      ELSE 
        MIN=MID 
      END IF 
      IF((MAX-MIN).GT.1) GOTO 10 
!                                                                       
      LH_JISEARCH=MIN 
!                                                                       
      RETURN 
      END                                           
                                                                        
!****************************************************                   
!                                                                       
!     polin2.f                                                          
!                                                                       
!     2D interpolation of arbitrary polinomial order                    
!     Uses polint                                                       
!     Given arrays x1a(1:m) and x2a(1:n) of independent variables,      
!     and an m by n array of function values ya(1:m,1:n) tabulated      
!     at the grid points defined by x1a,x2a; and given values x1,x2     
!     of the independent variable, this routine returns                 
!     an interpolated function value y with error dy                    
!                                                                       
!     Taken from NR fortran                                             
!                                                                       
!****************************************************                   
                                                                        
      subroutine lh_polin2(x1a,x2a,ya,m,n,x1,x2,y,dy) 
      implicit none 
!                                                                       
      integer m,n,nmax,mmax 
      integer j,k 
      parameter(nmax=1e3,mmax=1e3) 
                                                                        
      real*8 dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax) 
      real*8 ymtmp(nmax),yntmp(nmax) 
                                                                        
      do j=1,m 
         do k=1,n 
            yntmp(k)=ya(j,k) 
         enddo 
         call lh_polint(x2a,yntmp,n,x2,ymtmp(j),dy) 
      enddo 
      call lh_polint(x1a,ymtmp,m,x1,y,dy) 
!                                                                       
      return 
      END                                           
                                                                        
!**********************************************                         
!                                                                       
!     polint.f                                                          
!                                                                       
!     Order N polynomial interpolation using Lagrange's formula         
!     as descrived in Numerical Recipees:                               
!     Given arrays xa and ya each of length n, and given a value        
!     x, this routine returns a value y and an error estimate dy.       
!     If P(x) is the polynomial of degree N-1 such that                 
!     P(xa_i)=ya_i,i=1,...,n, then the returned value is y=P(x)         
!     The algorithm used is Neville's algorithm                         
!                                                                       
!******************************************************                 
                                                                        
      subroutine LH_POLINT(xa,ya,n,x,y,dy) 
      implicit none 
!                                                                       
      integer n,NMAX 
!     Largest anticipated value of n                                    
      parameter(nmax=1e3) 
      real*8 dy,x,y,xa(nmax),ya(nmax) 
      integer i,m,ns 
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax) 
      ns=1 
      dif=abs(x-xa(1)) 
      do 11 i=1,n 
         dift=abs(x-xa(i)) 
         if(dift.lt.dif) then 
            ns=i 
            dif=dift 
         endif 
         c(i)=ya(i) 
         d(i)=ya(i) 
   11 enddo 
      y=ya(ns) 
      ns=ns-1 
      do m=1,n-1 
         do i=1,n-m 
            ho=xa(i)-x 
            hp=xa(i+m)-x 
            w=c(i+1)-d(i) 
            den=ho-hp 
            if(den.eq.0) then 
               write(*,*)'failure in polint' 
               stop 
            endif 
            den=w/den 
            d(i)=hp*den 
            c(i)=ho*den 
         enddo 
         if(2*ns.lt.(n-m)) then 
            dy=c(ns+1) 
         else 
            dy=d(ns) 
            ns=ns-1 
         endif 
         y=y+dy 
      enddo 
                                                                        
      return 
      END                                           
