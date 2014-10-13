!***********************************************                        
!                                                                       
!                                                                       
!     wrapNNPDF20grid-lite.f
!     Special low-memory version with only a single ffit ....                                                   
!     Routine called by LHAPDF package for calculating                  
!     the value of all xpdfs at x and Q from replica KREP               
!     in a (x,Q) point as called by NNPDF.LHgrid file.                  
!                                                                       
!     In 'wrapevolve.f' the package calls:                              
!     IF(NAME(NSET).EQ.'NNPDF20int') call NNPDFINT20evolve(x,Q,f)           
!     IF(NAME(NSET).EQ.'NNPDF20int') call NNPDFINT20read(nset)              
!     IF(NAME(NSET).EQ.'NNPDF20int') call NNPDFINT20alfa(alfas,Q)           
!     IF(NAME(NSET).EQ.'NNPDF20int') call NNPDFINT20init(nset,Eorder,Q2fit)*
!                                                                       
!***********************************************                                                                                                
                                                                      
      subroutine NNPDFINT20evolve(XIN,QIN,XPDF) 
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
      common/nnpdf20CNREP/NREP 
!                                                                       
      INTEGER NX,NQ2,NPL 
      PARAMETER(NX=100,NQ2=50) 
      PARAMETER(NPL=5000) 
      INTEGER NXX,NQQ2 
      REAL*8 Q2MIN,Q2MAX,XPDFMIN,XPDFMAX 
      REAL*8 XG(NX,nmxset),Q2G(NQ2,nmxset),XPDFEV(NX,NQ2,-6:6,0:MXREP) 
      INTEGER IX,IQ2 
!     This common different from NNPDF1.X
      common/nnpdf20CPDFGR/XPDFEV,XG,Q2G,IX,IQ2 
      double precision XPDFMIN_INTER
      parameter(XPDFMIN_INTER=1d-7)
!     
      INTEGER ipt,imodev,ivfn,itmc 
!     This common is the same as for NNPDF1.0
      COMMON/nnpdf10EVFLAGS/ipt,imodev,ivfn,itmc 
      REAL*8 q0,alfas0 
      REAL*8 q20,qth(4:6)
!     This common is the same as for NNPDF1.0
      COMMON/nnpdf10EVSCALE/q20,q2 
      REAL*8 q2th(4:6),asref,q2ref 
!     This common is the same as for NNPDF1.0
      COMMON/nnpdf10vfns/q2th,asref,q2ref 
!                                                                       
      INTEGER I,J,K 
      INTEGER IPDF,KREP
      INTEGER IDUM,JDUM 
      REAL*8 X,XIN,Q,QIN,QAS,QQ,Q2,QQ2,XPDF(-6:6) 
      REAL*8 XCH 
      PARAMETER(XCH=1D-1) 
!                                                                       
      integer m,n,nmax,mmax,minq,maxq,midq,maxx,minx,midx
!     order of pol. interpolation           
      parameter(m=4,n=2) 
      parameter(nmax=1e3,mmax=1e3) 
      double precision dy,x1,x2,y,x1a(mmax),x2a(nmax),ya(mmax,nmax) 
      integer ix1a(m),ix2a(n) 

!-----------------newline------------
      character*512 filename
      common/lhafilename/filename
      double precision EPS_MC_set
      parameter(EPS_MC_set=1d-7)
!-----------------newline------------
!                                                                       
!      KREP = mmem
!  always use the 0 element in the lite version  
      KREP = 0 
      X = XIN
      Q = QIN
                                                                        
!     Set correct scale                                                 
      Q2=Q**2d0 
                                                                        
!     Check kinematic point is within allowed range                     
      call getnset(iset)                                                                  
      call GetXminM(iset,KREP,XPDFMIN) 
      call GetXmaxM(iset,KREP,XPDFMAX) 
      call GetQ2maxM(iset,KREP,Q2MAX) 
      call GetQ2minM(iset,KREP,Q2MIN) 
!                                                                       
      IF ( X.LT.XPDFMIN_INTER .OR. X.GT.XPDFMAX ) THEN 
         WRITE(6,2000) 
 2000    FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE -- STOP') 
         write(6,*) "X= ",X," XMAX, XMIN = ",XPDFMAX,XPDFMIN_INTER 
         IF ( X.LT.XPDFMIN_INTER ) X = XPDFMIN_INTER 
         IF ( X.GT.XPDFMAX ) X = XPDFMAX 
!         call exit(-10)
      ENDIF 
!                                                                       
      IF ( Q2.LT.Q2MIN .OR. Q2.GT.Q2MAX ) THEN 
         WRITE(6,2001) 
 2001    FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE -- STOP') 
         write(6,*) "Q2 ,Q2MIN, Q2MAX = ",Q2,Q2MIN,Q2MAX 
         IF ( Q2.LT.Q2MIN ) Q2 = Q2MIN 
         IF ( Q2.GT.Q2MAX ) Q2 = Q2MAX 
!         call exit(-10)
      ENDIF 

!     FIND NEAREST POINTS IN THE GRID        
      MINX = 1
      MAXX = NX+1
 10   CONTINUE
      MIDX = (MINX+MAXX)/2
      IF(X.LT.XG(MIDX,iset)) THEN
         MAXX=MIDX
      ELSE
         MINX=MIDX
      END IF
      IF((MAXX-MINX).GT.1) GO TO 10
      IX = MINX

      MINQ = 1
      MAXQ = NQ2+1
 20   CONTINUE
      MIDQ = (MINQ+MAXQ)/2
      IF(Q2.LT.Q2G(MIDQ,iset)) THEN
         MAXQ=MIDQ
      ELSE
         MINQ=MIDQ
      END IF
      IF((MAXQ-MINQ).GT.1) GO TO 20
      IQ2 = MINQ

!
!     POLYNOMIAL INTERPOLATION
!        
!     uncomment for 3rd order interp.                                                                        
!     IF((IX+(M-1)/2).GT.NX) IX = NX - (M-1)/2                      
!     IF((IX-(M-1)/2).LT.1) IX = (M+1)/2                            
!     IDUM = 0                                                      
!     DO I = -(M-1)/2,(M-1)/2,1                                     
!     IDUM = IDUM +1                                             
!     IX1A(IDUM) = IX + I                                        
!     ENDDO                                                         
!     IF((IQ2+(N-1)/2).GT.NQ2) IQ2 = NQ2 - (N-1)/2                  
!     IF((IQ2-(N-1)/2).LT.1) IQ2 = (N+1)/2                          
!     JDUM = 0                                                      
!     DO J = -(N-1)/2,(N-1)/2,1                                     
!     JDUM = JDUM +1                                             
!     IX2A(JDUM) = IQ2 + J                                       
!     ENDDO                                                         

!     Assign grid for interpolation. M, N -> order of polyN interpolation      
      do I=1,M
         if(IX.ge.M/2.and.IX.le.(NX-M/2)) IX1A(I) = IX - M/2 + I
         if(IX.lt.M/2) IX1A(I) = I
         if(IX.gt.(NX-M/2)) IX1A(I) = (NX - M) + I
         
!     Check grids
         if(IX1A(I).le.0.or.IX1A(I).gt.NX) then
            write(6,*) "Error in grids! "
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
            write(6,*) "Error in grids! "
            write(6,*) "J, IXIA(J) = ",J,IX2A(J)
            call exit(-10)
         endif
      enddo
            
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

!-----------------newline------------
!     Here we need a IF switch that only activates this
!     when the _mc.LHgrid files are used
!     
!     PDF positivity for NLO MC PDFs
      if(index(filename(1:len_trim(filename)),"_mc.LHgrid").gt.0) then
         if( XPDF(IPDF).le.0d0 ) XPDF(IPDF)=EPS_MC_set
      endif 
!-----------------newline-----------------------            
         
      enddo                 

      RETURN
                                                                        
!********************************************************               
      entry NNPDFINT20getgrid(nset,ngridx,ngridq,gridx,gridq)
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
                                                                        
      ENTRY NNPDFINT20read(nset) 
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
                                                                        
!     Dummy read in to get to End (stream 1 is still open)                       
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
                                                                        
      ENTRY NNPDFINT20alfa(alfas,QAS) 
!                                                                       
      QQ = QAS 
      alfas =  alphaNNPDF(QQ) 
!                                                                       
      RETURN 
                                                                        
!********************************************************               
                                                                        
      ENTRY NNPDFINT20init(nset,Eorder,Q2fit) 
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
      entry NNPDFINT20pdf(mem)
! have to reopen stream 1
      call getnset(iset)
      call setnmem(iset,mem)
      
      call getsetpath(setpath)
      
      open(1,file=setpath(1:len_trim(setpath)),action='READ')
      line = ''
      do while (line(3:12).ne.'Evolution:'.or.len_trim(line).gt.15)
          read(1,'(a)')line
      enddo
      read(1,'(a)')line
      read(1,'(a)')line
          
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
                                                                        
      READ(1,*) NREP 
                                                                        
!     Dummy read up to the member requested the number of replicas                                       
      DO K=0,mem-1 
         DO IX=1,NX 
            DO IQ2=1,NQ2 
               READ(1,'(a)') line 
            ENDDO 
         ENDDO 
      ENDDO 
 
 !    Read in the data of the requested member                                      
      DO IX=1,NX 
         DO IQ2=1,NQ2 
            READ(1,*) ( xpdfev(ix,iq2,ipdf,0+iset-1), ipdf=-6,6,1 ) 
         ENDDO 
      ENDDO 
      
      close(1)
      RETURN                                                                        
!
      END                                           
                                                                        

!     THE ROUTINES FOR THE POLYNOMIAL INTERPOLATION ARE INSIDE
!     THE wrapNNPDFgrid.f ROUTINE
