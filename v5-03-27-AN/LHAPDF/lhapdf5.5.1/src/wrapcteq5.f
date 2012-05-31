! -*- F90 -*-


      subroutine CTEQ5evolve(x,Q,pdf) 
      implicit real*8(a-h,o-z) 
      include 'parmsetup.inc' 
      character*16 name(nmxset) 
      integer nmem(nmxset),ndef(nmxset),mmem 
      common/NAME/name,nmem,ndef,mmem 
      integer nset 
      real*8 pdf(-6:6) 
      Character Line*80 
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6) 
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX) 
      PARAMETER (M= 2, M1 = M + 1) 
      Common                                                            &
     & / CtqPar1 /                                                      &
     & Al5(nmxset), XV5(0:MXX,nmxset), QL5(0:MXQ,nmxset),               &
     & UPD5(MXPQX,nmxset)                                               &
     & / CtqPar2 / Nx(nmxset), Nt(nmxset), NfMx(nmxset)                 &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             &
     & / QCDtable /  Alambda(nmxset), Nfl(nmxset), Iorder(nmxset)       &
     & / Masstbl / Amass(6,nmxset)                                      
      data pi / 3.141592653589793d0 / 
      save 
!                                                                       
      U =         X * CtLhCtq5Pdf(1,X,Q) 
      D =         X * CtLhCtq5Pdf(2,X,Q) 
      USEA =      X * CtLhCtq5Pdf(-1,X,Q) 
      DSEA =      X * CtLhCtq5Pdf(-2,X,Q) 
      STR =       X * CtLhCtq5Pdf(3,X,Q) 
      CHM =       X * CtLhCtq5Pdf(4,X,Q) 
      BOT =       X * CtLhCtq5Pdf(5,X,Q) 
      GLU  =      X * CtLhCtq5Pdf(0,X,Q) 
      UPV=U-USEA 
      DNV=D-DSEA 
!                                                                       
      pdf(0)  = glu 
      pdf(1)  = dnv+dsea 
      pdf(-1) = dsea 
      pdf(2)  = upv+usea 
      pdf(-2) = usea 
      pdf(3)  = str 
      pdf(-3) = str 
      pdf(4)  = chm 
      pdf(-4) = chm 
      pdf(5)  = bot 
      pdf(-5) = bot 
      pdf(6)  = 0.0d0 
      pdf(-6) = 0.0d0 
                                                                        
      return 
!                                                                       
      entry CTEQ5read(nset) 
                                                                        
      call CtLhbldat1 
      call CtLhbldat2 
                                                                        
      read(1,*)nmem(nset),ndef(nset) 
      Read  (1, '(A)') Line 
      Read  (1, '(A)') Line 
      Read  (1, *) Dr, Fl, Al5(nset), (Amass(I,nset),I=1,6) 
      Iorder(nset) = Nint(Dr) 
      Nfl(nset) = Nint(Fl) 
      Alambda(nset) = Al5(nset) 
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) NX(nset),  NT(nset), NfMx(nset) 
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) QINI(nset), QMAX(nset), (QL5(I,nset), I =0, NT(nset)) 
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) XMIN(nset), (XV5(I,nset), I =0, NX(nset)) 
                                                                        
      Do 11 Iq = 0, NT(nset) 
         QL5(Iq,nset) = Log (QL5(Iq,nset) /Al5(nset)) 
   11 Continue 
!                                                                       
!                  Since quark = anti-quark for nfl>2 at this stage,    
!                  we Read  out only the non-redundent data points      
!     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)              
                                                                        
      Nblk = (NX(nset)+1) * (NT(nset)+1) 
      Npts =  Nblk  * (NfMx(nset)+3) 
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *, IOSTAT=IRET) (UPD5(I,nset), I=1,Npts) 
      Return 
!                                                                       
                                                                        
      entry CTEQ5alfa(alfas,Qalfa) 
                                                                        
       alfas = pi*CtLhALPI(Qalfa) 
                                                                        
      return 
!                                                                       
      entry CTEQ5init(Eorder,Q2fit) 
                                                                        
      return 
!                                                                       
                                         !error corrected (jcp)         
      entry CTEQ5pdf(mem) 
!        imem = mem                                                     
        call getnset(iset) 
        call setnmem(iset,mem) 
      return 
!                                                                       
      END                                           
                                                                        
      FUNCTION CtLhPartonX5 (IPRTN, X, Q) 
!                                                                       
!   Given the parton distribution function in the array UPD in          
!   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value o
!   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lamb
!                                                                       
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
!                                                                       
      include 'parmsetup.inc' 
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6) 
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX) 
      PARAMETER (M= 2, M1 = M + 1) 
!                                                                       
      Common                                                            &
     & / CtqPar1 / Al5(nmxset), XV5(0:MXX,nmxset), QL5(0:MXQ,nmxset),   &
     & UPD5(MXPQX,nmxset)                                               &
     & / CtqPar2 / Nx(nmxset), Nt(nmxset), NfMx(nmxset)                 &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             
!                                                                       
      Dimension Fq(M1), Df(M1) 
                                                                        
! note, this routine doesn't have the features of not recalculating     
! the x or q point values if they have not changed since the last call. 
! that makes it slower than cteq6, but since cteq5 is already obsolete, 
! speed is not so important for it.                                     
                                                                        
      call getnset(iset) 
                                                                        
!                                                 Work with Log (Q)     
      QG  = LOG (Q/AL5(iset)) 
                                                                        
!                           Find lower end of interval containing X     
      JL = -1 
      JU = Nx(iset)+1 
   11 If (JU-JL .GT. 1) Then 
         JM = (JU+JL) / 2 
         If (X .GT. XV5(JM,iset)) Then 
            JL = JM 
         Else 
            JU = JM 
         Endif 
         Goto 11 
      Endif 
                                                                        
      Jx = JL - (M-1)/2 
! crude treatment if outside the defined interval.                      
      If (Jx .LT. 0) then 
         Jx = 0 
      Elseif (Jx .GT. Nx(iset)-M) Then 
         Jx = Nx(iset) - M 
      Endif 
!                                    Find the interval where Q lies     
      JL = -1 
      JU = NT(iset)+1 
   12 If (JU-JL .GT. 1) Then 
         JM = (JU+JL) / 2 
         If (QG .GT. QL5(JM,iset)) Then 
            JL = JM 
         Else 
            JU = JM 
         Endif 
         Goto 12 
      Endif 
                                                                        
      Jq = JL - (M-1)/2 
      If (Jq .LT. 0) Then 
         Jq = 0 
!         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))',                    
!     >     ' WARNING: Q << Qini, extrapolation used; Q,Qini=',Q,Qini   
      Elseif (Jq .GT. Nt(iset)-M) Then 
         Jq = Nt(iset) - M 
!         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))',                    
!     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif 
                                                                        
      If (Iprtn .GE. 3) Then 
         Ip = - Iprtn 
      Else 
         Ip = Iprtn 
      EndIf 
!                             Find the off-set in the linear array Upd  
      JFL = Ip + NfMx(iset) 
      J0  = (JFL * (NT(iset)+1) + Jq) * (NX(iset)+1) + Jx 
!                                                                       
!                                           Now interpolate in x for M1 
      Do 21 Iq = 1, M1 
         J1 = J0 + (Nx(iset)+1)*(Iq-1) + 1 
         Call CtLhPolint3                                               &
     &   (XV5(Jx,iset), Upd5(J1,iset), M1, X, Fq(Iq), Df(Iq))           
   21 Continue 
!                                          Finish off by interpolating i
      Call CtLhPolint3 (QL5(Jq,iset), Fq(1), M1, QG, Ftmp, Ddf) 
                                                                        
      CtLhPartonX5 = Ftmp 
!                                                                       
      RETURN 
                                                                        
!                        ****************************                   
      END                                           
                                                                        
      Function CtLhCtq5Pdf (Iparton, X, Q) 
      Implicit Double Precision (A-H,O-Z) 
      include 'parmsetup.inc' 
      Logical Warn 
      Common                                                            &
     & / CtqPar2 / Nx(nmxset), Nt(nmxset), NfMx(nmxset)                 &
     & / QCDtable /  Alambda(nmxset), Nfl(nmxset), Iorder(nmxset)       
                                                                        
      Data Warn /.true./ 
      save Warn 
                                                                        
      call getnset(iset) 
                                                                        
      If (X .lt. 0D0 .or. X .gt. 1D0) Then 
        Print *, 'X out of range in CtLhCtq5Pdf: ', X 
        Stop 
      Endif 
      If (Q .lt. Alambda(iset)) Then 
        Print *, 'Q out of range in CtLhCtq5Pdf: ', Q 
        Stop 
      Endif 
      If ((Iparton .lt. -NfMx(iset) .or. Iparton .gt. NfMx(iset))) Then 
         If (Warn) Then 
!        put a warning for calling extra flavor.                        
             Warn = .false. 
             Print *, 'Warning: Iparton out of range in CtLhCtq5Pdf: '  &
     &              , Iparton                                           
         Endif 
         CtLhCtq5Pdf = 0D0 
         Return 
      Endif 
                                                                        
      CtLhCtq5Pdf = CtLhPartonX5 (Iparton, X, Q) 
      if(CtLhCtq5Pdf.lt.0.D0)  CtLhCtq5Pdf = 0.D0 
                                                                        
      Return 
                                                                        
!                             ********************                      
      END                                           
