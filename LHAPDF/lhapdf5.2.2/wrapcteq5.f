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
      Common
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      data pi / 3.141592653589793d0 /
      save
c
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
c      
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
*
      entry CTEQ5read(nset)

      call CtLhbldat1           !this line was missing in previous releases (jcp)
      call CtLhbldat2           !this line was missing in previous releases (jcp)
 
      read(1,*)nmem(nset),ndef(nset)
      Read  (1, '(A)') Line     
      Read  (1, '(A)') Line
      Read  (1, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (1, '(A)') Line 
      Read  (1, *) NX,  NT, NfMx

      Read  (1, '(A)') Line
      Read  (1, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (1, '(A)') Line
      Read  (1, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
 
      Read  (1, '(A)') Line
      Read  (1, *, IOSTAT=IRET) (UPD(I), I=1,Npts)
      Return
*    

      entry CTEQ5alfa(alfas,Qalfa)

       alfas = pi*CtLhALPI(Qalfa)

      return
*
      entry CTEQ5init(Eorder,Q2fit)

      return
*
      entry CTEQ5pdf(mem)
        imem = mem
      return
* 
      end

      FUNCTION CtLhPartonX5 (IPRTN, X, Q)
C
C   Given the parton distribution function in the array UPD in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
c         Print '(A, 2(1pE12.4))', 
c     >     ' WARNING: X << Xmin, extrapolation used; X,Xmin=',X,Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
c         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
c     >     ' WARNING: Q << Qini, extrapolation used; Q,Qini=',Q,Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
c         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
c     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call CtLhPolint3 (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call CtLhPolint3 (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      CtLhPartonX5 = Ftmp
C
      RETURN
C                        ****************************
      END

      Function CtLhCtq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in CtLhCtq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in CtLhCtq5Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in CtLhCtq5Pdf: '
     >              , Iparton
         Endif
         CtLhCtq5Pdf = 0D0
         Return
      Endif

      CtLhCtq5Pdf = CtLhPartonX5 (Iparton, X, Q)
      if(CtLhCtq5Pdf.lt.0.D0)  CtLhCtq5Pdf = 0.D0

      Return

C                             ********************
      End
