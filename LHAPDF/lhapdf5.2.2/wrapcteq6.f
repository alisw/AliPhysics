      subroutine CTEQ6evolve(x,Q,pdf)
      implicit real*8(a-h,o-z)
      include 'parmsetup.inc'
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      real*8 pdf(-6:6)
      integer nset,k
      Character Line*80
      PARAMETER (MXX = 96, MXQ = 20, MXF = 5, nhess = 40)
      PARAMETER (MXPQX = (MXF + 3) * MXQ * MXX)
      Common
     > / CtqPar1nhess / Al, XV(0:MXX), TV(0:MXQ), UPD(0:nhess,MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      common/masses_LHA/cMass,bMass,tMass
      data pi / 3.141592653589793d0 /
      save
c
      call getnset(iset)
      call getnmem(iset,imem)
c
      U =         X * CtLhCtq6Pdf(imem,1,X,Q)
      D =         X * CtLhCtq6Pdf(imem,2,X,Q)
      USEA =      X * CtLhCtq6Pdf(imem,-1,X,Q)
      DSEA =      X * CtLhCtq6Pdf(imem,-2,X,Q)
      STR =       X * CtLhCtq6Pdf(imem,3,X,Q)
      CHM =       X * CtLhCtq6Pdf(imem,4,X,Q)
      BOT =       X * CtLhCtq6Pdf(imem,5,X,Q)
      GLU  =      X * CtLhCtq6Pdf(imem,0,X,Q)
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
      entry CTEQ6read(nset)

      call CtLhbldat1		!this line was missing in previous releases (jcp)
      call CtLhbldat2		!this line was missing in previous releases (jcp)

      read(1,*)nmem(nset),ndef(nset)	!*** nmem+1=number of members; ndef is not used for anything ***
      if(nmem(nset) .gt. nhess) then
         print *,'fatal error:  nmem=',nmem(nset),' > nhess=',nhess
         stop
      endif
      Read  (1, '(A)') Line
      Read  (1, '(A)') Line
      Read  (1, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      cMass = Amass(4)                  !missing in previous releases(jcp)
      bMass = Amass(5)                  !missing in previous releases(jcp)
      tMass = Amass(6)                  !missing in previous releases(jcp)

      Read  (1, '(A)') Line
      Read  (1, *) NX,  NT, NfMx

      Read  (1, '(A)') Line
      Read  (1, *) QINI, QMAX, (TV(I), I =0, NT)

      Read  (1, '(A)') Line
      Read  (1, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         TV(Iq) = Log(Log (TV(Iq) /Al))
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage,
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)


c **  do ihess = 0,nhess	!*** old version ***
      do ihess = 0,nmem(nset)		!*** new version: allows nmem < nhess ***
      
        Read  (1, '(A)') Line
        Read  (1, '(A)') Line
        Read  (1, *, IOSTAT=IRET) (UPD(ihess,I), I=1,Npts)

      enddo
      return
*    
      entry CTEQ6alfa(alfas,Qalfa)
      alfas = pi*CtLhALPI(Qalfa)
      return
*
      entry CTEQ6init(Eorder,Q2fit)

      return
*
      entry CTEQ6pdf(mem)
c        imem = mem
	call getnset(iset)
	call setnmem(iset,mem)
      return
* 
      end

      SUBROUTINE CtLhPOLINT (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes"
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=W/(HO-HP)
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

c===========================================================================
      Function CtLhPartonX6 (iset,IPRTN, XX, QQ)
c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      Parameter (MXX = 96, MXQ = 20, MXF = 5, nhess = 40)
      Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))

      Common
     > / CtqPar1nhess / Al, XV(0:MXX), TV(0:MXQ), UPD(0:nhess,MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx
     > / XQrange / Qini, Qmax, Xmin

	parameter(nqvec = 4)

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nxsave / -999 /
      Save nxsave, xvpow

 	save xlast, qlast
 	save jq, jx, JLx, JLq, ss, sy2, sy3, s23, ty2, ty3
 	save const1 , const2, const3, const4, const5, const6
 	save tt, t13, t12, t23, t34 , t24, tmp1, tmp2, tdet

c store the powers used for interpolation on first call...
      if(nx .ne. nxsave) then
         nxsave = nx
         xvpow(0) = 0.D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ

c skip the initialization in x if same as in the previous call.
	if(x .eq. xlast) goto 100
	xlast = x

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)','Severe error: x <= 0 in CtLhPartonX6 x=',x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)','Severe error: x > 1 in CtLhPartonX6 x=',x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

100	continue

c skip the initialization in q if same as in the previous call.
	if(q .eq. qlast) goto 110
	qlast = q

      tt = log(log(Q/Al))

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf

110	continue

c get the pdf function values at the lattice points...

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We cannot put the JLx.eq.1 bin into the "interior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(iset,J1+1) * XV(1)**2
         fij(3) = Upd(iset,J1+2) * XV(2)**2
         fij(4) = Upd(iset,J1+3) * XV(3)**2
C
C                 Use CtLhPolint which allows x to be anywhere w.r.t. the grid

         Call CtLhPolint4(XVpow(0), Fij(1), 4, ss, Fx, Dfx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

c** fix allow 4 consecutive elements with iset... mrw 19.9.2005
        fij(1) = Upd(iset,j1)
        fij(2) = Upd(iset,j1+1)
        fij(3) = Upd(iset,j1+2)
        fij(4) = Upd(iset,j1+3)
        Call CtLhPolint4 (XVpow(Nx-3), Fij(1), 4, ss, Fx, Dfx)

        Fvec(it) = Fx

       Else

C for all interior points, use Jon's in-line function
C This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
c (This is cubic spline interpolation, as used by cteq; it was 
c changed to polint in previous Durham releases (jcp).)
         sf2 = Upd(iset,J1+1)
         sf3 = Upd(iset,J1+2)

         Fvec(it) = (const5*(Upd(iset,J1) 
     &	                     - sf2*const1 + sf3*const2)
     &             + const6*(Upd(iset,J1+3) 
     &	                     + sf2*const3 - sf3*const4)
     &             + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call CtLhPolint4(TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call CtLhPolint4(TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      CtLhPartonX6 = ff

      Return
C                                       ********************
      End
c===========================================================================
      Function CtLhCtq6Pdf (iset,Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in CtLhCtq6Pdf: ', X
        Stop
      Endif
      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in CtLhCtq6Pdf: ', Q
        Stop
      Endif

c added to force pdf = 0.0 at x=1.0 exactly - mrw
      if(x .eq. 1.0d0) then
          CtLhCtq6Pdf = 0.0d0
          return
      endif
c
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in CtLhCtq6Pdf: '
     >              , Iparton
         Endif
         CtLhCtq6Pdf = 0D0
         Return
      Endif

      CtLhCtq6Pdf = CtLhPartonX6 (iset,Iparton, X, Q)
      if(CtLhCtq6Pdf.lt.0.D0)  CtLhCtq6Pdf = 0.D0

      Return

C                             ********************
      End
      SUBROUTINE CtLhPOLINT4 (XA,YA,N,X,Y,DY)
c fast version of polint, valid only for N=4
c Have explicitly unrolled the loops.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      PARAMETER (NMAX=4)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

	if(n .ne. 4) then
	   print *,'fatal CtLhPolint4 call',n
	   stop
	endif

      NS=1
      DIF=ABS(X-XA(1))

        DIFT=ABS(X-XA(1))
        IF (DIFT.LT.DIF) THEN
          NS=1
          DIF=DIFT
        ENDIF
        C(1)=YA(1)
        D(1)=YA(1)

        DIFT=ABS(X-XA(2))
        IF (DIFT.LT.DIF) THEN
          NS=2
          DIF=DIFT
        ENDIF
        C(2)=YA(2)
        D(2)=YA(2)

        DIFT=ABS(X-XA(3))
        IF (DIFT.LT.DIF) THEN
          NS=3
          DIF=DIFT
        ENDIF
        C(3)=YA(3)
        D(3)=YA(3)

        DIFT=ABS(X-XA(4))
        IF (DIFT.LT.DIF) THEN
          NS=4
          DIF=DIFT
        ENDIF
        C(4)=YA(4)
        D(4)=YA(4)


      Y=YA(NS)
      NS=NS-1


          HO=XA(1)-X
          HP=XA(2)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN


          HO=XA(2)-X
          HP=XA(3)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN


          HO=XA(3)-X
          HP=XA(4)-X
          W=C(4)-D(3)
          DEN=W/(HO-HP)
          D(3)=HP*DEN
          C(3)=HO*DEN

        IF (2*NS.LT.3)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY

          HO=XA(1)-X
          HP=XA(3)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN

          HO=XA(2)-X
          HP=XA(4)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN

        IF (2*NS.LT.2)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY

          HO=XA(1)-X
          HP=XA(4)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN

        IF (2*NS.LT.4-3)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY

      RETURN
      END
      SUBROUTINE CTLHPOLINT3 (XA,YA,N,X,Y,DY)
c fast version of polint, valid only for N=3
c Have explicitly unrolled the loops.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=3)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
	if(n .ne. 3) then
	   print *,'fatal CtLhPolint3 call',n
	   stop
	endif
      NS=1
      DIF=ABS(X-XA(1))
        DIFT=ABS(X-XA(1))
        IF (DIFT.LT.DIF) THEN
          NS=1
          DIF=DIFT
        ENDIF
        C(1)=YA(1)
        D(1)=YA(1)
        DIFT=ABS(X-XA(2))
        IF (DIFT.LT.DIF) THEN
          NS=2
          DIF=DIFT
        ENDIF
        C(2)=YA(2)
        D(2)=YA(2)
        DIFT=ABS(X-XA(3))
        IF (DIFT.LT.DIF) THEN
          NS=3
          DIF=DIFT
        ENDIF
        C(3)=YA(3)
        D(3)=YA(3)
      Y=YA(NS)
      NS=NS-1
          HO=XA(1)-X
          HP=XA(2)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN
          HO=XA(2)-X
          HP=XA(3)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN
        IF (2*NS.LT.2)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
          HO=XA(1)-X
          HP=XA(3)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN
        IF (2*NS.LT.1)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
      RETURN
      END
