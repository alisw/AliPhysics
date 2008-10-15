! -*- F90 -*-


      subroutine CTEQ6lgevolve(x,Q,pdf,gluino) 
      implicit real*8(a-h,o-z) 
      include 'parmsetup.inc' 
      character*16 name(nmxset) 
      integer nmem(nmxset),ndef(nmxset),mmem 
      common/NAME/name,nmem,ndef,mmem 
      real*8 pdf(-6:6) 
      integer nset 
      Character Line*80 
      PARAMETER (MXX = 105, MXQ = 51, MXF = 6, nhess = 8) 
!gino                  No of flavors = 2*NfMx + 1 gluon + 1 gluino      
      PARAMETER (MXPQX = (2*MXF + 2) * MXQ * MXX) 
      Common                                                            &
     & / CtqPar1nhesslg /                                               &
     &  Al(0:nhess,nmxset),                                             &
     &  XV(0:nhess,0:MXX,nmxset),                                       &
     &  TV(0:nhess,0:MXQ,nmxset),                                       &
     &  UPD(0:nhess,MXPQX,nmxset)                                       &
     & / CtqPar2lg / Nx(0:nhess,nmxset),                                &
     &   Nt(0:nhess,nmxset),                                            &
     &   NfMx(0:nhess,nmxset)                                           &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             &
     & / QCDtablelg /                                                   &
     &   Alambda(0:nhess,nmxset), Nfl(nmxset), Iorder(nmxset)           &
     & / Masstbllg / Amass(6,nmxset),Alg(0:nhess)                       
      common/masses_LHA/cMass(nmxset),bMass(nmxset),tMass(nmxset) 
      data pi / 3.141592653589793d0 / 
      save 
!                                                                       
      call getnset(iset) 
      call getnmem(iset,imem) 
!                                                                       
      U =         X * CtLhCtq6lgPdf(imem,1,X,Q) 
      D =         X * CtLhCtq6lgPdf(imem,2,X,Q) 
      USEA =      X * CtLhCtq6lgPdf(imem,-1,X,Q) 
      DSEA =      X * CtLhCtq6lgPdf(imem,-2,X,Q) 
      STR =       X * CtLhCtq6lgPdf(imem,3,X,Q) 
      CHM =       X * CtLhCtq6lgPdf(imem,4,X,Q) 
      BOT =       X * CtLhCtq6lgPdf(imem,5,X,Q) 
      GLU  =      X * CtLhCtq6lgPdf(imem,0,X,Q) 
      GLUINO   =  X * CtLhCtq6lgPdf(imem,-6,X,Q) 
!                                                                       
      pdf(0)  = glu 
      pdf(1)  = d 
      pdf(-1) = dsea 
      pdf(2)  = u 
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
      entry CTEQ6lgread(nset) 
                                                                        
      call CtLhbldat1 
      call CtLhbldat2 
                                                                        
      call LHct6set 
                                            !*** nmem+1=number of member
      read(1,*)nmem(nset),ndef(nset) 
      if(nmem(nset) .gt. nhess) then 
         print *,'fatal error:  nmem=',nmem(nset),' > nhess=',nhess 
         stop 
      endif 
                                             !*** new version: allows nm
      do ihess = 0,nmem(nset) 
        Read  (1, '(A)') Line 
        Read  (1, '(A)') Line 
        Read  (1, '(A)') Line 
        Read(1, *)Dr,Fl,Al(ihess,nset),(Amass(I,nset),I=1,6),Alg(ihess) 
        Iorder(nset) = Nint(Dr) 
        Nfl(nset) = Nint(Fl) 
        Alambda(ihess,nset) = Al(ihess,nset) 
                                                                        
        cMass(nset) = Amass(4,nset) 
        bMass(nset) = Amass(5,nset) 
        tMass(nset) = Amass(6,nset) 
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, *) NX(ihess,nset),  NT(ihess,nset), NfMx(ihess,nset) 
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, *)                                                    &
     &  QINI(nset), QMAX(nset), (TV(ihess,I,nset), I =0, NT(ihess,nset))
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, *) XMIN(nset), (XV(ihess,I,nset), I =0,NX(ihess,nset)) 
        Do 11 Iq = 0, NT(ihess,nset) 
           TV(ihess,Iq,nset) =Log(Log(TV(ihess,Iq,nset)/Al(ihess,nset))) 
   11   Continue 
!                                                                       
!                  Since quark = anti-quark for nfl>2 at this stage,    
!                  we Read  out only the non-redundent data points      
!     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)              
                                                                        
        Nblk = (NX(ihess,nset)+1) * (NT(ihess,nset)+1) 
        Npts =  Nblk  * (2*NfMx(ihess,nset)+2) 
                                                                        
                                                                        
!        Read  (1, '(A)') Line                                          
        Read  (1, '(A)') Line 
        Read  (1, *, IOSTAT=IRET) (UPD(ihess,I,nset), I=1,Npts) 
                                                                        
      enddo 
      return 
!                                                                       
      entry CTEQ6lgalfa(alfas,Qalfa) 
        alfas = CtLhalgluino(Qalfa) 
      return 
!                                                                       
      entry CTEQ6lginit(Eorder,Q2fit) 
                                                                        
      return 
!                                                                       
      entry CTEQ6lgpdf(mem) 
!        imem = mem                                                     
        call getnset(iset) 
        call setnmem(iset,mem) 
      return 
!                                                                       
      END                                           
                                                                        
!      subroutine LHct6set                                              
!      Implicit Double Precision (A-H,O-Z)                              
!      common /ctq6co/ xlast, qlast, nxsave                             
!      nxsave = -1000                                                   
!      xlast = -2.                                                      
!      qlast = -2.                                                      
!      return                                                           
!      end                                                              
                                                                        
!=======================================================================
      Function CtLhPartonX6lg (imem,IPRTN, XX, QQ) 
!  Given the parton distribution function in the array U in             
!  COMMON / PEVLDT / , this routine interpolates to find                
!  the parton distribution at an arbitray point in x and q.             
!                                                                       
      Implicit Double Precision (A-H,O-Z) 
                                                                        
      include 'parmsetup.inc' 
      PARAMETER (MXX = 105, MXQ = 51, MXF = 6, nhess = 8) 
!gino                  No of flavors = 2*NfMx + 1 gluon + 1 gluino      
      PARAMETER (MXPQX = (2*MXF + 2) * MXQ * MXX) 
                                                                        
      Common                                                            &
     & / CtqPar1nhesslg /                                               &
     &  Al(0:nhess,nmxset),                                             &
     &  XV(0:nhess,0:MXX,nmxset),                                       &
     &  TV(0:nhess,0:MXQ,nmxset),                                       &
     &  UPD(0:nhess,MXPQX,nmxset)                                       &
     & / CtqPar2lg / Nx(0:nhess,nmxset),                                &
     &   Nt(0:nhess,nmxset),                                            &
     &   NfMx(0:nhess,nmxset)                                           &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             
                                                                        
      common /ctq6co/ xlast, qlast, nxsave 
        parameter(nqvec = 4) 
                                                                        
      Dimension fvec(4), fij(4) 
      Dimension xvpow(0:mxx) 
      Data OneP / 1.00001 / 
                                !**** choice of interpolation variable  
      Data xpow / 0.3d0 / 
      Save xvpow 
      Data ixprint,iqprint/0,0/ 
      save ixprint,iqprint 
                                                                        
         save jq, jx, JLx, JLq, ss, sy2, sy3, s23, ty2, ty3 
         save const1 , const2, const3, const4, const5, const6 
         save tt, t13, t12, t23, t34 , t24, tmp1, tmp2, tdet 
                                                                        
      call getnset(iset) 
      call getnmem(iset,imem) 
                                                                        
                                                                        
! store the powers used for interpolation on first call...              
      if(nx(imem,iset).ne. nxsave) then 
         nxsave = nx(imem,iset) 
         xvpow(0) = 0.D0 
         do i = 1, nx(imem,iset) 
            xvpow(i) = xv(imem,i,iset)**xpow 
         enddo 
      endif 
                                                                        
      X = XX 
      Q = QQ 
                                                                        
      if((x.lt.xmin(iset)).or.(x.gt.1.d0)) then 
        ixprint=ixprint+1 
        if(ixprint.lt.11) print 98,x 
        if(ixprint.eq.10) print *,                                      &
     &  'more warning messages like the last suppressed.'               
      endif 
   98 format(' WARNING:  X=',e12.5,' OUT OF RANGE') 
      if((q.lt.qini(iset)).or.(q.gt.qmax(iset))) then 
        iqprint=iqprint+1 
        if(iqprint.lt.11) print 99,q 
        if(iqprint.eq.10) print *,                                      &
     &  'more warning messages like the last suppressed.'               
      endif 
   99 format(' WARNING:  Q=',e12.5,' OUT OF RANGE') 
                                                                        
! skip the initialization in x if same as in the previous call.         
        if(x .eq. xlast) goto 100 
        xlast = x 
                                                                        
                                                                        
                                                                        
!...new version                                                         
                                                                        
!      -------------    find lower end of interval containing x, i.e.,  
!                       get jx such that xv(jx) .le. x .le. xv(jx+1)... 
      JLx = -1 
      JU = Nx(imem,iset)+1 
   11 If (JU-JLx .GT. 1) Then 
         JM = (JU+JLx) / 2 
         If (X .Ge. XV(imem,JM,iset)) Then 
            JLx = JM 
         Else 
            JU = JM 
         Endif 
         Goto 11 
      Endif 
!                     Ix    0   1   2      Jx  JLx         Nx-2     Nx  
!                           |---|---|---|...|---|-x-|---|...|---|---|   
!                     x     0  Xmin               x                 1   
!                                                                       
      If     (JLx .LE. -1) Then 
        Print                                                           &
     &   '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6LG! x = ', x   
        Stop 
      ElseIf (JLx .Eq. 0) Then 
         Jx = 0 
      Elseif (JLx .LE. Nx(imem,iset)-2) Then 
                                                                        
!                For interrior points, keep x in the middle, as shown ab
         Jx = JLx - 1 
      Elseif (JLx.Eq.Nx(imem,iset)-1 .or. x.LT.OneP) Then 
                                                                        
!                  We tolerate a slight over-shoot of one (OneP=1.00001)
!              perhaps due to roundoff or whatever, but not more than th
!                                      Keep at least 4 points >= Jx     
         Jx = JLx - 2 
      Else 
        Print                                                           &
     &   '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6LG! x = ', x    
        Stop 
      Endif 
!          ---------- Note: JLx uniquely identifies the x-bin; Jx does n
                                                                        
!                       This is the variable to be interpolated in      
      ss = x**xpow 
                                                                        
      If (JLx.Ge.2 .and. JLx.Le.Nx(imem,iset)-2) Then 
                                                                        
!     initiation work for "interior bins": store the lattice points in s
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
                                                                        
! constants needed for interpolating in s at fixed t lattice points...  
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
                                                                        
  100 continue 
                                                                        
      if(q.eq.qlast) go to 110 
      qlast = q 
                                                                        
                                                                        
      tt = log(log(Q/Al(imem,iset))) 
                                                                        
!         --------------Now find lower end of interval containing Q, i.e
!                          get jq such that qv(jq) .le. q .le. qv(jq+1).
      JLq = -1 
      JU = NT(imem,iset)+1 
   12 If (JU-JLq .GT. 1) Then 
         JM = (JU+JLq) / 2 
         If (tt .GE. TV(imem,JM,iset)) Then 
            JLq = JM 
         Else 
            JU = JM 
         Endif 
         Goto 12 
       Endif 
                                                                        
      If     (JLq .LE. 0) Then 
         Jq = 0 
      Elseif (JLq .LE. Nt(imem,iset)-2) Then 
!                                  keep q in the middle, as shown above 
         Jq = JLq - 1 
      Else 
!                         JLq .GE. Nt-1 case:  Keep at least 4 points >=
        Jq = Nt(imem,iset) - 3 
      Endif 
                                                                        
      If (JLq.GE.1 .and. JLq.LE.Nt(imem,iset)-2) Then 
!                                        store the lattice points in t..
      tvec1 = Tv(imem,jq,iset) 
      tvec2 = Tv(imem,jq+1,iset) 
      tvec3 = Tv(imem,jq+2,iset) 
      tvec4 = Tv(imem,jq+3,iset) 
                                                                        
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
                                                                        
                                                                        
  110 continue 
                                                                        
! get the pdf function values at the lattice points...                  
                                                                        
      jtmp = ((IPRTN + NfMx(imem,iset) + 1)*(NT(imem,iset)+1)+(jq-1))   &
     & *(NX(imem,iset)+1)+jx+1                                          
                                                                        
                                                                        
      Do it = 1, nqvec 
                                                                        
         J1  = jtmp + it*(NX(imem,iset)+1) 
                                                                        
       If (Jx .Eq. 0) Then 
!                          For the first 4 x points, interpolate x^2*f(x
!                           This applies to the two lowest bins JLx = 0,
!            We can not put the JLx.eq.1 bin into the "interrior" sectio
!                           (as we do for q), since Upd(J1) is undefined
         fij(1) = 0 
         fij(2) = Upd(imem,J1+1,iset) * XV(imem,1,iset)**2 
         fij(3) = Upd(imem,J1+2,iset) * XV(imem,2,iset)**2 
         fij(4) = Upd(imem,J1+3,iset) * XV(imem,3,iset)**2 
!                                                                       
!                 Use Polint which allows x to be anywhere w.r.t. the gr
                                                                        
         Call CtLhPolint4lg (XVpow(0), Fij(1),  ss, Fx) 
                                                                        
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
!                                              Pdf is undefined for x.eq
       ElseIf  (JLx .Eq. Nx(imem,iset)-1) Then 
!                                                This is the highest x b
                                                                        
        Call CtLhPolint4lg (XVpow(Nx(imem,iset)-3), Upd(imem,J1,iset),  &
     &   ss, Fx)                                                        
                                                                        
        Fvec(it) = Fx 
                                                                        
       Else 
!                       for all interior points, use Jon's in-line funct
!                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx
         sf2 = Upd(imem,J1+1,iset) 
         sf3 = Upd(imem,J1+2,iset) 
                                                                        
         g1 =  sf2*const1 - sf3*const2 
         g4 = -sf2*const3 + sf3*const4 
                                                                        
         Fvec(it) = (const5*(Upd(imem,J1,iset)-g1)                      &
     &               + const6*(Upd(imem,J1+3,iset)-g4)                  &
     &               + sf2*sy3 - sf3*sy2) / s23                         
                                                                        
       Endif 
                                                                        
      enddo 
!                                   We now have the four values Fvec(1:4
!     interpolate in t...                                               
                                                                        
                                                                        
      If (JLq .LE. 0) Then 
!                         1st Q-bin, as well as extrapolation to lower Q
        Call CtLhPolint4lg (TV(imem,0,iset), Fvec(1),  tt, ff) 
                                                                        
      ElseIf (JLq .GE. Nt(imem,iset)-1) Then 
!                         Last Q-bin, as well as extrapolation to higher
        Call CtLhPolint4lg (TV(imem,Nt(imem,iset)-3,iset),              &
     &                      Fvec(1),  tt, ff)                           
      Else 
!                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2) 
!       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
!                         the full range QV(0:Nt)  (in contrast to XV)  
        tf2 = fvec(2) 
        tf3 = fvec(3) 
                                                                        
        g1 = ( tf2*t13 - tf3*t12) / t23 
        g4 = (-tf2*t34 + tf3*t24) / t23 
                                                                        
        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12                      &
     &          +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)                 
                                                                        
                                                                        
        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23 
      EndIf 
                                                                        
      CtLhPartonX6lg = ff 
                                                                        
      Return 
!                                       ********************            
      END                                           
!=======================================================================
                                                                        
      Function CtLhCtq6lgPdf (imem,Iparton, X, Q) 
      Implicit Double Precision (A-H,O-Z) 
      include 'parmsetup.inc' 
      PARAMETER (nhess = 8) 
      Logical Warn 
      Common                                                            &
     & / CtqPar2lg / Nx(0:nhess,nmxset), Nt(0:nhess,nmxset),            &
     &               NfMx(0:nhess,nmxset)                               &
     & / QCDtablelg /  Alambda(0:nhess,nmxset), Nfl(nmxset),            &
     &                 Iorder(nmxset)                                   
                                                                        
      Data Warn /.true./ 
      save Warn 
                                                                        
      call getnset(iset) 
                                                                        
      If (X .lt. 0D0 .or. X .gt. 1D0) Then 
        Print *, 'X out of range in CtLhCtq6Pdf: ', X 
        Stop 
      Endif 
      If (Q .lt. Alambda(imem,iset)) Then 
        Print *, 'Q out of range in CtLhCtq6Pdf: ', Q 
        Stop 
      Endif 
                                                                        
! added to force pdf = 0.0 at x=1.0 exactly - mrw                       
      if(x .eq. 1.0d0) then 
          CtLhCtq6lgPdf = 0.0d0 
          return 
      endif 
!                                                                       
      If ((Iparton .lt. -NfMx(imem,iset)-1 .or. Iparton .gt.            &
     &   NfMx(imem,iset))) Then                                         
         If (Warn) Then 
!        put a warning for calling extra flavor.                        
             Warn = .false. 
             Print *,'Warning: Iparton out of range in CtLhCtq6lgPdf: ' &
     &              , Iparton                                           
         Endif 
         CtLhCtq6lgPdf = 0D0 
         print *,' returning CtLhCtq6lgPdf = ', CtLhCtq6lgPdf 
         Return 
      Endif 
                                                                        
      CtLhCtq6lgPdf = CtLhPartonX6lg (imem,Iparton, X, Q) 
      if(CtLhCtq6lgPdf.lt.0.D0)  CtLhCtq6lgPdf = 0.D0 
      Return 
                                                                        
!                             ********************                      
      END                                           
!=======================================================================
                                                                        
      SUBROUTINE CtLhPOLINT4lg (XA,YA,X,Y) 
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
!  The POLINT4 routine is based on the POLINT routine from "Numerical Re
!  but assuming N=4, and ignoring the error estimation.                 
!  suggested by Z. Sullivan.                                            
      DIMENSION XA(*),YA(*) 
                                                                        
      H1=XA(1)-X 
      H2=XA(2)-X 
      H3=XA(3)-X 
      H4=XA(4)-X 
                                                                        
      W=YA(2)-YA(1) 
      DEN=W/(H1-H2) 
      D1=H2*DEN 
      C1=H1*DEN 
                                                                        
      W=YA(3)-YA(2) 
      DEN=W/(H2-H3) 
      D2=H3*DEN 
      C2=H2*DEN 
                                                                        
      W=YA(4)-YA(3) 
      DEN=W/(H3-H4) 
      D3=H4*DEN 
      C3=H3*DEN 
                                                                        
      W=C2-D1 
      DEN=W/(H1-H3) 
      CD1=H3*DEN 
      CC1=H1*DEN 
                                                                        
      W=C3-D2 
      DEN=W/(H2-H4) 
      CD2=H4*DEN 
      CC2=H2*DEN 
                                                                        
      W=CC2-CD1 
      DEN=W/(H1-H4) 
      DD1=H4*DEN 
      DC1=H1*DEN 
                                                                        
      If((H3+H4).lt.0D0) Then 
         Y=YA(4)+D3+CD2+DD1 
      Elseif((H2+H3).lt.0D0) Then 
         Y=YA(3)+D2+CD1+DC1 
      Elseif((H1+H2).lt.0D0) Then 
         Y=YA(2)+C2+CD1+DC1 
      ELSE 
         Y=YA(1)+C1+CC1+DC1                                             
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                           
