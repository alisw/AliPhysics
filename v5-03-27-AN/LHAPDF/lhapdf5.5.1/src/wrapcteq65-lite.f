! -*- F90 -*-


      subroutine CTEQ65evolve(x,Q,pdf) 
      implicit real*8(a-h,o-z) 
      include 'parmsetup.inc' 
      character*16 name(nmxset) 
      character*512 setpath 
      integer nmem(nmxset),ndef(nmxset),mmem 
      common/NAME/name,nmem,ndef,mmem 
      real*8 pdf(-6:6) 
      integer nset 
      Character Line*80 
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6, MaxVal=3, nhess = 0) 
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX) 
      Common                                                            &
     & / CtqPar1nhess65 /                                               &
     &  Al(nmxset), XV(0:MXX,nmxset), TV(0:MXQ,nmxset),                 &
     &  UPD(0:nhess,MXPQX,nmxset)                                       &
     & / CtqPar2 / Nx(nmxset), Nt(nmxset), NfMx(nmxset)                 &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             &
     & / QCDtable /  Alambda(nmxset), Nfl(nmxset), Iorder(nmxset)       &
     & / Masstbl / Amass(6,nmxset)                                      
                                                                        
      common/masses_LHA/cMass(nmxset),bMass(nmxset),tMass(nmxset) 
                                                                        
      data pi / 3.141592653589793d0 / 
      save 
!                                                                       
      call getnset(iset) 
      call getnmem(iset,imem) 
                                                                        
      U =         X * CtLhCtq65Pdf(0,1,X,Q) 
      D =         X * CtLhCtq65Pdf(0,2,X,Q) 
      USEA =      X * CtLhCtq65Pdf(0,-1,X,Q) 
      DSEA =      X * CtLhCtq65Pdf(0,-2,X,Q) 
      STR =       X * CtLhCtq65Pdf(0,-3,X,Q) 
      CHM =       X * CtLhCtq65Pdf(0,-4,X,Q) 
      BOT =       X * CtLhCtq65Pdf(0,-5,X,Q) 
      GLU  =      X * CtLhCtq65Pdf(0,0,X,Q) 
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
                                          !like cteq65evolve, but allows
      entry CTEQ65cevolve(x,Q,pdf) 
!                                                                       
      call getnset(iset) 
      call getnmem(iset,imem) 
                                                                        
      U =         X * CtLhCtq65Pdf(0,1,X,Q) 
      D =         X * CtLhCtq65Pdf(0,2,X,Q) 
      USEA =      X * CtLhCtq65Pdf(0,-1,X,Q) 
      DSEA =      X * CtLhCtq65Pdf(0,-2,X,Q) 
      STR =       X * CtLhCtq65Pdf(0,3,X,Q) 
      SBAR =      X * CtLhCtq65Pdf(0,-3,X,Q) 
      CHM =       X * CtLhCtq65Pdf(0,4,X,Q) 
      CBAR =      X * CtLhCtq65Pdf(0,-4,X,Q) 
      BOT =       X * CtLhCtq65Pdf(0,5,X,Q) 
      GLU  =      X * CtLhCtq65Pdf(0,0,X,Q) 
!                                                                       
      pdf(0)  = glu 
      pdf(1)  = d 
      pdf(-1) = dsea 
      pdf(2)  = u 
      pdf(-2) = usea 
      pdf(3)  = str 
      pdf(-3) = sbar 
      pdf(4)  = chm 
      pdf(-4) = cbar 
      pdf(5)  = bot 
      pdf(-5) = bot 
      pdf(6)  = 0.0d0 
      pdf(-6) = 0.0d0 
                                                                        
      return 
!                                                                       
                                          !like cteq65evolve, but allows
      entry CTEQ65sevolve(x,Q,pdf) 
!                                                                       
      call getnset(iset) 
      call getnmem(iset,imem) 
                                                                        
      U =         X * CtLhCtq65Pdf(0,1,X,Q) 
      D =         X * CtLhCtq65Pdf(0,2,X,Q) 
      USEA =      X * CtLhCtq65Pdf(0,-1,X,Q) 
      DSEA =      X * CtLhCtq65Pdf(0,-2,X,Q) 
      STR =       X * CtLhCtq65Pdf(0,3,X,Q) 
      SBAR =      X * CtLhCtq65Pdf(0,-3,X,Q) 
      CHM =       X * CtLhCtq65Pdf(0,-4,X,Q) 
      BOT =       X * CtLhCtq65Pdf(0,5,X,Q) 
      GLU  =      X * CtLhCtq65Pdf(0,0,X,Q) 
!                                                                       
      pdf(0)  = glu 
      pdf(1)  = d 
      pdf(-1) = dsea 
      pdf(2)  = u 
      pdf(-2) = usea 
      pdf(3)  = str 
      pdf(-3) = sbar 
      pdf(4)  = chm 
      pdf(-4) = chm 
      pdf(5)  = bot 
      pdf(-5) = bot 
      pdf(6)  = 0.0d0 
      pdf(-6) = 0.0d0 
                                                                        
      return 
!                                                                       
      entry CTEQ65read(nset) 
      call CtLhbldat1 
      call CtLhbldat2 
                                                                        
      call LHct65set 
                                                                        
                                            !*** nmem+1=number of member
      read(1,*)nmem(nset),ndef(nset) 
!      if(nmem(nset) .gt. nhess) then                                   
!         print *,'fatal error:  nmem=',nmem(nset),' > nhess=',nhess    
!         stop                                                          
!      endif                                                            
                                                                        
      MxVal = 3 
      Read  (1, '(A)') Line 
      Read  (1, '(A)') Line 
      Read  (1, *) Dr, Fl, Al(nset), (Amass(I,nset),I=1,6) 
      Iorder(nset) = Nint(Dr) 
      Nfl(nset) = Nint(Fl) 
      Alambda(nset) = Al(nset) 
                                                                        
      cMass(nset) = Amass(4,nset) 
      bMass(nset) = Amass(5,nset) 
      tMass(nset) = Amass(6,nset) 
                                                                        
      Read  (1, '(A)') Line 
!                                               This is the .pds (WKT) f
      Read  (1, *) N0, N0, N0, NfMx(nset), N0, N0 
      Read  (1, '(A)') Line 
      Read  (1, *) NX(nset),  NT(nset), N0, N0, N0 
      Read  (1, '(A)') (Line,I=1,4) 
      Read  (1, *)                                                      &
     &  QINI(nset), QMAX(nset), (aa,TV(I,nset), I =0, NT(nset))         
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) XMIN(nset), aa, (XV(I,nset), I =1, NX(nset)) 
      XV(0,nset)=0D0 
!                                                                       
!                  Since quark = anti-quark for nfl>2 at this stage,    
!                  we Read  out only the non-redundent data points      
!     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)              
                                                                        
                                                                        
      Nblk = (NX(nset)+1) * (NT(nset)+1) 
      Npts =  Nblk  * (NfMx(nset)+1+MxVal) 
                                                                        
                                             !*** new version: allows nm
      do ihess = 0,nmem(nset) 
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, '(A)') Line 
        Read  (1, *, IOSTAT=IRET) (UPD(0,I,nset), I=1,Npts) 
                                                                        
      enddo 
      return 
!                                                                       
      entry CTEQ66read(nset) 
      call CtLhbldat1 
      call CtLhbldat2 
                                                                        
      call LHct65set 
                                                                        
                                            !*** nmem+1=number of member
      read(1,*)nmem(nset),ndef(nset) 
!      if(nmem(nset) .gt. nhess) then                                   
!         print *,'fatal error:  nmem=',nmem(nset),' > nhess=',nhess    
!         stop                                                          
!      endif                                                            
                                                                        
      MxVal = 2 
      Read  (1, '(A)') Line 
      Read  (1, '(A)') Line 
      Read  (1, *) Dr, Fl, Al(nset), (Amass(I,nset),I=1,6) 
      Iorder(nset) = Nint(Dr) 
      Nfl(nset) = Nint(Fl) 
      Alambda(nset) = Al(nset) 
                                                                        
      cMass(nset) = Amass(4,nset) 
      bMass(nset) = Amass(5,nset) 
      tMass(nset) = Amass(6,nset) 
                                                                        
      Read  (1, '(A)') Line 
!                                               This is the .pds (WKT) f
      Read  (1, *) N0, N0, N0, NfMx(nset), N0, N0 
      Read  (1, '(A)') Line 
      Read  (1, *) NX(nset),  NT(nset), N0, N0, N0 
      Read  (1, '(A)') (Line,I=1,4) 
      Read  (1, *)                                                      &
     &  QINI(nset), QMAX(nset), (aa,TV(I,nset), I =0, NT(nset))         
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) XMIN(nset), aa, (XV(I,nset), I =1, NX(nset)) 
      XV(0,nset)=0D0 
!                                                                       
!                  Since quark = anti-quark for nfl>2 at this stage,    
!                  we Read  out only the non-redundent data points      
!     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)              
                                                                        
                                                                        
      Nblk = (NX(nset)+1) * (NT(nset)+1) 
      Npts =  Nblk  * (NfMx(nset)+1+MxVal) 
                                                                        
                                             !*** new version: allows nm
      do ihess = 0,nmem(nset) 
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, '(A)') Line 
        Read  (1, *, IOSTAT=IRET) (UPD(0,I,nset), I=1,Npts) 
                                                                        
      enddo 
      return 
!                                                                       
                                             !like CTEQ65read, but c.ne.
      entry CTEQ65cread(nset) 
      call CtLhbldat1 
      call CtLhbldat2 
      call LHct65set 
                                                                        
                                            !*** nmem+1=number of member
      read(1,*)nmem(nset),ndef(nset) 
!      if(nmem(nset) .gt. nhess) then                                   
!         print *,'fatal error:  nmem=',nmem(nset),' > nhess=',nhess    
!         stop                                                          
!      endif                                                            
                                                                        
                                               !one more species free   
      MxVal = 4 
      Read  (1, '(A)') Line 
      Read  (1, '(A)') Line 
      Read  (1, *) Dr, Fl, Al(nset), (Amass(I,nset),I=1,6) 
      Iorder(nset) = Nint(Dr) 
      Nfl(nset) = Nint(Fl) 
      Alambda(nset) = Al(nset) 
                                                                        
      cMass(nset) = Amass(4,nset) 
      bMass(nset) = Amass(5,nset) 
      tMass(nset) = Amass(6,nset) 
                                                                        
      Read  (1, '(A)') Line 
!                                               This is the .pds (WKT) f
      Read  (1, *) N0, N0, N0, NfMx(nset), N0, N0 
      Read  (1, '(A)') Line 
      Read  (1, *) NX(nset),  NT(nset), N0, N0, N0 
      Read  (1, '(A)') (Line,I=1,4) 
      Read  (1, *)                                                      &
     &  QINI(nset), QMAX(nset), (aa,TV(I,nset), I =0, NT(nset))         
                                                                        
      Read  (1, '(A)') Line 
      Read  (1, *) XMIN(nset), aa, (XV(I,nset), I =1, NX(nset)) 
      XV(0,nset)=0D0 
!                                                                       
! we Read only the non-redundent data points                            
                                                                        
      Nblk = (NX(nset)+1) * (NT(nset)+1) 
      Npts =  Nblk  * (NfMx(nset)+1+MxVal) 
                                                                        
                                             !*** new version: allows nm
      do ihess = 0,nmem(nset) 
                                                                        
        Read  (1, '(A)') Line 
        Read  (1, '(A)') Line 
        Read  (1, *, IOSTAT=IRET) (UPD(0,I,nset), I=1,Npts) 
                                                                        
      enddo 
      return 
!                                                                       
      entry CTEQ65alfa(alfas,Qalfa) 
      alfas = pi*CtLhALPI(Qalfa) 
      return 
!                                                                       
      entry CTEQ65init(Eorder,Q2fit) 
      return 
!                                                                       
      entry CTEQ65pdf(mem) 
!        imem = mem                                                     
        call getnset(iset) 
        call setnmem(iset,mem) 
                                                                        
! have to reopen stream 1                                               
        call getsetpath(setpath) 
        open(1,file=setpath(1:len_trim(setpath))) 
                                                                        
        line = '' 
        do while (line(1:3).ne.'==>') 
          read(1,'(a)'),line 
        enddo 
! - backspace by one record                                             
        backspace(1) 
! - dummy read up to the member requested                               
        do n=0,mem-1 
          read(1,'(a)'),line 
          read(1,'(a)'),line 
          Read  (1, *, IOSTAT=IRET) (UPD(0,I,iset), I=1,Npts) 
        enddo 
!- read in the data of the requested member                             
          Read  (1, '(A)') Line 
          Read  (1, '(A)') Line 
          Read  (1, *, IOSTAT=IRET) (UPD(0,I,iset), I=1,Npts) 
          close(1) 
      return 
!                                                                       
      END                                           
                                                                        
      subroutine LHct65set 
      Implicit Double Precision (A-H,O-Z) 
      common /ctq65co/ xlast, qlast, nxsave 
      nxsave = -1000 
      xlast = -2. 
      qlast = -2. 
      return 
      END                                           
                                                                        
!=======================================================================
      Function CtLhPartonX65 (imem,IPRTN, XX, QQ) 
!  Given the parton distribution function in the array U in             
!  COMMON / PEVLDT / , this routine interpolates to find                
!  the parton distribution at an arbitray point in x and q.             
!                                                                       
      Implicit Double Precision (A-H,O-Z) 
                                                                        
      include 'parmsetup.inc' 
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6, MaxVal=3, nhess = 0) 
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX) 
                                                                        
      Common                                                            &
     & / CtqPar1nhess65 /                                               &
     &  Al(nmxset), XV(0:MXX,nmxset), TV(0:MXQ,nmxset),                 &
     &  UPD(0:nhess,MXPQX,nmxset)                                       &
     & / CtqPar2 / Nx(nmxset), Nt(nmxset), NfMx(nmxset)                 &
     & / XQrange / Qini(nmxset), Qmax(nmxset), Xmin(nmxset)             
      common /ctq65co/ xlast,qlast, nxsave 
                                                                        
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
!      call getnmem(iset,imem)                                          
                                                                        
! store the powers used for interpolation on first call...              
      if(nx(iset) .ne. nxsave) then 
         nxsave = nx(iset) 
         xvpow(0) = 0.D0 
         do i = 1, nx(iset) 
            xvpow(i) = xv(i,iset)**xpow 
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
                                                                        
!      -------------    find lower end of interval containing x, i.e.,  
!                       get jx such that xv(jx) .le. x .le. xv(jx+1)... 
      JLx = -1 
      JU = Nx(iset)+1 
   11 If (JU-JLx .GT. 1) Then 
         JM = (JU+JLx) / 2 
         If (X .Ge. XV(JM,iset)) Then 
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
        Print '(A,1pE12.4)','Severe error: x <= 0 in CtLhPartonX65 x=',x 
        Stop 
      ElseIf (JLx .Eq. 0) Then 
         Jx = 0 
      Elseif (JLx .LE. Nx(iset)-2) Then 
                                                                        
!                For interior points, keep x in the middle, as shown abo
         Jx = JLx - 1 
      Elseif (JLx.Eq.Nx(iset)-1 .or. x.LT.OneP) Then 
                                                                        
!                  We tolerate a slight over-shoot of one (OneP=1.00001)
!              perhaps due to roundoff or whatever, but not more than th
!                                      Keep at least 4 points >= Jx     
         Jx = JLx - 2 
      Else 
        Print '(A,1pE12.4)','Severe error: x > 1 in CtLhPartonX65 x=',x 
        Stop 
      Endif 
!          ---------- Note: JLx uniquely identifies the x-bin; Jx does n
                                                                        
!                       This is the variable to be interpolated in      
      ss = x**xpow 
                                                                        
      If (JLx.Ge.2 .and. JLx.Le.Nx(iset)-2) Then 
                                                                        
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
                                                                        
  100      continue 
                                                                        
! skip the initialization in q if same as in the previous call.         
        if(q .eq. qlast) goto 110 
        qlast = q 
                                                                        
      tt = log(log(Q/Al(iset))) 
                                                                        
!         --------------Now find lower end of interval containing Q, i.e
!                          get jq such that qv(jq) .le. q .le. qv(jq+1).
      JLq = -1 
      JU = NT(iset)+1 
   12 If (JU-JLq .GT. 1) Then 
         JM = (JU+JLq) / 2 
         If (tt .GE. TV(JM,iset)) Then 
            JLq = JM 
         Else 
            JU = JM 
         Endif 
         Goto 12 
       Endif 
                                                                        
      If     (JLq .LE. 0) Then 
         Jq = 0 
      Elseif (JLq .LE. Nt(iset)-2) Then 
!                                  keep q in the middle, as shown above 
         Jq = JLq - 1 
      Else 
!                         JLq .GE. Nt-1 case:  Keep at least 4 points >=
        Jq = Nt(iset) - 3 
                                                                        
      Endif 
!                                   This is the interpolation variable i
                                                                        
      If (JLq.GE.1 .and. JLq.LE.Nt(iset)-2) Then 
!                                        store the lattice points in t..
      tvec1 = Tv(jq,iset) 
      tvec2 = Tv(jq+1,iset) 
      tvec3 = Tv(jq+2,iset) 
      tvec4 = Tv(jq+3,iset) 
                                                                        
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
                                                                        
  110      continue 
                                                                        
! get the pdf function values at the lattice points...                  
! In this code, we store 10 flavors: u,ubar,d,dbar,s,sbar,c,cbar,b=bbar,
! hence Iprtn=5 (b) is obtained from -5 (bbar)                          
                                                                        
      If (Iprtn .GE. 5) Then 
         Ip = - Iprtn 
      Else 
         Ip = Iprtn 
      EndIf 
      jtmp = ((Ip + NfMx(iset))*(NT(iset)+1)+(jq-1))*(NX(iset)+1)+jx+1 
                                                                        
      Do it = 1, nqvec 
                                                                        
         J1  = jtmp + it*(NX(iset)+1) 
                                                                        
       If (Jx .Eq. 0) Then 
!                          For the first 4 x points, interpolate x^2*f(x
!                           This applies to the two lowest bins JLx = 0,
!            We cannot put the JLx.eq.1 bin into the "interior" section 
!                           (as we do for q), since Upd(J1) is undefined
         fij(1) = 0 
         fij(2) = Upd(imem,J1+1,iset) * XV(1,iset)**2 
         fij(3) = Upd(imem,J1+2,iset) * XV(2,iset)**2 
         fij(4) = Upd(imem,J1+3,iset) * XV(3,iset)**2 
!                                                                       
!                 Use CtLhPolint which allows x to be anywhere w.r.t. th
                                                                        
         Call CtLhPolint4(XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
                                                                        
         If (x .GT. 0D0)  fvec(it) =  Fx / x**2 
!                                              Pdf is undefined for x.eq
       ElseIf  (JLx .Eq. Nx(iset)-1) Then 
!                                                This is the highest x b
                                                                        
!** fix allow 4 consecutive elements with iset... mrw 19.9.2005         
        fij(1) = Upd(imem,j1,iset) 
        fij(2) = Upd(imem,j1+1,iset) 
        fij(3) = Upd(imem,j1+2,iset) 
        fij(4) = Upd(imem,j1+3,iset) 
        Call CtLhPolint4 (XVpow(Nx(iset)-3), Fij(1), 4, ss, Fx, Dfx) 
                                                                        
        fvec(it) = Fx 
                                                                        
       Else 
                                                                        
! for all interior points, use Jon's in-line function                   
! This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)                          
! (This is cubic spline interpolation, as used by cteq; it was          
! changed to polint in previous Durham releases (jcp).)                 
         sf2 = Upd(imem,J1+1,iset) 
         sf3 = Upd(imem,J1+2,iset) 
                                                                        
         Fvec(it) = (const5*(Upd(imem,J1,iset)                          &
     &                             - sf2*const1 + sf3*const2)           &
     &             + const6*(Upd(imem,J1+3,iset)                        &
     &                             + sf2*const3 - sf3*const4)           &
     &             + sf2*sy3 - sf3*sy2) / s23                           
                                                                        
       Endif 
                                                                        
      enddo 
!                                   We now have the four values Fvec(1:4
!     interpolate in t...                                               
                                                                        
      If (JLq .LE. 0) Then 
!                         1st Q-bin, as well as extrapolation to lower Q
        Call CtLhPolint4(TV(0,iset), Fvec(1), 4, tt, ff, Dfq) 
                                                                        
      ElseIf (JLq .GE. Nt(iset)-1) Then 
!                         Last Q-bin, as well as extrapolation to higher
        Call CtLhPolint4(TV(Nt(iset)-3,iset), Fvec(1), 4, tt, ff, Dfq) 
      Else 
!                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2) 
!       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
!                         the full range QV(0:Nt)  (in contrast to XV)  
        tf2 = fvec(2) 
        tf3 = fvec(3) 
                                                                        
        g1 = ( tf2*t13 - tf3*t12) / t23 
        g4 = (-tf2*t34 + tf3*t24) / t23 
                                                                        
        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12                      &
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)                       
                                                                        
        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23 
      EndIf 
                                                                        
      CtLhPartonX65 = ff 
                                                                        
      Return 
      END                                           
!=======================================================================
      Function CtLhCtq65Pdf (imem,Iparton, X, Q) 
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
        Print *, 'X out of range in CtLhCtq65Pdf: ', X 
        Stop 
      Endif 
      If (Q .lt. Alambda(iset)) Then 
        Print *, 'Q out of range in CtLhCtq65Pdf: ', Q 
        Stop 
      Endif 
                                                                        
! added to force pdf = 0.0 at x=1.0 exactly - mrw                       
      if(x .eq. 1.0d0) then 
          CtLhCtq65Pdf = 0.0d0 
          return 
      endif 
!                                                                       
      If ((Iparton .lt. -NfMx(iset) .or. Iparton .gt. NfMx(iset))) Then 
         If (Warn) Then 
!        put a warning for calling extra flavor.                        
             Warn = .false. 
             Print *, 'Warning: Iparton out of range in CtLhCtq65Pdf: ' &
     &              , Iparton                                           
         Endif 
         CtLhCtq65Pdf = 0D0 
         Return 
      Endif 
                                                                        
      CtLhCtq65Pdf = CtLhPartonX65 (imem,Iparton, X, Q) 
      if(CtLhCtq65Pdf.lt.0.D0)  CtLhCtq65Pdf = 0.D0 
                                                                        
      Return 
                                                                        
!                             ********************                      
      END                                           
