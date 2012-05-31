! -*- F90 -*-


double precision function alphasPDF(Q)
  implicit none
  integer nset
  double precision Q,a
  call getnset(nset)
  call evolveAs(nset,Q,a)
  alphasPDF=a
  return
end function alphasPDF


double precision function alphasPDFM(nset,Q)
  implicit none
  integer nset
  double precision Q,a
  call evolveAs(nset,Q,a)
  alphasPDFM=a
  return
end function alphasPDFM


subroutine GetQmass(nnf,mass)
  implicit none
  double precision mass
  integer nnf,order,nset
  nset = 1
  call GetQmassM(nset,nnf,mass)
  return
  entry GetOrderAs(order)
  nset = 1
  call GetOrderAsM(nset,order)
  return
end subroutine GetQmass

 
subroutine evolveAs(nset,Q,alphas)
  implicit none
  include 'parmsetup.inc'
  integer iset,imem
  common/SET/iset,imem
  integer nset
  character*16 s1,s2,s3
  double precision Q,alphas,alfas0,scale0
  double precision Q0(nmxset)
  double precision AlfasQ(nmxset)
  double precision b0,b1,b2,L,As,mass
  double precision CtLhALPI,CtLhAlphaNew,CtLhALPInew,CtLhalgluino
  parameter (b0=1.2202,b1=0.4897,b2=0.1913)
  integer n!,k
  integer order
  integer EvlOrd(nmxset)
  integer parm(nmxset)
  integer Etype(nmxset)
  integer Method(nmxset)
  integer nnf,naf
  double precision cmass(nmxset),bmass(nmxset),tmass(nmxset)
  common/masses_LHA/cmass,bmass,tmass
  save EvlOrd,parm,Etype,alfasQ,Q0,Method
  
  call setnset(nset)
  
  if (method(nset).eq.0) then
     if ((Etype(nset).eq.1).or.(Etype(nset).eq.2)) then
        L=log(Q/Q0(nset))
        As=alfasQ(nset)
        if (Etype(nset).eq.2) call GetParmPDF(nset,parm(nset),As)
        if (EvlOrd(nset).eq.0) L=b0*L
        if (EvlOrd(nset).eq.1) L=(b0+As*b1)*L
        if (EvlOrd(nset).eq.2) L=(b0+As*b1+As**2*b2)*L -0.5*As**2*b0*b1/2d0*L**2
        alphas=As/(1.0+As*L)
     endif
  endif
  if (method(nset).eq.1) then
     call alfasevolve(nset,alphas,Q)
  elseif (method(nset).eq.2) then
     call alphamrs(5,alphas,Q)
  elseif (method(nset).eq.3) then
     if(Etype(nset).eq.1) alphas = 3.1415926535898d0*CtLhALPI(Q)
     if(Etype(nset).eq.2) alphas = 3.1415926535898d0*CtLhALPInew(Q)
  elseif (method(nset).eq.4) then
     alphas = CtLhAlphaNew(Q)
  elseif (method(nset).eq.5) then
     call alphamrs(3,alphas,Q)
  elseif (method(nset).eq.6) then
     call alphamrs(4,alphas,Q)
  elseif (method(nset).eq.7) then
     call alphacteq5f34(3,alphas,Q)
  elseif (method(nset).eq.8) then
     call alphacteq5f34(4,alphas,Q)
  elseif (method(nset).eq.9) then
     alphas = CtLhalgluino(Q)
  endif
  return
  
  entry GetQmassM(nset,nnf,mass)
  n=abs(nnf)
  mass=0d0
  if (n.eq.4) mass=cmass(nset)
  if (n.eq.5) mass=bmass(nset)
  if (n.eq.6) mass=tmass(nset)
  return
  
  entry GetNactive(naf,q)
  ! compute nnfn = number of active flavors at scale qin.
  n = 3
  if(q .ge. cmass(nset)) n = 4
  if(q .ge. bmass(nset)) n = n + 1
  if(q .ge. tmass(nset)) n = n + 1
  naf = n
  return
  
  entry GetAlfas(nset,alfas0,scale0)
  scale0=Q0(nset)
  alfas0=alfasQ(nset)
  if (Etype(nset).eq.2) call GetParmPDF(nset,parm(nset),alfas0)
  return
  
  entry GetOrderAsM(nset,order)
  order=EvlOrd(nset)
  return
  
  entry InitAlphasPDF(nset)
  Etype(nset)=-1
  EvlOrd(nset)=-1
  read(1,*) s1,s2,s3
  if (index(s2,'lo').eq.1) EvlOrd(nset)=0
  if (index(s2,'nlo').eq.1) EvlOrd(nset)=1
  if (index(s2,'nnlo').eq.1) EvlOrd(nset)=2
  if (EvlOrd(nset).lt.0) then
     write(*,*) 'File description error:'
     write(*,*) 'Unknown alpha_s evolution order ',s2
     stop
  endif
  if (index(s1,'Fixed').eq.1) then
     Etype(nset)=1
     parm(nset)=-1
     read(1,*) alfasQ(nset),Q0(nset),cmass(nset),bmass(nset),tmass(nset)
  endif
  if (index(s1,'Variable').eq.1) then
     Etype(nset)=2
     alfasQ(nset)=0d0
     read(1,*) parm(nset),Q0(nset),cmass(nset),bmass(nset),tmass(nset)
  endif
  if (Etype(nset).lt.0) then
     write(*,*) 'File description error:'
     write(*,*) 'Unknown alpha_s evolution method ',s1
     stop
  endif
  Method(nset)=-1
  if (index(s3,'Internal').eq.1) Method(nset)=0
  if (index(s3,'EvolCode').eq.1) Method(nset)=1
  if (index(s3,'MRSTalfa').eq.1) Method(nset)=2
  if (index(s3,'CTEQalfa').eq.1) Method(nset)=3
  if (index(s3,'CTEQABalfa').eq.1) Method(nset)=4
  if (index(s3,'MRST3alfa').eq.1) Method(nset)=5
  if (index(s3,'MRST4alfa').eq.1) Method(nset)=6
  if (index(s3,'CTEQ5F3alfa').eq.1) Method(nset)=7
  if (index(s3,'CTEQ5F4alfa').eq.1) Method(nset)=8
  if (index(s3,'CTEQ6LGalfa').eq.1) Method(nset)=9
  if (Method(nset).lt.0) then
     write(*,*) 'File description error:'
     write(*,*) 'Unknown alpha_s method ',s3
     stop
  endif
  return
end subroutine evolveAs


subroutine alphamrs(nflav,alpha,Q)
  IMPLICIT double precision (a-h,o-z)
  include 'parmsetup.inc'
  common/SET/iset,imem
  common/masses_LHA/cmass(nmxset),bmass(nmxset),tmass(nmxset)
  dimension parms(31)
  DATA PI/3.14159/
  DATA TOL/.0005/
  data memold/-1/
  integer nset
  save alambda,memold 
  call getnset(nset)
  
  ! Find alambda corresponding to alphas given in first parameter
  call getnmem(nset,imem)
     call listPDF(nset,imem,parms)
     alphas = parms(1)
     memold=imem
     qz2=8315.
     qz = dsqrt(qz2)
     alambda = 0.3000
     astep = 0.010
     tol2 = 0.0000001
     idir = +1
10   continue
     alambda = alambda + idir*astep
     call mrslambda(nflav,alpha,qz,alambda)
     
     if(idir*(alphas-alpha).gt.0.0) then
        go to 20
     else
        astep = 0.5*astep
        idir = -1*idir
        go to 20
     endif
20   continue
     if(abs(alpha-alphas).gt.tol2) go to 10
  ! alambda found  -- save it !!!!
  ! next call mrslambda to get alphas at q with the correct alambda
  call mrslambda(nflav,alpha,q,alambda)
  RETURN
END subroutine alphamrs


subroutine rgras(alpha,Q2)
  IMPLICIT double precision (a-h,o-z)
  !double precision mc,mb,mt
  include 'parmsetup.inc'
  character*16 name(nmxset)
  integer nmem(nmxset),ndef(nmxset),mem
  common/NAME/name,nmem,ndef,mem
  q=dsqrt(q2)
  call getnset(nset)
  nflav=5
  if(name(nset).eq.'QCDNUM_MRST3') nflav=3     
  if(name(nset).eq.'QCDNUM_MRST4') nflav=4     
  call alphamrs(nflav,alpha,Q)
  return
end subroutine rgras


! new vewrsion of mrslambda 13/5/2004 - includes lo, nlo, and nnlo with flags
subroutine mrslambda(nflav,alpha,Q,alambda)
  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
  DATA PI/3.14159/
  DATA TOL/.0005/
  !  The value of Lambda required corresponds to nflav=4
  !  iord=0 gives leading order, iord=1 gives NLO, iord=2 gives NNLO
  qsdt=8.18    !!  This is the value of 4m_c^2
  qsct=74.0    !!  This is the value of 4m_b^2
  al2=alambda*alambda
  q2=q*q
  t=dlog(q2/al2)
  
  ! CHECK: explicitly initialising ALFQC{3,4,5} (by AB)
  ALFQC3 = 0
  ALFQC4 = 0
  ALFQC5 = 0
  
  call getnset(nset)
  call GetOrderAsM(nset,iord)
  ITH=0
  TT=T
  qsdtt=qsdt/4.
  qsctt=qsct/4.
  AL=ALAMBDA
  AL2=AL*AL
  FLAV=4.
  if(nflav.eq.3) flav=3.
  QS=AL2*dEXP(T)
  
  if(qs.lt.0.5d0) then   !!  running stops below 0.5
     qs=0.5d0
     t=dlog(qs/al2)
     tt=t
  endif
  
  IF(QS.gt.QSCTT.and.nflav.gt.4) GOTO 12  
  IF(QS.lt.QSDTT.and.nflav.gt.3) GOTO 312  
11 CONTINUE
  B0=11-2.*FLAV/3. 
  X1=4.*PI/B0
  IF(IORD.eq.0) then
     ALPHA=X1/T
  ELSE
     if(iord.gt.1) then
        alpha=qwikalf(t,iord,flav)
        go to 51
     endif
     B1=102.-38.*FLAV/3.
     X2=B1/B0**2
     AS2=X1/T*(1.-X2*dLOG(T)/T)
5    AS=AS2
     F=-T+X1/AS-X2*dLOG(X1/AS+X2)
     FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
     AS2=AS-F/FP
     DEL=ABS(F/FP/AS)
     IF((DEL-TOL).GT.0.) go to 5
     ALPHA=AS2
51   continue
  ENDIF
  IF(ITH.EQ.0) RETURN
  GO TO (13,14,15) ITH
  ! GO TO 5
12 ITH=1
  T=dLOG(QSCTT/AL2)
  GO TO 11
13 ALFQC4=ALPHA
  FLAV=5.   
  ITH=2
  
  GO TO 11
14 ALFQC5=ALPHA
  ITH=3
  T=TT
  GO TO 11
15 ALFQS5=ALPHA
  ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
  ALPHA=1./ALFINV
  RETURN
  
311 CONTINUE
  B0=11-2.*FLAV/3. 
  X1=4.*PI/B0
  IF(IORD.eq.0) then
     ALPHA=X1/T
     
  ELSE
     if(iord.gt.1) then
        alpha=qwikalf(t,iord,flav)
        go to 351
     endif
     B1=102.-38.*FLAV/3.
     X2=B1/B0**2
     AS2=X1/T*(1.-X2*dLOG(T)/T)
35   AS=AS2
     F=-T+X1/AS-X2*dLOG(X1/AS+X2)
     FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
     AS2=AS-F/FP
     DEL=ABS(F/FP/AS)
     IF((DEL-TOL).GT.0.) go to 35
     ALPHA=AS2
     
351  continue
  endif
  IF(ITH.EQ.0) RETURN
  GO TO (313,314,315) ITH
312 ITH=1
  T=dLOG(QSDTT/AL2)
  GO TO 311
313 ALFQC4=ALPHA
  FLAV=3.   
  ITH=2
  GO TO 311
314 ALFQC3=ALPHA
  ITH=3
  T=TT
  GO TO 311
315 ALFQS3=ALPHA
  ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
  ALPHA=1./ALFINV
  RETURN
END subroutine mrslambda


double precision function qwikalf(t,iord,flav)
  implicit double precision(a-h,o-z)
  dimension z3(6),z4(6),z5(6),zz3(6),zz4(6),zz5(6)
  data z3/ -.161667E+01,0.954244E+01,0.768623E+01,0.101523E+00,-.360127E-02,0.457867E-04/
  data z4/ -.172239E+01,0.831185E+01,0.721463E+01,0.835531E-01,-.285436E-02,0.349129E-04/
  data z5/ -.872190E+00,0.572816E+01,0.716119E+01,0.195884E-01,-.300199E-03,0.151741E-05/
  data zz3/-.155611E+02,0.168406E+02,0.603014E+01,0.257682E+00,-.970217E-02,0.127628E-03/
  data zz4/-.106762E+02,0.118497E+02,0.664964E+01,0.112996E+00,-.317551E-02,0.302434E-04/
  data zz5/-.531860E+01,0.708503E+01,0.698352E+01,0.274170E-01,-.426894E-03,0.217591E-05/
  
  data pi/3.14159/
  nfm2=flav-2.
  x=dsqrt(t)
  x2=x*x
  x3=x*x2
  x4=x*x3
  x5=x*x4
  go to (1,2) iord
1 go to (3,4,5) nfm2
3 y=z3(1)+z3(2)*x+z3(3)*x2+z3(4)*x3+z3(5)*x4+z3(6)*x5
  go to 10
4 y=z4(1)+z4(2)*x+z4(3)*x2+z4(4)*x3+z4(5)*x4+z4(6)*x5
  go to 10
5 y=z5(1)+z5(2)*x+z5(3)*x2+z5(4)*x3+z5(5)*x4+z5(6)*x5
  go to 10
2 go to (6,7,8) nfm2
6 y=zz3(1)+zz3(2)*x+zz3(3)*x2+zz3(4)*x3+zz3(5)*x4+zz3(6)*x5
  go to 10
7 y=zz4(1)+zz4(2)*x+zz4(3)*x2+zz4(4)*x3+zz4(5)*x4+zz4(6)*x5
  go to 10
8 y=zz5(1)+zz5(2)*x+zz5(3)*x2+zz5(4)*x3+zz5(5)*x4+zz5(6)*x5
  go to 10
10 qwikalf=4.*pi/y
  return
end function qwikalf


!=====================================================================
! alphas routine from PDFLIB
! D OUBLE PRECISION FUNCTION ALPHAS2(SCALE,qcdl5)
subroutine aspdflib(alphas2,SCALE,iord,qcdl5)
  implicit double precision (a-h,o-z)
  double precision NF
  DATA XMC/1.43D0/,XMB/4.30D0/,XMT/100.D0/
  DATA ZEROD/0.D0/,PONED/0.001D0/,ONED/1.D0/,TWOD/2.D0/
  
  LO = 1
  if(iord.ne.0) LO = 2
  
  TMAS = 180.0d0
  
  ALPHAS2 = ZEROD
  PI=4.0D0*ATAN(ONED)
  B6  = (33.D0-2.D0*6.D0)/PI/12.D0
  BP6 = (153.D0 - 19.D0*6.D0) / PI / TWOD / (33.D0 - 2.D0*6.D0)
  B5  = (33.D0-2.D0*5.D0)/PI/12.D0
  BP5 = (153.D0 - 19.D0*5.D0) / PI / TWOD / (33.D0 - 2.D0*5.D0)
  B4  = (33.D0-2.D0*4.D0)/PI/12.D0
  BP4 = (153.D0 - 19.D0*4.D0) / PI / TWOD / (33.D0 - 2.D0*4.D0)
  B3  = (33.D0-2.D0*3.D0)/PI/12.D0
  BP3 = (153.D0 - 19.D0*3.D0) / PI / TWOD / (33.D0 - 2.D0*3.D0)
  XLC = TWOD * LOG( XMC/QCDL5)
  XLB = TWOD * LOG( XMB/QCDL5)
  XLT = TWOD * LOG( XMT/QCDL5 * TMAS/XMT)
  XLLC = LOG( XLC)
  XLLB = LOG( XLB)
  XLLT = LOG( XLT)
  C65  =  ONED/( ONED/(B5 * XLT) - XLLT*BP5/(B5 * XLT)**2 ) - &
       ONED/( ONED/(B6 * XLT) - XLLT*BP6/(B6 * XLT)**2 )
  C45  =  ONED/( ONED/(B5 * XLB) - XLLB*BP5/(B5 * XLB)**2 ) - &
       ONED/( ONED/(B4 * XLB) - XLLB*BP4/(B4 * XLB)**2 )
  C35  =  ONED/( ONED/(B4 * XLC) - XLLC*BP4/(B4 * XLC)**2 ) - &
       ONED/( ONED/(B3 * XLC) - XLLC*BP3/(B3 * XLC)**2 ) + C45
  
  Q = SCALE
  XLQ = TWOD *  LOG( Q/QCDL5 )
  XLLQ =  LOG( XLQ )
  
  ! IF  ( NF .LT. ZEROD) THEN
  IF ( Q .GT. XMT * TMAS/XMT) THEN
     NF = 6.D0
  ELSEIF  ( Q .GT. XMB ) THEN
     NF = 5.D0
  ELSEIF  ( Q .GT. XMC ) THEN
     NF = 4.D0
  ELSE
     NF = 3.D0
  ENDIF
  IF(NF .GT. 6.D0) NF = 6.D0
  IF      ( NF .EQ. 6.D0 ) THEN
     ALF = ONED/(ONED/(ONED/(B6*XLQ)- BP6/(B6*XLQ)**2*XLLQ) + C65)
     IF (LO.EQ.1) ALF = ONED/B6/XLQ
  ELSEIF  ( NF .EQ. 5.D0 ) THEN
     ALF = ONED/(B5 * XLQ) -  BP5/(B5 * XLQ)**2 * XLLQ
     IF (LO.EQ.1) ALF = ONED/B5/XLQ
  ELSEIF  ( NF .EQ. 4.D0 ) THEN
     ALF = ONED/(ONED/(ONED/(B4*XLQ)- BP4/(B4*XLQ)**2*XLLQ) + C45)
     IF (LO.EQ.1) ALF = ONED/B4/XLQ
  ELSEIF  ( NF .EQ. 3.D0 ) THEN
     ALF = ONED/(ONED/(ONED/(B3*XLQ)- BP3/(B3*XLQ)**2*XLLQ) + C35)
     IF (LO.EQ.1) ALF = ONED/B3/XLQ
  ELSE
     WRITE(*,*) 'Error in Alphas2'
     STOP
  ENDIF
  ALPHAS2 = ALF
  RETURN
END subroutine aspdflib


! ========================================================================
SUBROUTINE CtLhAlphaNewSET(MC,MB,MT,Q0,ALPHA0,IORDER,IMODE)
  ! call to set quark masses for alpha_s, and choose lambda or its 
  ! equivalent to make alpha_s take the value alpha0 at scale Q0.
  IMPLICIT NONE
  include 'parmsetup.inc'
  DOUBLE PRECISION MC, MB, MT, Q0, ALPHA0
  INTEGER IORDER, IMODE
  
  DOUBLE PRECISION UDSCBT, QQ0, AALPHA0
  INTEGER IIORDER, IIMODE
  COMMON /QMASSES/ UDSCBT(6), IIMODE, IIORDER
  COMMON /ALSCALE/ QQ0(nmxset), AALPHA0(nmxset)
  integer iset
  double precision Adummy(nmxset)
  integer Nfl(nmxset), Idummy(nmxset)
  common / QCDtable /  Adummy, Nfl, Idummy
  
  call getnset(iset)
  
  if((imode .lt. 1) .or. (imode .gt. 3)) then
     print *,'CtLhAlphaNewSET: fatal imode=',imode
     stop
  endif
  
  IIMODE = IMODE
  QQ0(iset) = Q0
  AALPHA0(iset) = ALPHA0
  IIORDER = IORDER
  
  UDSCBT(1) = .005D0
  UDSCBT(2) = .010D0
  UDSCBT(3) = .300D0
  UDSCBT(4) = MC
  UDSCBT(5) = MB
  UDSCBT(6) = MT
  
  ! set artificial quark masses, if necessary, in alpha_s to enforce 
  ! the requested maximum number of flavors...
  if(Nfl(iset) .le. 5) UDSCBT(6) = 6.d99
  if(Nfl(iset) .le. 4) UDSCBT(5) = 5.d99
  if(Nfl(iset) .le. 3) UDSCBT(4) = 4.d99
  if(Nfl(iset) .le. 2) UDSCBT(3) = 3.d99
  if(Nfl(iset) .le. 1) UDSCBT(2) = 2.d99
  if(Nfl(iset) .le. 0) UDSCBT(1) = 1.d99
  
  RETURN
END SUBROUTINE CtLhAlphaNewSET


FUNCTION CtLhAlphaNew(Q)
  IMPLICIT NONE
  DOUBLE PRECISION CtLhAlphaNew, Q
  DOUBLE PRECISION Q2, CtLhQALPHAS, UDSCBT
  integer nf,ier
  INTEGER IIMODE, IIORDER
  COMMON /QMASSES/ UDSCBT(6), IIMODE, IIORDER
  if ((iimode .ge. 1) .and. (iimode .le. 3)) then
     Q2=Q*Q
     nf=6
     if (Q2.lt. UDSCBT(6)**2) nf=5
     if (Q2.lt. UDSCBT(5)**2) nf=4
     if (Q2.lt. UDSCBT(4)**2) nf=3
     if (Q2.lt. UDSCBT(3)**2) nf=2
     if (Q2.lt. UDSCBT(2)**2) nf=1
     if (Q2.lt. UDSCBT(1)**2) nf=0
     
     ! external maximum number of flavors -- typically, nfmax=5 so 
     ! top quark is not a parton even at large Q...
     CtLhAlphaNew=CtLhQALPHAS(Q2,nf,ier)
     if(ier .ne. 0) then
        print *,'warning in CtLhAlphaNew, Q=',Q,' nf=',nf,' ier=',ier,' CtLhAlphaNew=',CtLhAlphaNew
     endif
  else
     print *,'CtLhAlphaNew: undefined mode=',iimode 
     stop
  endif
  
  return
end FUNCTION CtLhAlphaNew


FUNCTION CtLhQALPHAS(QQ2,NF,IERR)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  include 'parmsetup.inc' 
  COMMON /ALSCALE/ Q0(nmxset), ALPHA0(nmxset)
  COMMON /QMASSES/ UDSCBT(6), IIMODE, IIORDER
  
  call getnset(iset)
  
  Q02 = Q0(iset)**2
  ALP0 = ALPHA0(iset)
  IOR = IIORDER
  
  NFFF = NF
  CtLhQALPHAS = CtLhA0TOA1(QQ2,Q02,ALP0,IOR,NFFF,IERR)
  
  RETURN
END FUNCTION CtLhQALPHAS


FUNCTION CtLhA0TOA1(QSU,QS0,AS0,IORD,NFF,IERR)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  COMMON /QMASSES/ UDSCBT(6), IIMODE, IIORDER
  QS1   = QSU
  QMU0  = SQRT(QS0)
  QMU1  = SQRT(QS1)
  
  ! CHECK: explicit initialisation of NF{0,1} (by AB)
  NF0 = 0
  NF1 = 0
  DO I = 1, 6
     IF(QMU0.GE.UDSCBT(I)) NF0 = I
     IF(QMU1.GE.UDSCBT(I)) NF1 = I
  ENDDO
  
  IF(NF1.LT.NF0) THEN
     IST = -1
     JST =  0
  ELSE
     IST = 1
     JST = 1
  ENDIF
  
  ALFA0 = AS0
  Q00   = QS0
  
  IERR = 0
  DO NF = NF0,NF1,IST
     IF(NF.NE.NF1) THEN
        Q21 = UDSCBT(NF+JST)**2
     ELSE
        Q21 = QS1
     ENDIF
     IMODE = IIMODE
     ALFA1 = CtLhALPHAR(Q21,Q00,ALFA0,NF,IORD,IMODE,JERR)
     IERR = IERR + JERR!IERR is sum of all errors
     ALFA0 = ALFA1
     Q00   = Q21
  END DO
  
  CtLhA0TOA1 = ALFA0
  NFF    = NF1
  
  RETURN
END FUNCTION CtLhA0TOA1


FUNCTION CtLhALPHAR(QSQ,QS0,AS0,NF,IORD,IMODE,IERR)
  ! calculate ALPHAS FROM RGE GIVEN AS0 AT QS0.
  !
  ! IORD=1: LEADING ORDER DEFINED BY 
  !                 Q*d(alpha)/d(Q) = c1*alpha**2 
  !
  ! IORD=2,IMODE=1: QCD NUM CHOICE DEFINED BY 
  !                 Q*d(alpha)/d(Q) = c1*alpha**2 + c2*alpha**3
  !
  ! IORD=2,IMODE=2: AD HOC ALTERNATIVE DEFINED BY 
  !                 Q*d(alpha)/d(Q) = c1*alpha**2 / (1 - (c2/c1)*alpha)
  !
  ! IORD=2,IMODE=3: TRADITIONAL CTEQ CHOICE DEFINED BY 
  !                 ALPHA = c3*(1 - c4*log(L)/L)/L, WHERE L=log((Q/lambda)**2)
  !
  ! c1 = -beta0/(2*pi)     where beta0 =  11. -  (2./3.)*nf
  ! c2 = -beta1/(8*pi**2)  where beta1 = 102. - (38./3.)*nf
  !
  ! c3 = -2/c1
  ! c4 = -2*c2/c1**2
  !
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  DATA PI / 3.14159265358979d0 /
  
  BET0 = 11.d0 -(2*NF)/3.d0
  BET1 = 102.d0 - (38*NF)/3.d0
  B0   = BET0/(4.d0*PI)
  B1   = BET1/(4.d0*PI*BET0)
  IERR = 0
  
  TERM0 = 1.d0/AS0+B0*LOG(QSQ/QS0)
  IF (TERM0.LE.0.) THEN
     CtLhALPHAR = 100.
     IERR   = 1
     PRINT *,'CtLhALPHAR WARNING: RETURN 100.'
     RETURN
  ENDIF
  ALFA0 = 1.d0/TERM0

  ! ORDER=1 IS LEADING ORDER, WHICH IS SAME FOR ALL IMODE.
  IF (IORD.EQ.1) THEN
     CtLhALPHAR = ALFA0
     RETURN
  ELSEIF(IORD.NE.2) THEN
     PRINT *,'FATAL ERROR: UNDEFINED ORDER IN CtLhALPHAR'
     STOP
  ENDIF

  
  ! QCD NUM CHOICE: Q*d(alpha)/d(Q) = c1*alpha**2 + c2*alpha**3
  IF(IMODE .EQ. 1) THEN
     
     ! use Newton's method to solve the equation, instead of the 
     ! simple iterative method used in qcdnum (jcp 9/01)
     DO ITER = 1, 20
        ARG   = (1.d0/ALFA0+B1)/(1.d0/AS0+B1)
        IF(ARG.LE.0.) THEN
           CtLhALPHAR = 10.
           IERR   = 1
           PRINT *,'CtLhALPHAR WARNING: RETURN 10.'
           RETURN
        ENDIF
        
        TERM  = TERM0 + B1*LOG(ARG) - 1.d0/ALFA0
        ALFA1 = ALFA0/(1.d0 + ALFA0*(1.d0 + B1*ALFA0)*TERM)
        
        IF(ABS(ALFA1-ALFA0).LT.1.E-12) GOTO 20
        ALFA0 = ALFA1
     ENDDO
     
     CtLhALPHAR = 10.
     IERR   = 1
     RETURN
     
20   CONTINUE
     CtLhALPHAR = ALFA1
     RETURN

     ! AD HOC ALTERNATIVE: Q*d(alpha)/d(Q) = c1*alpha**2 / (1 - (c2/c1)*alpha)
  ELSEIF(IMODE .EQ. 2) THEN
     
     ! first get a good starting point, to be sure Newton's method doesn't go 
     ! to the wrong root.
     BEST = 9.d99
     ! CHECK: adding explicit initialisation of ALBST (by AB)
     ALBST = ALFA0
     DO ITRY = 0, 20
        IF(ITRY.NE.0) ALFA0 = (1./B1)*(ITRY-0.5D0)/20.D0
        F = -1.d0/ALFA0 + TERM0 - B1*LOG(ALFA0/AS0) 
        IF(ABS(F) .LT. BEST) THEN
           BEST = ABS(F)
           ALBST = ALFA0
        ENDIF
     ENDDO
     
     ALFA0 = ALBST
     DO ITER=1, 20
        F = -1.d0/ALFA0 + TERM0 - B1*LOG(ALFA0/AS0) 
        ALFA1 = ALFA0/(1.d0 + ALFA0*F/(1.d0 - B1*ALFA0))
        IF(ABS(ALFA1-ALFA0) .LT. 1.E-12) GOTO 30
        ALFA0 = ALFA1
     ENDDO
     
     CtLhALPHAR = 10.
     IERR   = 1
     RETURN
     
30   CONTINUE
     CtLhALPHAR = ALFA1
     RETURN
     
     ! TRADITIONAL CTEQ CHOICE: ALPHA = c3*(1 - c4*log(L)/L)/L, WHERE L=log((Q/lambda)**2)
  ELSEIF(IMODE .EQ. 3) THEN
     
     Z = -LOG(B0*AS0)
     TMP = BET1/BET0**2
     
     DO ITER = 1, 20
        F = EXP(Z) - (1.D0 - TMP*Z*EXP(-Z))/(B0*AS0)
        FPRI = EXP(Z)  + TMP*(1.D0-Z)*EXP(-Z)/(B0*AS0)
        ZNEW = Z - F/FPRI
        IF(ABS(Z-ZNEW) .LT. 1.E-10) GOTO 40
        Z = ZNEW
     ENDDO
     
     CtLhALPHAR = 10.
     IERR   = 1
     RETURN
     
40   CONTINUE
     XLAMSQ = QS0 * EXP(-EXP(ZNEW))
     
     XL = LOG(QSQ/XLAMSQ)
     
     ! return a fixed value if no solution...
     IF (XL .LE. 0.D0) THEN
        CtLhALPHAR = 10.D0
        IERR = 1
        RETURN
     ENDIF
     
     CtLhALPHAR = (1.d0 - TMP*LOG(XL)/XL)/(B0*XL)
     
     ! place a cutoff if comes out very large...
     if (CtLhALPHAR .gt. 10.d0) then
        CtLhALPHAR = 10.d0
        IERR = 1
     endif
     
     RETURN
     
  ELSE
     PRINT *,'FATAL UNDEFINED IMODE=',IMODE
     STOP
  ENDIF
  
  RETURN
END FUNCTION CtLhALPHAR


FUNCTION CtLhALPInew (AMU)
  ! Returns effective g**2/(4pi**2) = alpha/pi.
  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
  data pi / 3.14159265358979d0 /
  q = amu
  alpha = CtLhAlphaNew(q)
  CtLhalpinew = alpha/pi
  RETURN
END FUNCTION CtLhALPInew


subroutine CTEQ6NewAlpha(nset,mem)
  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
  include 'parmsetup.inc'
  common / QCDtable /  Alambda(nmxset), Nfl(nmxset), Iorder(nmxset)
  dimension parms(2)
  q0 = 91.188d0
  call listPDF(nset,mem,parms)
  alpha0 = parms(1)
  ximode = parms(2)
  imode = int(ximode)
  iiorder = iorder(nset)
  
  ! **********************************************************
  ! the input data file should probably specify a very large 
  ! value for xmt, because we never allow top quark as a parton
  ! in the PDF fitting, so there are never more than 5 active 
  ! flavors.  Presume it is most consistent therefore to 
  ! also keep top quark loop out of running of alpha.
  ! **********************************************************
  
  call GetQmassM(nset,4,cmass)
  call GetQmassM(nset,5,bmass)
  call GetQmassM(nset,6,tmass)
  call CtLhAlphaNewSET(cmass,bmass,tmass,Q0,ALPHA0,IIORDER,IMODE)
  RETURN
END subroutine CTEQ6NewAlpha


subroutine getnset(nset)
  integer iset,nset
  save iset
  nset = iset      
  return
  
  entry setnset(nset)
  iset = nset
  return      
end subroutine getnset


subroutine getnmem(nset,nmem)
  include 'parmsetup.inc'
  integer nmem,nset,member(nmxset)
  save member
  nmem = member(nset)      
  return
  
  entry setnmem(nset,nmem)
  member(nset) = nmem
  return      
end subroutine getnmem


subroutine alphacteq5f34(nflav,alphas,q)
  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
  parameter(pi=3.14159265358979323846d0)
  
  alphas = pi*CtLhALPI34(nflav,q)
  return
end subroutine alphacteq5f34


FUNCTION CtLhALPI34 (nflav,AMU)
  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
  COMMON / LhCtCWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
  COMMON / LhCtQCDPAR_LHA / AL, NF, NORDER, SET
  LOGICAL SET
  PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
  DATA IW1, IW2 / 2*0 /
  IF(.NOT.SET) CALL CtLhLAMCWZ
  NEFF = LhCtNFL(AMU)
  
  if(neff.gt.nflav) neff=nflav
  
  ALM  = ALAM(NEFF)
  if(neff.eq.3) alm = 0.395
  if(neff.eq.4) alm = 0.309
  CtLhALPI34 = CtLhALPQCD (NORDER, NEFF, AMU/ALM, IRT)
  IF (IRT .EQ. 1) THEN
     CALL CtLhWARNR (IW1, 'AMU < ALAM in CtLhALPI34', 'AMU', AMU,ALM, BIG, 1)
  ELSEIF (IRT .EQ. 2) THEN
     CALL CtLhWARNR(IW2,'CtLhALPI34 > 3; Be aware!','CtLhALPI34',CtLhALPI34, D0, D1, 0)
  ENDIF
  RETURN
END FUNCTION CtLhALPI34


function CtLhalgluino(Q)
  !     51  mg=50  GeV  WITH TOP MASS OF 174
  !     52  mg=15  GeV  
  !     53  mg=25  GeV  
  !     54  mg=35  GeV  
  !     55  mg=50  GeV  
  !     56  mg=90  GeV  
  !     57  mg=100  GeV  
  !     58  mg=200  GeV  
  !     59  mg=inf  GeV  

  IMPLICIT DOUBLE PRECISION (A-H, O-Z)
  Parameter (nsets=8)
  Dimension xmass(6), gArray(nsets), nfArray(nsets)
  Dimension xlamStrange(nsets),xlamCharm(nsets),xlamBottom(nsets),xlamGluino(nsets),xlamTop(nsets)
  data iorder /2/
  Data xmass / 0.0, 0.0, 0.5d0, 1.3d0, 4.5d0, 174.0d0/
  Data gArray / 15.0d0, 25.0d0, 35.0d0, 50.0d0, 90.0d0, 100.0d0, 200.0d0, 1.0d99/  !*** Last value = Infinity
  Data  xlamStrange /  0.240257, 0.272782, 0.296194, 0.322873,  0.371374,  0.372537, 0.372537, 0.372537/
  Data  xlamCharm   /  0.197153, 0.227834, 0.250338, 0.276386, 0.324772, 0.325947, 0.325947, 0.325947/
  Data  xlamBottom  /  0.128981, 0.151553, 0.168364, 0.188074,  0.225349, 0.226264, 0.226264, 0.226264/
  Data  xlamGluino  / 0.005,  0.005,  0.005,  0.005,  0.005, 0.00480354,  0.00355387,  0.0/
  
  Data  xlamTop     /  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
  Data nfArray /5,5,5,5,5,5,5,5/ 
  
  ! 3456789012345678901234567890123456789012345678901234567890123456789012
  
  zmass=91.2
  
  pi=4.0d0*atan(1.0d0)
  
  call getnset(iset)
  call getnmem(iset,imem)
  
  nfMax=nfArray(imem)
  nf=3
  if(Q.ge.xmass(4))  nf=4
  if(Q.ge.xmass(5))  nf=5
  if(Q.ge.xmass(6))  nf=6  !*** Top for Iset=51 only
  if(nf.gt.nfMax) nf=nfMax
  
  gmass=gArray(imem+1)
  ng=0
  if(Q.ge.gmass) ng=1
  
  xlam = xlamStrange(imem+1)
  if(Q.ge.xmass(4))  xlam = xlamCharm(imem+1)
  if(Q.ge.xmass(5))  xlam = xlamBottom(imem+1)
  if(Q.ge.gmass   )  xlam = xlamGluino(imem+1)
  
  if((Q.ge.xmass(6)).and.(nfMax.ge.6))  xlam = xlamTop(imem+1)
  
  nsquark=0
  b0=11. - 2./3.*nf - 2.*ng- 1./6.*nsquark
  b1=102.-38./3.*nf-48.*ng-11./3.*nsquark+13./3.*ng*nsquark
  if(iorder.eq.1) b1=0
  
  xlog=Log(Q**2/xlam**2)
  
  tmp= 4.*pi*( 1./(b0*xlog) - b1*Log(xlog)/(b0**3*xlog**2))
  
  CtLhalgluino=tmp
  Return
End function CtLhalgluino
