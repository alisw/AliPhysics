! -*- F90 -*-


!****************************************************************************
!
!                         eps08.f
!
! An interface for calculating the scale dependent nuclear modifications
!                 R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 
! where f_A is the distribution of the parton flavour f for a PROTON in a
! nucleus A, and f_p is the corresponding parton distribution in the 
! free proton.
!
! The EPS08 set of modifications R_f^A(x,Q) which we provide here,
! corresponds to the LO DGLAP evolution of the nPDFs f_A(x,Q).
!
! When using this interface, please refer to:
!  
!  K.J. Eskola, H. Paukkunen and C.A. Salgado,
!  "An improved global analysis of nuclear parton distribution functions
!  including RHIC data"
!  arXiv:0802.0139 [hep-ph].
!
! Questions & comments to:
!           hannu.paukkunen@phys.jyu.fi
!           carlos.salgado@cern.ch
!           kari.eskola@phys.jyu.fi
! 
! ***************************************************************************
! Instructions:
!
! To obtain the nuclear corrections for given input
! values of momentum fraction x, scale Q (in GeV),
! and atomic number A (all in double precision), use the command
!
!    Call eps08(x,Q,A,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
!
! which returns the bound proton nuclear corrections R_f^A(x,Q)
! (in double precision) for
!        
!        ruv = up valence
!        rdv = down valence
!        ru  = up sea
!        rd  = down sea
!        rs  = strange
!        rc  = charm
!        rb  = bottom (=0 for Q < mb=4.64)
!       rt  = top (=1 for all x and Q)
!        rg  = gluons
!
! The nuclear corrections for bound neutrons can be obtained
! by the isospin symmetry (see Sec.2 of article quoted above),
! so that e.g. the total up quark distribution per nucleon in a
! nucleus A with Z protons is
!
!  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
!            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]
!
! Note that the parametrization should only be applied at the region
!
!             1e-6 <= x < 1
!              1.3 <= Q <= 1000 GeV
!                2 <= A <= 208
!
! No warning message is displayed if these limits are
! exceeded, and outside these boundaries the modifications
! are frozen to the boundary values, i.e
!
!   for Q > 1000, the modifications at Q=1000 are returned,
!   for Q < 1.3,  the modifications at Q=1.3 are returned,
!   for x < 1e-6, the modifications at x=1e-6 are returned,
!   for A > 208,  the modifications at A=208 are returned,
!   for A <= 2 or x >= 1, the function returns 1.
!
! The program uses the data file eps08dta.all,
! which must be located in current working directory.
!
! *********************************************************
! *********************************************************

subroutine eps08(xxx,q,aaa,ruv,rdv,ru,rd,rs,rc,rb,rt,rg)
  
  implicit double precision(a-h,o-z)
  dimension rateps08(8)
  dimension pint08(5), pint082(5)
  common/eps0800/qq0, qqb
  common/eps0803/ kty, ind
  common/eps0810/qq1
  common/eps0811/pp08(5)
  external eps08gv, eps08se
  data readFR/0/
  
  aa = aaa
  xx = xxx
  qq = q**2
  rt = 1.d0   
  
  if (readFR.ne.1) then
     call eps08init
     readFR=1
  endif
  if ((aa.le.2.d0).or.(1.d0.le.xx)) then
     ruv=1.d0
     rdv=1.d0
     ru=1.d0
     rd=1.d0
     rs=1.d0
     rc=1.d0
     rb=1.d0
     rg=1.d0
     return
  endif
  
  if (qq .lt. qq0) then
     qq = qq0
  end if
  
  if (qq .gt. 1.d6) then
     qq = 1.d6
  end if
  
  if (aa .gt. 208.d0) then
     aa = 208.d0
  end if
  
  if (xx .lt. 1.d-6) then
     xx = 1.d-6
  end if
  
  if (xx .gt. 0.964) then
     xx = 0.964
  end if
  
  x2=xx+10.d0**(dlog10(xx)-2.5)
  ind=1+int(20.d0*(dlog10(x2)+6.d0))
  x0=10.d0**(-6.d0+(ind-1)/20.d0)
  x1=10.d0**(-6.d0+ind/20.d0)
  if(x2.gt.0.1d0) then
     ind=int((x2-0.1d0)/0.009d0+101.d0)
     x0=0.1d0+(ind-101)*0.009d0
     x1=0.1d0+(ind-100)*0.009d0
  endif
  do ikty=1,3
     kty=ikty
     if(ikty.eq.3) kty=8
     qq1=qq0 
     ! The xx arg was not used anywhere!
     !call ppeps08(xx,aa,pint08,pint082) 
     call ppeps08(aa,pint08,pint082) 
     do ii=1,5
        pp08(ii)=pint08(ii)
     enddo
     tt1=eps08gv(qq)
     do ii=1,5
        pp08(ii)=pint082(ii)
     enddo
     tt2=eps08gv(qq)
     rateps08(kty)=tt1+(tt2-tt1)/(x1-x0)*(xx-x0)
  enddo
  do ikty=3,7
     kty=ikty
     qq1=qq0 
     
     if(kty.eq.7) then
        qq1=qqb 
        if(qq .lt. qqb) then 
           rateps08(7)=0.d0
           Goto 125
        end if
     end if
     
     ! The xx arg was not used anywhere!
     !call ppeps08(xx,aa,pint08,pint082) 
     call ppeps08(aa,pint08,pint082) 
     do ii=1,5
        pp08(ii)=pint08(ii)
     enddo
     tt1=eps08se(qq)
     do ii=1,5
        pp08(ii)=pint082(ii)
     enddo
     tt2=eps08se(qq)
     rateps08(ikty)=tt1+(tt2-tt1)/(x1-x0)*(xx-x0)
  enddo
  
125 continue  
  
  ruv=rateps08(1)
  rdv=rateps08(2)
  rs=rateps08(3)
  ru=rateps08(4)
  rd=rateps08(5)
  rc=rateps08(6)
  rb=rateps08(7)
  rg=rateps08(8)
  return
end subroutine eps08

subroutine ppeps08(aa,pint08,pint082)
  !subroutine ppeps08(xx,aa,pint08,pint082)
  ! The xx arg was not used anywhere!
  implicit double precision(a-h,o-z)
  dimension pint08(5), pint082(5), paint081(13), paint082(13) 
  common/eps0801/ peps08i(197,13,8,5)
  common/eps0802/ amx(13)
  common/eps0803/ kty, ind
  do ip=1,5
     do ia=1,13
        paint081(ia)=peps08i(ind,ia,kty,ip)
        paint082(ia)=peps08i(ind+1,ia,kty,ip)
        if(ind.eq.197) paint082(ia)=paint081(ia)
     enddo
     pint08(ip)=eps08ddv(paint081,amx,13,aa,1)
     pint082(ip)=eps08ddv(paint082,amx,13,aa,1)
  enddo
end subroutine ppeps08


function eps08gv(qq)
  implicit double precision(a-h,o-z)
  common/eps0810/qq1
  common/eps0811/yy1,p1,p2,p3,p4
  z=dlog(qq/qq1)
  zz=dlog(1.d0+z)
  eps08gv=yy1+p1*zz+p2*zz**2+p3*zz**0.5+p4/(1.d0+zz)**4-p4
  return
end function eps08gv

function eps08se(qq)
  implicit double precision(a-h,o-z)
  common/eps0810/qq1
  common/eps0811/yy1,p1,p2,p3,p4
  z=dlog(qq/qq1)
  zz=dlog(1.d0+z)
  eps08se=yy1+p1*zz+p2*zz**2+p3*zz**3+p4/(1.d0+zz)-p4
  return
end function eps08se

subroutine eps08init
  implicit double precision(a-h,o-z)
  include 'commonlhapdfc.inc'
  include 'commonlhacontrol.inc'
  !
  ! the LHAPDF bits
  !
  character*512 dirpath,setpath
  common/eps0800/qq0, qqb
  common/eps0801/peps08i(197,13,8,5)
  common/eps0802/ amx(13)
  data readFR2/0/
  
  ! Initialise common blocks  
  call commoninit()

  qq0 = 1.69d0
  qqb = 21.516d0
  
  amx(1) = 4.d0
  amx(2) = 6.d0
  amx(3) = 9.d0
  amx(4) = 12.d0
  amx(5) = 27.d0
  amx(6) = 40.d0
  amx(7) = 56.d0
  amx(8) = 64.d0
  amx(9) = 108.d0
  amx(10) = 117.d0
  amx(11) = 184.d0
  amx(12) = 197.d0
  amx(13) = 208.d0
  
  if (readFR2.ne.1) then
     !
     ! sort out the path to the input file
     !
     call getdirpath(dirpath)
     
     setpath=dirpath(:LEN_TRIM(dirpath))//"/"//'eps08dta.LHgrid'
     open(11,file=setpath(:LEN_TRIM(setpath)),status='unknown')
     do ik=1,8
        do ia=1,13
           do ix=1,197
              read(11,137) peps08i(ix,ia,ik,1), peps08i(ix,ia,ik,2), &
                   peps08i(ix,ia,ik,3), peps08i(ix,ia,ik,4), peps08i(ix,ia,ik,5)
           enddo
        enddo
     enddo
     close(11)
     readFR2=1
  endif
  return
137 format(5f15.8)
end subroutine eps08init


!        divdif routine from CERNLIB adapted to double precision
!
! $Id: divdif.F,v 1.1.1.1 1996/02/15 17:48:36 mclareni Exp $
!
! $Log: divdif.F,v $
! Revision 1.1.1.1  1996/02/15 17:48:36  mclareni
! Kernlib
!
!
FUNCTION EPS08ddv(F,A,NN,X,MM)
  implicit double precision (a-h,o-z)
  DIMENSION A(NN),F(NN),T(20),D(20)
  LOGICAL EXTRA
  LOGICAL MFLAG !,RFLAG
  DATA MMAX/10/
  !
  !  TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
  !
  !  START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
  !carlos      IF( (NN.LT.2) .OR. (MM.LT.1) ) GO TO 20
  N=NN
  M=MIN0(MM,MMAX,N-1)
  MPLUS=M+1
  IX=0
  IY=N+1
  IF(A(1).GT.A(N)) GO TO 4
  !     (SEARCH INCREASING ARGUMENTS.)
1 MID=(IX+IY)/2
  IF(X.GE.A(MID)) GO TO 2
  IY=MID
  GO TO 3
  !        (IF TRUE.)
2 IX=MID
3 IF(IY-IX.GT.1) GO TO 1
  GO TO 7
  !     (SEARCH DECREASING ARGUMENTS.)
4 MID=(IX+IY)/2
  IF(X.LE.A(MID)) GO TO 5
  IY=MID
  GO TO 6
  !        (IF TRUE.)
5 IX=MID
6 IF(IY-IX.GT.1) GO TO 4
  !
  !  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
  !  *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
7 NPTS=M+2-MOD(M,2)
  IP=0
  L=0
  GO TO 9
8 L=-L
  IF(L.GE.0) L=L+1
9 ISUB=IX+L
  IF((1.LE.ISUB).AND.(ISUB.LE.N)) GO TO 10
  !        (SKIP POINT.)
  NPTS=MPLUS
  GO TO 11
  !        (INSERT POINT.)
10 IP=IP+1
  T(IP)=A(ISUB)
  D(IP)=F(ISUB)
11 IF(IP.LT.NPTS) GO TO 8
  EXTRA=NPTS.NE.MPLUS
  !
  !  REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, SUP-
  !  PLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
  DO L=1,M
     IF(.NOT.EXTRA) GO TO 12
     ISUB=MPLUS-L
     D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
12   I=MPLUS
     DO J=L,M
        ISUB=I-L
        D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
        I=I-1
     ENDDO
  ENDDO
  !
  !  EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
  !  OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
  SUM=D(MPLUS)
  IF(EXTRA) SUM=0.5*(SUM+D(M+2))
  J=M
  DO L=1,M
     SUM=D(J)+(X-T(J))*SUM
     J=J-1
  ENDDO
  EPS08ddv=SUM
  RETURN
  !
  !carlos   20 CALL KERMTR('E105.1',LGFILE,MFLAG,RFLAG)
  EPS08ddv=0
  IF(MFLAG) THEN
     IF(LGFILE.EQ.0) THEN
        IF(MM.LT.1) WRITE(*,101) MM
        IF(NN.LT.2) WRITE(*,102) NN
     ELSE
        IF(MM.LT.1) WRITE(LGFILE,101) MM
        IF(NN.LT.2) WRITE(LGFILE,102) NN
     ENDIF
  ENDIF
  !carlos      IF(.NOT.RFLAG) CALL ABEND
  RETURN
101 FORMAT( 7X, 'FUNCTION EPS08ddv ... M =',I6,' IS LESS THAN 1')
102 FORMAT( 7X, 'FUNCTION EPS08ddv ... N =',I6,' IS LESS THAN 2')
END FUNCTION EPS08ddv
