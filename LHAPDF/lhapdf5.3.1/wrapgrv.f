      subroutine GRVevolve(xin,qin,pdf)
      implicit real*8 (a-h,o-z)
      include 'parmsetup.inc'
      PARAMETER(ngrid=2)
      PARAMETER (NPART=6, NX=68, NQ=27, NARG=2)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      DIMENSION XXUVF(0:ngrid,NX,NQ), XXDVF(0:ngrid,NX,NQ), 
     +          XXDEF(0:ngrid,NX,NQ), XXUDF(0:ngrid,NX,NQ),
     1          XXSF(0:ngrid,NX,NQ),   XXGF(0:ngrid,NX,NQ),
     +          XUVF(NX,NQ), XDVF(NX,NQ), 
     +          XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ),   XGF(NX,NQ),
     +          PARTON (NPART,NQ,NX-1), 
     2          QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      CHARACTER*80 LINE
      dimension pdf(-6:6)
      save 
      x=xin
      q2=qin*qin      
*...CHECK OF X AND Q2 VALUES : 
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
         WRITE(6,91) 
  91     FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
         STOP
      ENDIF
      IF ( (Q2.LT.0.799) .OR. (Q2.GT.1.01E6) ) THEN
         WRITE(6,92) 
  92     FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
         STOP
      ENDIF
*...INTERPOLATION :
      DO IQ=1,NQ
      DO IX=1,NX
      xuvf(ix,iq)=xxuvf(imem,ix,iq)
      xdvf(ix,iq)=xxdvf(imem,ix,iq)
      xdef(ix,iq)=xxdef(imem,ix,iq)
      xudf(ix,iq)=xxudf(imem,ix,iq)
      xsf(ix,iq)=xxsf(imem,ix,iq)
      xgf(ix,iq)=xxgf(imem,ix,iq)
      enddo
      enddo
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      X1 = 1.- X
      XV = X**0.5
      XS = X**(-0.2)
      UV = FINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = FINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = FINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = FINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      US = 0.5 * (UD - DE)
      DS = 0.5 * (UD + DE)
      SS = FINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = FINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS 
*

      pdf(-6) = 0.0d0
       pdf(6) = 0.0d0
      pdf(-5) = 0.0d0
       pdf(5) = 0.0d0
      pdf(-4) = 0.0d0
       pdf(4) = 0.0d0
      pdf(-3) = ss
       pdf(3) = ss
      pdf(-2) = us
       pdf(2) = uv+us
      pdf(-1) = ds
       pdf(1) = dv+ds
       pdf(0) = gl
      return
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVread(nset)
c      
c      print *,'calling grvread'
      read(1,*)nmem(nset),ndef(nset)

      read(1,93)xb
c      print *,xb
      read(1,93)qs
c      print *,qs
  93  format(8e8.2)
c
      do ng=0,nmem(nset)
c      
      READ(1,89) LINE
  89  FORMAT(A80)
c      print *,line
      DO 15 M = 1, NX-1 
      DO 15 N = 1, NQ
      READ(1,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1           PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M) 
  90  FORMAT (6(1PE10.3))
  15  CONTINUE     
*
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0V = XB(IX)**0.5 
        XB0S = XB(IX)**(-0.2) 
        XB1 = 1.-XB(IX)
        xXUVF(ng,IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0V)
        xXDVF(ng,IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0V)
        xXDEF(ng,IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0V) 
        xXUDF(ng,IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0S)
        xXSF(ng,IX,IQ)  = PARTON(5,IQ,IX) / (XB1**7 * XB0S)
        xXGF(ng,IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0S)
  20  CONTINUE
        xXUVF(ng,NX,IQ) = 0.E0
        xXDVF(ng,NX,IQ) = 0.E0
        xXDEF(ng,NX,IQ) = 0.E0
        xXUDF(ng,NX,IQ) = 0.E0
        xXSF(ng,NX,IQ)  = 0.E0
        xXGF(ng,NX,IQ)  = 0.E0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
      
      enddo

      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRValfa(alfas,qalfa)
      call getnset(iset)
      q2alfa = qalfa*qalfa
      call GetOrderAsM(iset,iord)
      nord=iord+1
      alfas=grvals(q2alfa,nord)
      alfas = 4.0d0*3.14159d0*alfas
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry GRVpdf(mem)
      imem = mem
      return
c
 1000 format(5e13.5)
      end
c
      FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE)
*********************************************************************
*                                                                   *
*   THE INTERPOLATION ROUTINE (CERN LIBRARY ROUTINE E104)           *
*                                                                   *
*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ARG(5),NENT(5),ENT(10),TABLE(10)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT=0.
   10 FAC=1.
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT=FINT+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
           
c      FUNCTION ALPHAS (Q2, NAORD)
      FUNCTION grvals (Q2, NAORD)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE.                                            *
*                                                                   *
*   INPUT :  Q2    =  scale in GeV**2  (not too low, of course);    *
*            NAORD =  1 (LO),  2 (NLO).                             *
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) for use with the GRV(98) partons.      *  
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A - Z)
      INTEGER NF, K, I, NAORD
      DIMENSION LAMBDAL (3:6),  LAMBDAN (3:6), Q2THR (3)
*
*...HEAVY QUARK THRESHOLDS AND LAMBDA VALUES :
      DATA Q2THR   /  1.960,  20.25,  30625. /
      DATA LAMBDAL / 0.2041, 0.1750, 0.1320, 0.0665 /
      DATA LAMBDAN / 0.2994, 0.2460, 0.1677, 0.0678 /
*
*...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
      NF = 3
      DO 10 K = 1, 3
      IF (Q2 .GT. Q2THR (K)) THEN
         NF = NF + 1
      ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
         ALP  = 1./(B0 * DLOG (Q2/LAM2))
         GO TO 1
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
         WRITE (6,91)
  91     FORMAT ('INVALID CHOICE FOR ORDER IN ALPHA_S')
         STOP
       END IF
*
*...START VALUE FOR NLO ITERATION :
       LQ2 = DLOG (Q2 / LAM2)
       ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)
*
*...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE :
       DO 2 I = 1, 6
       XL  = DLOG (1./(B0*ALP) + B10)
       XLP = DLOG (1./(B0*ALP*1.01) + B10)
       XLM = DLOG (1./(B0*ALP*0.99) + B10)
       Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
       Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP
     1       + 1./ (B0*ALP*0.99) - B10 * XLP) / (0.02D0*ALP)
       ALP = ALP - Y/Y1
  2    CONTINUE
*
*...OUTPUT :
c  1    ALPHAS = ALP
  1    grvals = ALP
       RETURN
       END


 
