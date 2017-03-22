*CMZ :          26/11/93  11.31.36  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE DELSLV(N,FMAJOR,FMINOR,LMAX,FRACT,X,XLOW,XUP,VOL,NCUT,N
     1CDIM,ICUT,DELPLS,DELNEG,REGTOL,FTOL,FORIG,DFORIG,FNEW,FNLIN,FNLROW
     2,DIAGJ,SPDIAG,SOL,Z,NFCNT)
      INTEGER N, NCUT, NCDIM, NFCNT
      INTEGER ICUT(NCDIM)
      DOUBLE PRECISION FMAJOR, FMINOR, FRACT, REGTOL, FTOL
      DOUBLE PRECISION X(N), XUP(N), XLOW(N), DELPLS(N), DELNEG(N)
      DOUBLE PRECISION FORIG(NCDIM), FNLIN(NCDIM)
      DOUBLE PRECISION DFORIG(NCDIM), FNEW(NCDIM), FNLROW(NCDIM)
      DOUBLE PRECISION DIAGJ(NCDIM), SPDIAG(NCDIM), SOL(NCDIM), Z(N)
      LOGICAL LMAX
      INTEGER I, IAB, IBACK, II, ISAVE, ITRY,
     1        NCUTM1, NEAR
      DOUBLE PRECISION BIG, DEL, DELMAX, DELMIN, DFNEW, DFUN,
     1        FDIF, FGAM, FNORM, FNRMNW, FOMX, FRAT, FZ, GAMMA,
     2        GAMNEW, PROD, RATGAM, REGINV, SINGTL, TSTVAL,
     3        VAL, VOL, XMULT, YDI, YDMIN
CMM   INTEGER IABS
      DATA SINGTL/ 1.0D-4/
      DATA BIG/ 1.0D+10/  
      SAVE
*      write(*,*) 'SOL',sol
*      write(*,*) 'FNLIN',fnlin
      NFCNT=0
      FDIF=FMAJOR-FMINOR
      PROD=1.0D+0
      DO 10 I=1,N
C ** JMB
      IF (ABS(DELPLS(I)).LT.1.D-16) DELPLS(I)=0.D0
      IF (ABS(DELNEG(I)).LT.1.D-16) DELNEG(I)=0.D0
C -- JMB
      PROD=PROD*(DELPLS(I)+DELNEG(I))
 10   CONTINUE
C ** JMB
      IF (VOL.NE.0.D0) THEN
        GAMMA=PROD/VOL
        FGAM=GAMMA*FMAJOR+(1.0D+0-GAMMA)*FMINOR
      ELSE
        GAMMA=1.0D0
        FGAM = 0.0D0
      ENDIF
C -- JMB
      DO 20 I=1,N
      Z(I)=X(I)
 20   CONTINUE
      DO 50 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.LT.0) GOTO 30
      Z(IAB)=X(IAB)+DELPLS(IAB)
      GOTO 40
 30   Z(IAB)=X(IAB)-DELNEG(IAB)
 40   FNEW(I)=DFUN(N,Z)
      Z(IAB)=X(IAB)
      NFCNT=NFCNT+1
 50   CONTINUE
   60 YDMIN= 1.0D+30
      NEAR=0
      DO 70 I=1,NCUT
      YDI=ABS(FMAJOR-FNEW(I))
      IF(YDI.GT.YDMIN) GOTO 70
      YDMIN=YDI
      NEAR=I
 70   CONTINUE
      IF(NEAR.EQ.0) RETURN
      IF((LMAX.AND.FNEW(NEAR).LT.FGAM).OR.
     1(.NOT.LMAX.AND.FNEW(NEAR).GT.FGAM)) GOTO 130
      ISAVE=ICUT(NEAR)
      IF(NEAR.EQ.NCUT.OR.NCUT.EQ.1) GOTO 90
      NCUTM1=NCUT-1
      DO 80 I=NEAR,NCUTM1
      ICUT(I)=ICUT(I+1)
      FNEW(I)=FNEW(I+1)
 80   CONTINUE
 90   NCUT=NCUT-1
      IF(ISAVE.LT.0) GOTO 100
      DELPLS(ISAVE)=XUP(ISAVE)-X(ISAVE)
*      write(*,*) 'A:delpls(isave),isave,xup(isave),x(isave)'
*     &     ,delpls(isave),isave,xup(isave),x(isave)
      GOTO 110
 100  IAB=ABS(ISAVE)
      DELNEG(IAB)=X(IAB)-XLOW(IAB)
*      write(*,*) 'delneg(iab),iab,x(iab),xlow(iab)',delneg(iab),iab
*     &     ,x(iab),xlow(iab)
 110  PROD=1.0D+0
      DO 120 I=1,N
C ** JMB
      IF (ABS(DELPLS(I)).LT.1.D-16) DELPLS(I)=0.D0
      IF (ABS(DELNEG(I)).LT.1.D-16) DELNEG(I)=0.D0
C -- JMB
      PROD=PROD*(DELPLS(I)+DELNEG(I))
 120  CONTINUE
C ** JMB
      IF (VOL.NE.0.D0) THEN
        GAMMA=PROD/VOL
        FGAM=GAMMA*FMAJOR+(1.0D+0-GAMMA)*FMINOR
      ELSE
        GAMMA=1.0D0
        FGAM = 0.0D0
      ENDIF
C -- JMB
      IF(NCUT.EQ.0) RETURN
      GOTO 60
 130  DO 160 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.LT.0) GOTO 140
      DEL=DELPLS(IAB)
      GOTO 150
 140  DEL=DELNEG(IAB)
 150  FORIG(I)=FNEW(I)
      IF (DEL.NE.0.D0) THEN
         DFORIG(I)=(FNEW(I)-FMAJOR)/DEL
      ELSE
         DFORIG(I)=0.D0         
      ENDIF
 160  CONTINUE
      CALL FEQN(NCUT,FORIG,FGAM,FNLIN)
      CALL RLEN(NCUT,FNLIN,FNORM)
 170  DO 180 I=1,NCUT
      FNLIN(I)=-FNLIN(I)
 180  CONTINUE
      IF(NCUT.EQ.1) GOTO 200
      DIAGJ(1)=DFORIG(1)
      SPDIAG(1)=-DFORIG(2)
      NCUTM1=NCUT-1
      DO 190 I=1,NCUTM1
      DIAGJ(I)=DFORIG(I)
      SPDIAG(I)=-DFORIG(I+1)
 190  CONTINUE
 200  DO 210 I=1,NCUT
      II=ICUT(I)
      II=ABS(II)
      IF (GAMMA.NE.0) THEN
        IF (ABS(DELPLS(II)+DELNEG(II)).NE.0.D0) THEN
           FNLROW(I)=-GAMMA*FDIF/(DELPLS(II)+DELNEG(II))
        ELSE
*           write(*,440) 'DELPLS(II)+DELNEG(II)='
*           write(*,*) DELPLS(II)+DELNEG(II)
           FNLROW(I)=0
        ENDIF
      ELSE
        FNLROW(I)=0.0
      ENDIF
 210  CONTINUE
      FNLROW(1)=DFORIG(1)+FNLROW(1)
      IF(NCUT.EQ.1) GOTO 230
      DO 220 I=1,NCUTM1
      XMULT=0.0D+0
      IF (ABS(FNLROW(I)).LT.BIG*ABS(DIAGJ(I))) THEN
         IF (DIAGJ(I).NE.0.D0) THEN
            XMULT=FNLROW(I)/DIAGJ(I)
         ELSE
*            write(*,440) 'diag(j)=0!'
            XMULT=0.D0
         ENDIF
      ENDIF
      FNLROW(I+1)=FNLROW(I+1)-XMULT*SPDIAG(I)
      FNLIN(NCUT)=FNLIN(NCUT)-XMULT*FNLIN(I)
 220  CONTINUE
 230  SOL(NCUT)=FNLIN(NCUT)
      IF (ABS(FNLROW(NCUT)).LT.BIG*ABS(FNLIN(NCUT))) THEN
         IF (FNLROW(NCUT).NE.0.D0) THEN
            SOL(NCUT)=FNLIN(NCUT)/FNLROW(NCUT)
         ELSE
            SOL(NCUT)=0.D0
         ENDIF
      ENDIF
*      write(*,*) 'FNLIN(NCUT),ncut,fnlrow(ncut)',FNLIN(NCUT),ncut
*     &     ,fnlrow(ncut)
      IF(NCUT.EQ.1) GOTO 250
      DO 240 I=2,NCUT
      IBACK=NCUT-I+1
      VAL=FNLIN(IBACK)-SOL(IBACK+1)*SPDIAG(IBACK)
      SOL(IBACK)=VAL
*     write(*,*) 'sol(iback),iback,val',sol(iback),iback,val
      IF(ABS(DIAGJ(IBACK)).LT.BIG*ABS(VAL)) THEN 
         IF (DIAGJ(IBACK).NE.0.D0) THEN
            SOL(IBACK)=VAL/DIAGJ(IBACK)
         ELSE
            SOL(IBACK)=0.D0
         ENDIF
      ENDIF
 240  CONTINUE
 250  ITRY=0
      DO 280 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.LT.0) GOTO 260
      DELMAX=FRACT*(XUP(IAB)-X(IAB)-DELPLS(IAB))
      DELMIN=-DELPLS(IAB)
      GOTO 270
 260  DELMAX=FRACT*(X(IAB)-XLOW(IAB)-DELNEG(IAB))
      DELMIN=-DELNEG(IAB)
 270  IF(SOL(I).GT.DELMAX) SOL(I)=0.75D+0*DELMAX
      IF(SOL(I).LT.DELMIN) SOL(I)=0.75D+0*DELMIN
 280  CONTINUE
 290  DO 310 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.LT.0) GOTO 300
      DELPLS(IAB)=DELPLS(IAB)+SOL(I)
*      write(*,*) 'A:delpls(iab),iab,sol(i),i,delmax,delmin,',delpls(iab)
*     &     ,iab,sol(i),i,delmax,delmin
      GOTO 310
 300  DELNEG(IAB)=DELNEG(IAB)+SOL(I)
*      write(*,*) 'A:delneg(iab),iab,sol(i),i',delneg(iab),iab
*     &    ,sol(i),i
 310  CONTINUE
      DO 340 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.LT.0) GOTO 320
      Z(IAB)=X(IAB)+DELPLS(IAB)
      GOTO 330
 320  Z(IAB)=X(IAB)-DELNEG(IAB)
 330  CONTINUE
*      write(*,*) 'iab,x(iab),z(iab),delpls,delneg',iab,x(iab),z(iab)
*     &     ,delpls,delneg
      FZ=DFUN(N,Z)
      NFCNT=NFCNT+1
      FNEW(I)=FZ
      Z(IAB)=X(IAB)
 340  CONTINUE
      PROD=1.0D+0
      DO 350 I=1,N
C ** JMB
      IF (ABS(DELPLS(I)).LT.1.D-16) DELPLS(I)=0.D0
      IF (ABS(DELNEG(I)).LT.1.D-16) DELNEG(I)=0.D0
C -- JMB
      PROD=PROD*(DELPLS(I)+DELNEG(I))
 350  CONTINUE
      IF (VOL.NE.0.D0) THEN
         GAMNEW=PROD/VOL
      ELSE
         GAMNEW=0.D0
      ENDIF
      FGAM=GAMNEW*FMAJOR+(1.0D+0-GAMNEW)*FMINOR
      CALL FEQN(NCUT,FNEW,FGAM,FNLIN)
      CALL RLEN(NCUT,FNLIN,FNRMNW)
      IF(FNRMNW.GT.FNORM) GOTO 380
      FOMX=0.0D+0
      DO 360 I=1,NCUT
      IF(ABS(FNEW(I)).GT.FOMX) FOMX=ABS(FNEW(I))
      DFNEW=FNEW(I)-FORIG(I)
      TSTVAL=1.0D+0
      IF (ABS(SOL(I)).LT.BIG*ABS(DFNEW)) THEN
         IF (SOL(I).NE.0.D0) THEN
            TSTVAL=DFNEW/SOL(I)
         ELSE
            TSTVAL=0.D0
         ENDIF
      ENDIF
      IF(ABS(TSTVAL).LT.SINGTL*ABS(DFORIG(I))) TSTVAL=SINGTL*DFORIG(I
     1)
      DFORIG(I)=TSTVAL
      FORIG(I)=FNEW(I)
 360  CONTINUE
      FNORM=FNRMNW
      FOMX=MAX(FOMX,ABS(FGAM))
      FRAT=FNORM/(1.0D+0+FOMX)
      REGINV=1.0D+0/REGTOL
      RATGAM=GAMNEW/GAMMA
      GAMMA=GAMNEW
C---  Activate to do debugging
C     WRITE(6,420) FNORM,FOMX,FRAT
C     WRITE(6,430) FGAM,RATGAM
      FGAM=GAMMA*FMAJOR+(1.0D+0-GAMMA)*FMINOR
      IF(RATGAM.GT.REGTOL.AND.RATGAM.LT.REGINV) GOTO 370
      IF(FRAT.LT.FTOL) GOTO 370
      GOTO 170
 370  RETURN
 380  ITRY=ITRY+1
      IF(ITRY.GT.2) RETURN
      DO 410 I=1,NCUT
      II=ICUT(I)
      IAB=ABS(II)
      IF(II.GE.0) THEN
         DELPLS(IAB)=DELPLS(IAB)-SOL(I)
*         write(*,*) 'B:delpls(iab),iab,sol(i),i',delpls(iab),iab
*     &        ,sol(i),i
      ELSE
         DELNEG(IAB)=DELNEG(IAB)-SOL(I)
*         write(*,*) 'B:delneg(iab),iab,sol(i),i',delneg(iab),iab
*     &        ,sol(i),i
      ENDIF
      SOL(I)=SOL(I)*0.25D+0
 410  CONTINUE
      GOTO 290
 420  FORMAT(' FNORM, FOMX, FRAT', 3(1PD15.5))
 430  FORMAT(' FGAM, RATGAM', 2(1PD15.5))
 440  FORMAT(A)
      END


