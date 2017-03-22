*CMZ :          23/08/93  13.30.10  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE TSTEXT (N,XEXTR,XLOW,XUP,EDGFAC,FRACT,NCUT,NCDIM,ICUT,D
     1ELPLS,DELNEG)
      INTEGER N, NCDIM, NCUT
      INTEGER ICUT(NCDIM)
      DOUBLE PRECISION EDGFAC, FRACT
      DOUBLE PRECISION XEXTR(N), XLOW(N), XUP(N), DELPLS(N), DELNEG(N)
      INTEGER I
      DOUBLE PRECISION DIFX
      SAVE
      NCUT=0
      DO 10 I=1,N
      DIFX=XUP(I)-XEXTR(I)
      DELPLS(I)=DIFX
*      write(*,*) 'delpls(i),i,xup(i),xextr(i)',delpls(i),i,xup(i)
*     &     ,xextr(i)
      IF(DIFX.LT.EDGFAC*(XUP(I)-XLOW(I))) GOTO 10
      DELPLS(I)=FRACT*DIFX
*      write(*,*) 'delpls(i),i,fract',delpls(i),i,fract
      NCUT=NCUT+1
      ICUT(NCUT)=I
 10   CONTINUE
      DO 20 I=1,N
      DIFX=XEXTR(I)-XLOW(I)
      DELNEG(I)=DIFX
*      write(*,*) 'delneg(i),i,xup(i),xextr(i)',delneg(i),i,xup(i)
*     &     ,xextr(i)
      IF(DIFX.LT.EDGFAC*(XUP(I)-XLOW(I))) GOTO 20
      DELNEG(I)=FRACT*DIFX
*      write(*,*) 'delpls(i),i,fract',delpls(i),i,fract
      NCUT=NCUT+1
      ICUT(NCUT)=-I
 20   CONTINUE
      RETURN
      END








