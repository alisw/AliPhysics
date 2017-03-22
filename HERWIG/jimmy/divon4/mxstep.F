*CMZ :          23/08/93  13.30.12  by  Jonathan Butterworth
*-- Author :
      SUBROUTINE MXSTEP (X,XLOW,XHI,P,STPPOS,STPNEG)
      DOUBLE PRECISION X, XLOW, XHI, P, STPPOS, STPNEG
      SAVE
      IF(P.EQ.0.0D+0) RETURN
      IF(P.GT.0.0D+0) GOTO 10
      STPPOS=-(X-XLOW)/P
      STPNEG=-(XHI-X)/P
      RETURN
 10   STPPOS=(XHI-X)/P
      STPNEG=(X-XLOW)/P
      RETURN
      END
