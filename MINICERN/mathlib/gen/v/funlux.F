      SUBROUTINE FUNLXP (FUNC,XFCUM,X2LOW,X2HIGH)
C         F. JAMES,   Sept, 1994
C
C         Prepares the user function FUNC for FUNLUX
C         Inspired by and mostly copied from FUNPRE and FUNRAN
C         except that 
C    1. FUNLUX uses RANLUX underneath, 
C    2. FUNLXP expands the first and last bins to cater for
C              functions with long tails on left and/or right,
C    3. FUNLXP calls FUNPCT to do the actual finding of percentiles.
C    4. both FUNLXP and FUNPCT use RADAPT for Gaussian integration.
C
      EXTERNAL FUNC
      COMMON/FUNINT/TFTOT
      DIMENSION XFCUM(200)
      PARAMETER (RTEPS=0.0002) 
      SAVE IFUNC
      DATA IFUNC/0/
      IFUNC = IFUNC + 1
C         FIND RANGE WHERE FUNCTION IS NON-ZERO.
      CALL FUNLZ(FUNC,X2LOW,X2HIGH,XLOW,XHIGH)
      XRANGE = XHIGH-XLOW
      IF(XRANGE .LE. 0.)  THEN
        WRITE (6,'(A,2G15.5)') ' FUNLXP finds function range .LE.0',
     +  XLOW,XHIGH
        GO TO 900
      ENDIF
      CALL RADAPT(FUNC,XLOW,XHIGH,1,RTEPS,0.,TFTOT ,UNCERT)
      WRITE(6,1003) IFUNC,XLOW,XHIGH,TFTOT
 1003 FORMAT(' FUNLXP: integral of USER FUNCTION',
     +  I3,' from ',E12.5,' to ',E12.5,' is ',E14.6)
C
C      WRITE (6,'(A,A)') ' FUNLXP preparing ',
C     + 'first the whole range, then left tail, then right tail.'
      CALL FUNPCT(FUNC,IFUNC,XLOW,XHIGH,XFCUM,1,99,TFTOT,IERR)
      IF (IERR .GT. 0)  GO TO 900
      X2 = XFCUM(3)
      CALL RADAPT(FUNC,XLOW,X2,1,RTEPS,0.,TFTOT1 ,UNCERT)
      CALL FUNPCT(FUNC,IFUNC,XLOW,X2 ,XFCUM,101,49,TFTOT1,IERR)
      IF (IERR .GT. 0)  GO TO 900
      X3 = XFCUM(98)
      CALL RADAPT(FUNC,X3,XHIGH,1,RTEPS,0.,TFTOT2 ,UNCERT)
      CALL FUNPCT(FUNC,IFUNC,X3,XHIGH,XFCUM,151,49,TFTOT2,IERR)
      IF (IERR .GT. 0)  GO TO 900
      WRITE(6,1001) IFUNC,XLOW,XHIGH
 1001 FORMAT(' FUNLXP has prepared USER FUNCTION',I3,
     + ' between',G12.3,' and',G12.3,' for FUNLUX')
      RETURN
  900 CONTINUE
      WRITE (6,'(A)') ' Fatal error in FUNLXP. FUNLUX will not work.'
      END
C
      SUBROUTINE FUNPCT(FUNC,IFUNC,XLOW,XHIGH,XFCUM,NLO,NBINS,TFTOT,
     +                  IERR)
C        Array XFCUM is filled from NLO to NLO+NBINS, which makes
C        the number of values NBINS+1, or the number of bins NBINS
      EXTERNAL FUNC
      DIMENSION XFCUM(*) 
      PARAMETER (RTEPS=0.005, NZ=10, MAXZ=20, NITMAX=6,PRECIS=1.E-6)
CC      DOUBLE PRECISION TPCTIL, TZ, TCUM, XINCR, DTABS,
CC     +  TINCR, TZMAX, XBEST, DTBEST, DTPAR2
C
      IERR = 0
      IF (TFTOT .LE. 0.) GO TO 900
      TPCTIL = TFTOT/NBINS
      TZ = TPCTIL/NZ
      TZMAX = TZ * 2.
      XFCUM(NLO) = XLOW
      XFCUM(NLO+NBINS) = XHIGH
      X = XLOW
      F = FUNC(X)
      IF (F .LT. 0.) GO TO 900
C         Loop over percentile bins
      DO 600 IBIN = NLO, NLO+NBINS-2
      TCUM = 0.
      X1 = X
      F1 = F
      DXMAX = (XHIGH -X) / NZ
      FMIN = TZ/DXMAX
      FMINZ = FMIN
C         Loop over trapezoids within a supposed percentile
      DO 500 IZ= 1, MAXZ
      XINCR = TZ/MAX(F1,FMIN,FMINZ)
  350 X = X1 + XINCR
      F = FUNC(X)
      IF (F .LT. 0.) GO TO 900
      TINCR = (X-X1) * 0.5 * (F+F1)
      IF (TINCR .LT. TZMAX) GO TO 370
      XINCR = XINCR * 0.5
      GO TO 350
  370 CONTINUE
      TCUM = TCUM + TINCR
      IF (TCUM .GE. TPCTIL*0.99) GO TO 520
      FMINZ = TZ*F/ (TPCTIL-TCUM)
      F1 = F
      X1 = X
  500 CONTINUE
      WRITE(6,'(A)') ' FUNLUX:  WARNING. FUNPCT fails trapezoid.'
C         END OF TRAPEZOID LOOP
C         Adjust interval using Gaussian integration with 
C             Newton corrections since F is the derivative 
  520 CONTINUE
      X1 = XFCUM(IBIN)
      XBEST = X
      DTBEST = TPCTIL
      TPART = TPCTIL
C         Allow for maximum NITMAX more iterations on RADAPT
      DO 550 IHOME= 1, NITMAX
  535 XINCR = (TPCTIL-TPART) / MAX(F,FMIN)
      X = XBEST + XINCR
      X2 = X
        IF (IHOME .GT. 1 .AND. X2 .EQ. XBEST) THEN
        WRITE (6,'(A,G12.3)')
     +  ' FUNLUX: WARNING from FUNPCT: insufficient precision at X=',X
        GO TO 580
        ENDIF
      REFX = ABS(X)+PRECIS
      CALL RADAPT(FUNC,X1,X2,1,RTEPS,0.,TPART2,UNCERT)
      DTPAR2 = TPART2-TPCTIL
      DTABS = ABS(DTPAR2)
      IF(ABS(XINCR)/REFX .LT. PRECIS) GOTO 545
      IF(DTABS .LT. DTBEST) GOTO 545
      XINCR = XINCR * 0.5
      GOTO 535
  545 DTBEST = DTABS
      XBEST = X
      TPART = TPART2
      F = FUNC(X)
      IF(F .LT. 0.) GOTO 900
      IF(DTABS .LT. RTEPS*TPCTIL) GOTO 580
  550 CONTINUE
      WRITE (6,'(A,I4)') 
     +   ' FUNLUX: WARNING from FUNPCT: cannot converge, bin',IBIN
C
  580 CONTINUE
      XINCR = (TPCTIL-TPART) / MAX(F,FMIN)
      X = XBEST + XINCR
      XFCUM(IBIN+1) = X
      F = FUNC(X)
      IF(F .LT. 0.) GOTO 900
  600 CONTINUE
C         END OF LOOP OVER BINS
      X1 = XFCUM(NLO+NBINS-1)
      X2 = XHIGH
      CALL RADAPT(FUNC,X1,X2,1,RTEPS,0.,TPART ,UNCERT)
      ABERR = ABS(TPART-TPCTIL)/TFTOT
C      WRITE(6,1001) IFUNC,XLOW,XHIGH
      IF(ABERR .GT. RTEPS)  WRITE(6,1002) ABERR
      RETURN
  900 WRITE(6,1000) X,F
      IERR = 1
      RETURN
 1000 FORMAT(/' FUNLUX fatal error in FUNPCT: function negative:'/
     + ,' at X=',E15.6,', F=',E15.6/)
C 1001 FORMAT(' FUNPCT has prepared USER FUNCTION',I3,
C     + ' between',G12.3,' and',G12.3,' for FUNLUX.')
 1002 FORMAT(' WARNING: Relative error in cumulative distribution',
     + ' may be as big as',F10.7)
      END
      SUBROUTINE FUNLUX(ARRAY,XRAN,LEN)
C         Generation of LEN random numbers in any given distribution,
C         by 4-point interpolation in the inverse cumulative distr.
C         which was previously generated by FUNLXP
      COMMON/FUNINT/XUNI
      DIMENSION ARRAY(200)
      DIMENSION XRAN(LEN)
C  Bin width for main sequence, and its inverse
      PARAMETER (GAP= 1./99.,  GAPINV=99.)
C  Top of left tail, bottom of right tail (each tail replaces 2 bins)
      PARAMETER (TLEFT= 2./99.,BRIGHT=97./99.)
C  Bin width for minor sequences (tails), and its inverse
      PARAMETER (GAPS=TLEFT/49.,  GAPINS=1./GAPS)
C
C   The array ARRAY is assumed to have the following structure: 
C        ARRAY(1-100) contains the 99 bins of the inverse cumulative
C                     distribution of the entire function.
C        ARRAY(101-150) contains the 49-bin blowup of main bins 
C                       1 and 2 (left tail of distribution)
C        ARRAY(151-200) contains the 49-bin blowup of main bins 
C                       98 and 99 (right tail of distribution)
C
      CALL RANLUX(XRAN,LEN)
      
      DO 500 IBUF= 1, LEN
      X = XRAN(IBUF)
      J = INT(  X    *GAPINV) + 1
      IF (J .LT. 3)  THEN
         J1 = INT( X *GAPINS)
             J = J1 + 101
             J = MAX(J,102)
             J = MIN(J,148)
         P = (   X -GAPS*(J1-1)) * GAPINS
         A = (P+1.0) * ARRAY(J+2) - (P-2.0)*ARRAY(J-1)
         B = (P-1.0) * ARRAY(J) - P * ARRAY(J+1)
         XRAN(IBUF) = A*P*(P-1.0)*0.16666667  + B*(P+1.0)*(P-2.0)*0.5
      ELSE IF (J .GT. 97)  THEN
         J1 = INT((X-BRIGHT)*GAPINS) 
             J = J1 + 151
             J = MAX(J,152)
             J = MIN(J,198)
         P = (X -BRIGHT -GAPS*(J1-1)) * GAPINS
         A = (P+1.0) * ARRAY(J+2) - (P-2.0)*ARRAY(J-1)
         B = (P-1.0) * ARRAY(J) - P * ARRAY(J+1)
         XRAN(IBUF) = A*P*(P-1.0)*0.16666667  + B*(P+1.0)*(P-2.0)*0.5
      ELSE
CC      J = MAX(J,2)
CC      J = MIN(J,98)
         P = (   X -GAP*(J-1)) * GAPINV
         A = (P+1.0) * ARRAY(J+2) - (P-2.0)*ARRAY(J-1)
         B = (P-1.0) * ARRAY(J) - P * ARRAY(J+1)
         XRAN(IBUF) = A*P*(P-1.0)*0.16666667  + B*(P+1.0)*(P-2.0)*0.5
      ENDIF
  500 CONTINUE
      XUNI = X
      RETURN
      END
      SUBROUTINE FUNLZ(FUNC,X2LOW,X2HIGH,XLOW,XHIGH)
C         FIND RANGE WHERE FUNC IS NON-ZERO.
C         WRITTEN 1980, F. JAMES
C         MODIFIED, NOV. 1985, TO FIX BUG AND GENERALIZE
C         TO FIND SIMPLY-CONNECTED NON-ZERO REGION (XLOW,XHIGH)
C         ANYWHERE WITHIN THE GIVEN REGION (X2LOW,H2HIGH).
C            WHERE 'ANYWHERE' MEANS EITHER AT THE LOWER OR UPPER
C            EDGE OF THE GIVEN REGION, OR, IF IN THE MIDDLE,
C            COVERING AT LEAST 1% OF THE GIVEN REGION.
C         OTHERWISE IT IS NOT GUARANTEED TO FIND THE NON-ZERO REGION.
C         IF FUNCTION EVERYWHERE ZERO, FUNLZ SETS XLOW=XHIGH=0.
      EXTERNAL FUNC
      XLOW = X2LOW
      XHIGH = X2HIGH
C         FIND OUT IF FUNCTION IS ZERO AT ONE END OR BOTH
      XMID = XLOW
      IF (FUNC(XLOW) .GT. 0.) GO TO 120
      XMID = XHIGH
      IF (FUNC(XHIGH) .GT. 0.)  GO TO 50
C         FUNCTION IS ZERO AT BOTH ENDS,
C         LOOK FOR PLACE WHERE IT IS NON-ZERO.
      DO 30 LOGN= 1, 7
      NSLICE = 2**LOGN
      DO 20 I= 1, NSLICE, 2
      XMID = XLOW + I * (XHIGH-XLOW) / NSLICE
      IF (FUNC(XMID) .GT. 0.)  GO TO 50
   20 CONTINUE
   30 CONTINUE
C         FALLING THROUGH LOOP MEANS CANNOT FIND NON-ZERO VALUE
      WRITE(6,554)
      WRITE(6,555) XLOW, XHIGH
      XLOW = 0.
      XHIGH = 0.
      GO TO 220
C
   50 CONTINUE
C         DELETE 'LEADING' ZERO RANGE
      XH = XMID
      XL = XLOW
      DO 70 K= 1, 20
      XNEW = 0.5*(XH+XL)
      IF (FUNC(XNEW) .EQ. 0.) GO TO 68
      XH = XNEW
      GO TO 70
   68 XL = XNEW
   70 CONTINUE
      XLOW = XL
      WRITE(6,555) X2LOW,XLOW
  120 CONTINUE
      IF (FUNC(XHIGH) .GT. 0.) GO TO 220
C         DELETE 'TRAILING' RANGE OF ZEROES
      XL = XMID
      XH = XHIGH
      DO 170 K= 1, 20
      XNEW = 0.5*(XH+XL)
      IF (FUNC(XNEW) .EQ. 0.) GO TO 168
      XL = XNEW
      GO TO 170
  168 XH = XNEW
  170 CONTINUE
      XHIGH = XH
      WRITE(6,555) XHIGH, X2HIGH
C
  220 CONTINUE
      RETURN
  554 FORMAT('0CANNOT FIND NON-ZERO FUNCTION VALUE')
  555 FORMAT(' FUNCTION IS ZERO FROM X=',E12.5,' TO ',E12.5)
      END
