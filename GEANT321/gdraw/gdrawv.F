*
* $Id$
*
* $Log$
* Revision 1.2  1996/09/30 13:37:24  ravndal
* Backward compatibility for view banks
*
* Revision 1.1.1.1  1995/10/24 10:20:27  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRAWV(U,V,NP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       GDRAWV stores sparse segments ( NP = 2 ) in a buffer.    *
C.    *       When the buffer is full it is scanned to collect         *
C.    *       segments into polylines ( NP > 2 ), that will be         *
C.    *       drawn by the GDRAWP routine.                             *
C.    *                                                                *
C.    *       Variables used in the routine :                          *
C.    *                                                                *
C.    *       EPS    = Error range of user coordinates until which     *
C.    *                they are considered as coincident points        *
C.    *                                                                *
C.    *       I0     = Pointer to starting point of polyline (in X,Y)  *
C.    *       I1     = Pointer to starting point of segment (in X,Y)   *
C.    *       I2     = Pointer to ending point of segment (in X,Y)     *
C.    *       INEW   = Number of points of polyline (in XNEW,YNEW)     *
C.    *                                                                *
C.    *       NSEG   = Current number of segments                      *
C.    *       MAXSEG = Max number of segments                          *
C.    *       MAXNEW = Max length of polyline                          *
C.    *                                                                *
C.    *       X      = Vector of x-coordinates of input segments       *
C.    *       Y      = Vector of y-coordinates of input segments       *
C.    *       IEMPTY = Vector to flag empty points; if IEMPTY(I)=1     *
C.    *                then X(I),Y(I) has been already linked          *
C.    *                                                                *
C.    *                (X,Y,IEMPTY are dimensioned to MAXSEG*2)        *
C.    *                                                                *
C.    *       IATTL  = Vector of attributes of input segments,         *
C.    *                dimensioned to MAXSEG                           *
C.    *                                                                *
C.    *       XNEW   = Vector of x-coordinates of output polyline      *
C.    *       YNEW   = Vector of y-coordinates of output polyline      *
C.    *                                                                *
C.    *                (XNEW,YNEW are dimensioned to MAXNEW)           *
C.    *                                                                *
C.    *       LUDEB  = If < 0 then no debug printout will be done      *
C.    *                                                                *
C.    *       ISTF   = Flag for statistics printout;                   *
C.    *                if ISTF=1 a statistics of polylines length      *
C.    *                will be done (to be printed interactively)      *
C.    *                                                                *
C.    *       LINBUF = Flag to determine if GDRAWV has to perform      *
C.    *                line buffering logic (LINBUF=1) or not          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDAHIT, GDARC, GDARR,       *
C.    *                   GDCIRC, GDCIRR, GDCIRY, GDCUT, GDHEAD,       *
C.    *                   GDLINE, GDMAN, GDPRTR, GDRAW, GDRECT,        *
C.    *                   GDSCAL, GDSHOW, GDSURF, GDTREE               *
C.    *       Author : P.Zanarini   S.Giani 1992  *********            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcspee.inc"
      CHARACTER*80 ILINE
      DIMENSION U(*),V(*)
      DIMENSION X(160),Y(160),IEMPTY(160)
      DIMENSION IATTL(80)
      DIMENSION XNEW(45),YNEW(45)
      DIMENSION IDBIST(2,20),IDBERR(10)
      SAVE X,Y,IEMPTY,IATTL,IDBIST,IDBERR,XNEW,YNEW
      SAVE MAXSEG,MAXNEW,EPS,LUDEB,ISTF,NSEG,INEW
      DATA MAXSEG/80/
      DATA MAXNEW/45/
      DATA EPS/0.00001/
      DATA LUDEB/-1/
      DATA ISTF/0/
      DATA IDBERR,IDBIST/50*0/
C.
C.    ------------------------------------------------------------------
C.
      IF (NP.EQ.-1) THEN
C
C             Initialize (open) the line buffer
C
         LINBUF=1
         NSEG=0
         INEW=0
         GO TO 999
C
      ELSE IF (NP.EQ.0) THEN
C
C             Empties (close) the line buffer
C
         LINBUF=0
         GO TO 10
C
      ELSE IF (NP.GT.2) THEN
C
C             Call directly GDRAWP
C
         CALL GDRAWP(U,V,NP)
         GO TO 999
C
C
      ELSE IF (LINBUF.EQ.0) THEN
C
C             Call directly GDRAWP
C
         CALL GDRAWP(U,V,NP)
         GO TO 999
C
      ENDIF
C
      LTY=IBITS(LINATT,10,3)
      CALL UCTOH('ON  ',IFLH,4,4)
      IF(LTY.NE.6.AND.IHIDEN.EQ.IFLH)THEN
         CALL GDRAWP(U,V,NP)
         GO TO 999
      ENDIF
C             Fill the line buffer
C
      NSEG=NSEG+1
      IATTL(NSEG)=LINATT
C
      NPOINT=NSEG*2
      X(NPOINT-1)=U(1)
      Y(NPOINT-1)=V(1)
      X(NPOINT)=U(2)
      Y(NPOINT)=V(2)
C
      IF (NSEG.NE.MAXSEG) GO TO 999
C
C             Line buffer full : check for joined lines
C
C             Initialisation phase
C
   10 CONTINUE
C
      I0=1
      INEW=0
      NPOINT=NSEG*2
C
      DO 20 I=1,NPOINT
         IEMPTY(I)=0
   20 CONTINUE
C
C             Loop phase
C
   30 CONTINUE
C
      NPOINT=NSEG*2
C
      DO 40 K=I0,NPOINT
         IF (IEMPTY(K).EQ.0) GO TO 50
   40 CONTINUE
C
C             Buffer empty
C
      NSEG=0
      GO TO 999
C
   50 CONTINUE
C
      KDIV2=(K+1)/2
      IF (INEW.EQ.0) LINATT=IATTL(KDIV2)
C
      I1=K
      INEW=INEW+1
      XNEW(INEW)=X(I1)
      YNEW(INEW)=Y(I1)
      IEMPTY(I1)=1
      I0=I1
C
   60 CONTINUE
C
      IF ((I1/2)*2.EQ.I1) THEN
         I2=I1-1
      ELSE
         I2=I1+1
      ENDIF
C
      I2DIV2=(I2+1)/2
      IF (INEW.EQ.0) LINATT=IATTL(I2DIV2)
C
      INEW=INEW+1
      XNEW(INEW)=X(I2)
      YNEW(INEW)=Y(I2)
      IEMPTY(I2)=1
C
      IF (INEW.EQ.MAXNEW) THEN
         IF (LUDEB.GE.0) THEN
            WRITE (CHMAIL,10000) NAMES(NLEVEL),INEW
            CALL GMAIL(0,0)
         ENDIF
         IDBERR(1)=IDBERR(1)+1
         GO TO 70
      ENDIF
C
C             Link is closed ?
C
      IF ((X(I2).LT.X(I0)-EPS).OR.(X(I2).GT.X(I0)+EPS)) GO TO 140
      IF ((Y(I2).LT.Y(I0)-EPS).OR.(Y(I2).GT.Y(I0)+EPS)) GO TO 140
C
C             Link closed
C
   70 CONTINUE
C
      IF (INEW.GT.3) THEN
         CALL GDRAWP(XNEW,YNEW,INEW)
      ELSE
         IF (INEW.EQ.3) CALL GDRAWP(XNEW,YNEW,INEW)
         IF (LUDEB.GE.0) THEN
            WRITE (CHMAIL,10100) NAMES(NLEVEL),INEW
            CALL GMAIL(0,0)
         ENDIF
         IDBERR(2)=IDBERR(2)+1
      ENDIF
C
   80 CONTINUE
C
C             Do some statistics if in statistics mode
C
      IF (ISTF.EQ.1) THEN
         IF (INEW.LE.3) GO TO 120
         LENPOL=INEW-1
   90    CONTINUE
         DO 100 I=1,20
            IF (LENPOL.EQ.IDBIST(1,I)) THEN
               IDBIST(2,I)=IDBIST(2,I)+1
               GO TO 120
            ENDIF
  100    CONTINUE
         DO 110 I=1,20
            IF (IDBIST(1,I).EQ.0) THEN
               IDBIST(1,I)=LENPOL
               GO TO 90
            ENDIF
  110    CONTINUE
         WRITE (CHMAIL,10800) LENPOL
         CALL GMAIL(0,0)
  120    CONTINUE
      ENDIF
C
      IF (LUDEB.LT.0) GO TO 130
C
C             Do some printout if in debugging mode
C
      LENPOL=INEW-1
      WRITE (CHMAIL,10400) NAMES(NLEVEL),LENPOL
      CALL GMAIL(0,0)
      IC0=I0/2
      IF ((I0/2)*2.NE.I0) IC0=IC0+1
      IC1=I1/2
      IF ((I1/2)*2.NE.I1) IC1=IC1+1
      ILINE=' '
      ILINE(IC1:IC1)='1'
      ILINE(IC0:IC0)='0'
      WRITE (CHMAIL,10500) I0,I1,I2,INEW,NSEG
      CALL GMAIL(1,0)
      WRITE (CHMAIL,10600) ILINE(1:NSEG)
      CALL GMAIL(0,0)
      NPOINT=NSEG*2
      NPOIN1=NPOINT-1
      WRITE (CHMAIL,10700) (IEMPTY(JJ),JJ=1,NPOIN1,2)
      CALL GMAIL(0,0)
      WRITE (CHMAIL,10700) (IEMPTY(JJ),JJ=2,NPOINT,2)
      CALL GMAIL(0,0)
C
  130 CONTINUE
      INEW=0
      GO TO 30
C
  140 CONTINUE
C
      NPOINT=NSEG*2
C
C             Compute starting point of next segment
C
      NEXT=0
C
      DO 150 K=I0,NPOINT
         IF (K.EQ.I2) GO TO 150
         IF ((X(I2).LT.X(K)-EPS).OR.(X(I2).GT.X(K)+EPS)) GO TO 150
         IF ((Y(I2).LT.Y(K)-EPS).OR.(Y(I2).GT.Y(K)+EPS)) GO TO 150
         NEXT=K
         GO TO 160
  150 CONTINUE
C
  160 CONTINUE
C
      IF (NEXT.NE.0) THEN
         I1=NEXT
         IEMPTY(I1)=1
         GO TO 60
      ENDIF
C
      IF (I0.EQ.1) THEN
         CALL GDRAWP(XNEW,YNEW,INEW)
         IF (LUDEB.GE.0) THEN
            WRITE (CHMAIL,10200) NAMES(NLEVEL),INEW
            CALL GMAIL(0,0)
         ENDIF
         IDBERR(3)=IDBERR(3)+1
         GO TO 80
      ENDIF
C
C             Reset flag : next call to GDRAWV will re-process
C             the points now forming an open polyline, trying then
C             (with the new points stored) to form a closed polyline
C
      DO 170 I=I0,NPOINT
         IEMPTY(I)=0
  170 CONTINUE
C
      NSEGOL=NSEG
C
C             Try to see if the buffers can be shifted
C
  180 CONTINUE
C
      DO 190 NLIN=1,NSEG
         IF (IEMPTY(NLIN*2).EQ.0) GO TO 210
  190 CONTINUE
C
      IF (NSEG.NE.NSEGOL) THEN
C
C             The buffers have been shifted enough
C
         GO TO 230
C
      ELSE
C
C             The buffers cannot be shifted at all
C
         IF (LUDEB.GE.0) THEN
            WRITE (CHMAIL,10300) NAMES(NLEVEL)
            CALL GMAIL(0,0)
         ENDIF
         IDBERR(4)=IDBERR(4)+1
         LINATT=IATTL(1)
         DO 200 I=1,NSEG
            J=I*2-1
            CALL GDRAWP(X(J),Y(J),2)
  200    CONTINUE
         INEW=0
         GO TO 30
C
      ENDIF
C
  210 CONTINUE
C
C             Shift the buffers
C
      IF (NLIN.EQ.1) GO TO 230
      NSEG=NSEG-NLIN+1
      N=NSEG*2
      J=NLIN*2-2
      DO 220 I=1,N
         X(I) = X(J+I)
         Y(I) = Y(J+I)
         IEMPTY(I) = IEMPTY(J+I)
  220 CONTINUE
      CALL UCOPY(IATTL(NLIN),IATTL(1),NSEG)
C
      IF (NSEG.GT.0) GO TO 180
C
      GO TO 999
C
  230 CONTINUE
C
C             If last call to GDRAWV (NP=0, i.e. close the buffer)
C             and some lines have still to be processed (NSEG>0)
C             then don't RETURN but re-execute the body of GDRAWV
C
      IF (NP.EQ.0.AND.NSEG.GT.0) GO TO 10
C
10000 FORMAT (' *** GDRAWV: ERROR 1 -  NAME = ',A4,' INEW = ',I3)
10100 FORMAT (' *** GDRAWV: ERROR 2 -  NAME = ',A4,' INEW = ',I3)
10200 FORMAT (' *** GDRAWV: ERROR 3 -  NAME = ',A4,' INEW = ',I3)
10300 FORMAT (' *** GDRAWV: ERROR 4 -  NAME = ',A4)
10400 FORMAT (' *** GDRAWV: GDRAWP called for ',A4,' (',I2,' segments)')
10500 FORMAT (' LOOP : I0=',I3,' - I1=',I3,' - I2=',I3,
     +        ' - INEW=',I3,' - NSEG=',I3)
10600 FORMAT (1X,A)
10700 FORMAT (1X,80I1)
10800 FORMAT (' *** STATISTICS OVERFLOW ***   LENPOL =',I3)
  999 END
