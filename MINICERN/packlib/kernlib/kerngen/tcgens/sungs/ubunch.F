*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:22  mclareni
* Kernlib
*
*
      SUBROUTINE UBUNCH (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBUNCH          .VERSION KERNSUN  1.06  920511
C ORIG. 22/09/88, JZ
C

      DIMENSION    MS(99), MT(99), NCHP(9)
C-  - PARAMETER    (IBLAN1 = X'20202020')
      PARAMETER    (IBLAN1 =  538976288 )
C-  - PARAMETER    (MASK1  = X'FF000000')
      PARAMETER    (MASK1  =  -16777216 )


      NCH = NCHP(1)
      IF   (NCH)             91,39,11
   11 NWT    = RSHIFT (NCH,2)
      NTRAIL = AND (NCH,3)
      JS     = 0
      IF (NWT.EQ.0)          GO TO 31

C--                Pack the initial complete words

      DO 24  JT=1,NWT
      MT(JT) = OR (OR (OR (
     +                   AND(MS(JS+1),MASK1),
#if defined(CERNLIB_BUGLRSHFT)
     +            ishft (AND(MS(JS+2),MASK1), -8)),
     +            ishft (AND(MS(JS+3),MASK1),-16)),
     +            ishft     (MS(JS+4),       -24) )
#endif
#if !defined(CERNLIB_BUGLRSHFT)
     +           lrshft (AND(MS(JS+2),MASK1), 8)),
     +           lrshft (AND(MS(JS+3),MASK1),16)),
     +           lrshft     (MS(JS+4),       24) )
#endif
   24 JS = JS + 4

      IF (NTRAIL.EQ.0)       RETURN

C--                Pack the trailing word

   31 MWD = IBLAN1
      JS  = NCH

      DO 34 JT=1,NTRAIL
#if defined(CERNLIB_BUGLRSHFT)
      MWD = OR (ishft(MWD,-8), AND(MS(JS),MASK1))
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MWD = OR (lrshft(MWD,8), AND(MS(JS),MASK1))
#endif
   34 JS  = JS - 1
      MT(NWT+1) = MWD
   39 RETURN

   91 CALL ABEND
      END
