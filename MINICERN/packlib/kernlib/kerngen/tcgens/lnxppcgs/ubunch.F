*
* $Id$
*
* $Log$
* Revision 1.1  1998/09/25 09:32:14  mclareni
* Modifications for the Mklinux port flagged by CERNLIB_PPC
*
* Revision 1.1.1.1  1996/02/15 17:50:20  mclareni
* Kernlib
*
*
      SUBROUTINE UBUNCH (MS,MT,NCHP)
C
C CERN PROGLIB# M409    UBUNCH          .VERSION KERNLNX  1.02  940511
C ORIG. 03/02/89 K.M.STORR
C

      DIMENSION    MS(99), MT(99), NCHP(9)
      data iblan1/x'20202020'/
      data mask1/x'ff000000'/

      NCH = NCHP(1)
      IF   (NCH)             91,39,11
   11 NWT    = ishftr (NCH,2)
      NTRAIL = AND (NCH,3)
      JS     = 0
      IF (NWT.EQ.0)          GO TO 31

C--                Pack the initial complete words

      DO 24  JT=1,NWT
      MT(JT) = OR (OR (OR (
     +                  AND(MS(JS+1),MASK1),
     +          ISHFTR (AND(MS(JS+2),MASK1), 8)),
     +          ISHFTR (IAND(MS(JS+3),MASK1),16)),
     +          ISHFTR      (MS(JS+4),       24) )
   24 JS = JS + 4

      IF (NTRAIL.EQ.0)       RETURN

C--                Pack the trailing word

   31 MWD = IBLAN1
      JS  = NCH

      DO 34 JT=1,NTRAIL
      MWD = OR (ISHFTR(MWD,8), AND(MS(JS),MASK1))
   34 JS  = JS - 1
      MT(NWT+1) = MWD
   39 RETURN

   91 CALL ABEND
      END
