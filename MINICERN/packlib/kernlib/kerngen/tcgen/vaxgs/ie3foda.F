*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:59  mclareni
* Kernlib
*
*
      SUBROUTINE IE3FOD (MS,MT,NDPN,JBAD)
C
C CERN PROGLIB# M220    IE3FOD          .VERSION KERNVAX  2.44  931108
C     Alfred Nathaniel - 12-Oct-93
C-    Convert double precision for output with copy
C-    from source in native to target in IEEE data format

      INTEGER MS(2,*), MT(2,*), IDBL(2)
      DOUBLE PRECISION DBL
      EQUIVALENCE(DBL,IDBL(1))

      INCLUDE '($CVTDEF)'
      INTEGER CVT$CONVERT_FLOAT

      JBAD = 0
      DO 10 I=1,NDPN
         IDBL(1) = MS(1,I)
         IDBL(2) = MS(2,I)
         ISTAT = CVT$CONVERT_FLOAT(
     +        DBL, %VAL(CVT$K_VAX_G),
     +        DBL, %VAL(CVT$K_IEEE_T),
     +        %VAL(CVT$M_ROUND_TO_NEAREST) )
         IF(.NOT.ISTAT) THEN
            JBAD = I
         ENDIF
         MT(1,I) = IDBL(2)
         MT(2,I) = IDBL(1)
 10   CONTINUE
      END
